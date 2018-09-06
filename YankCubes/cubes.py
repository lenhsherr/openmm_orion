# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import traceback
from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube

from cuberecord.ports import RecordInputPort, RecordOutputPort

from datarecord import (Types,
                        Meta,
                        OEField,
                        OERecord,
                        OEFieldMeta)

import MDCubes.OpenMMCubes.utils as omm_utils

from openeye import oechem

from simtk.openmm import (app,
                          unit,
                          XmlSerializer)

from tempfile import TemporaryDirectory

import tarfile

import os

import itertools

from YankCubes.yank_templates import (yank_solvation_template,
                                      yank_binding_template)

from YankCubes import utils as yankutils

from yank.analyze import ExperimentAnalyzer

from Standards import (MDStageNames,
                       Fields,
                       MDRecords)

import copy

import textwrap

import subprocess

from oeommtools import utils as oeommutils

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

from YankCubes.yank_templates import (resources,
                                      max_cube_running_time)


class YankSolvationFECube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "YankSolvationFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = [["Alchemical free energy calculations"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Total Number of Yank iterations for the entire floe. "
                  "A Yank iteration is 500 MD steps")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of MD steps per iteration")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    ligand_res_name = parameter.StringParameter(
        'ligand_res_name',
        required=True,
        default='LIG',
        help_text='Ligand residue name')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text="Print verbose YANK logging output")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='On enables Hydrogen Mass Repartitioning')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            self.log.warn(">>>>>>> {} verbose {}".format(self.title, self.opt['verbose']))
            self.log.warn(">>>>>>> {} iterations {}".format(self.title, self.opt['iterations']))
            self.log.warn(">>>>>>> {} pressure {}".format(self.title, self.opt['pressure']))
            self.log.warn(">>>>>>> {} temperature {}".format(self.title, self.opt['temperature']))

            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            current_iteration_field = OEField("current_iterations", Types.Int)

            # Logger string
            str_logger = '-' * 32 + ' YANK SOLV CUBE PARAMETERS ' + '-' * 32
            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing the Primary Molecule field")

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            opt['system_title'] = system_title

            if not record.has_value(Fields.id):
                raise ValueError("Missing the ID field")

            opt['system_id'] = record.get_value(Fields.id)

            if opt['iterations'] <= 0:
                raise ValueError("The number of iterations cannot be a non-negative number: {}".format(opt['iterations']))

            if not record.has_value(current_iteration_field):
                raise ValueError("The current number of iterations has not been defined")

            # Current number of iterations
            current_iterations = record.get_value(current_iteration_field)
                
            if current_iterations == 0:
                minimize = True
                resume_sim = False
                resume_setup = False

            total_time_per_iteration = resources['k80']['w29']['intercept'] + resources['k80']['w29']['slope'] * system.NumAtoms()

            iterations_per_cube = int(max_cube_running_time/total_time_per_iteration)
            # # TODO DEBUGGING REMOVE NEXT LINE
            # iterations_per_cube = 5

            # Calculate the new number of iterations to run
            if current_iterations + iterations_per_cube > opt['iterations']:
                new_iterations = opt['iterations']
            else:
                new_iterations = current_iterations + iterations_per_cube
                
            # Checkpoint interval
            checkpoint_interval = new_iterations - current_iterations
                        
            self.log.warn(">>>>>>> {} current iterations {}".format(self.title, current_iterations))
            self.log.warn(">>>>>>> {} iterations per cube {}".format(self.title, iterations_per_cube))
            self.log.warn(">>>>>>> {} new_iterations {}".format(self.title, new_iterations))
            self.log.warn(">>>>>>> {} checkpoint_interval {}".format(self.title, checkpoint_interval))

            prot_split, lig_split, water, excipients = oeommutils.split(system, ligand_res_name=opt['ligand_res_name'])

            fchg_lig = 0
            for at in lig_split.GetAtoms():
                fchg_lig += at.GetFormalCharge()

            if fchg_lig != 0:
                alchemical_pme_treatment = 'exact'
            else:
                alchemical_pme_treatment = 'direct-space'

            if fchg_lig != 0:
                raise ValueError("The ligand is not neutral: formal charge = {}".format(fchg_lig))

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(system) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(system.GetTitle(), new_args))
                opt.update(new_args)

            if not record.has_value(Fields.md_stages):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the the Parmed structure
            parmed_structure = md_system_record.get_value(Fields.structure)

            # Extract the MD data
            mdData = omm_utils.MDData(parmed_structure)
            solvated_structure = mdData.structure

            # Extract the ligand parmed structure
            solute_structure = solvated_structure.split()[0][0]
            solute_structure.box = None

            solvent_res_names = set()
            for res in solvated_structure.residues:
                solvent_res_names.add(res.name)
            solvent_res_names.remove(solute_structure.residues[0].name)

            solvent_res_names = list(solvent_res_names)

            for i in range(0, len(solvent_res_names)):
                if '+' in solvent_res_names[i]:
                    solvent_res_names[i] = "'"+solvent_res_names[i]+"'"
                if '-' in solvent_res_names[i]:
                    solvent_res_names[i] = "'"+solvent_res_names[i]+"'"

            solvent_str_names = ' '.join(solvent_res_names)
            
            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                solvated_structure_fn = os.path.join(output_directory, "solvated.pdb")
                solute_structure_fn = os.path.join(output_directory, "solute.pdb")

                solvated_omm_sys_serialized_fn = os.path.join(output_directory, "solvated.xml")
                solute_omm_sys_serialized_fn = os.path.join(output_directory, "solute.xml")

                # Restarting
                if current_iterations != 0:
                    yank_files = md_stage_record.get_value(Fields.trajectory)
                    filename = omm_utils.download(yank_files, delete=True)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)

                    # Remove file after extraction
                    os.remove(filename)

                    # Disable minimization if restart is enabled
                    minimize = False
                    # Enable Yank Restarting
                    resume_sim = True
                    resume_setup = True
                else:
                    with open(solvated_structure_fn, 'w') as f:
                        app.PDBFile.writeFile(solvated_structure.topology, solvated_structure.positions, file=f)

                    with open(solute_structure_fn, 'w') as f:
                        app.PDBFile.writeFile(solute_structure.topology, solute_structure.positions, file=f)

                    # Create the solvated and vacuum system
                    solvated_omm_sys = solvated_structure.createSystem(nonbondedMethod=app.PME,
                                                                       nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                       constraints=app.HBonds,
                                                                       removeCMMotion=False)

                    solute_omm_sys = solute_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                                   constraints=app.HBonds,
                                                                   removeCMMotion=False)

                    solvated_omm_sys_serialized = XmlSerializer.serialize(solvated_omm_sys)
                    with open(solvated_omm_sys_serialized_fn, 'w') as solvated_f:
                        solvated_f.write(solvated_omm_sys_serialized)

                    solute_omm_sys_serialized = XmlSerializer.serialize(solute_omm_sys)
                    with open(solute_omm_sys_serialized_fn, 'w') as solute_f:
                        solute_f.write(solute_omm_sys_serialized)

                yank_template = yank_solvation_template.format(
                                                 verbose='yes' if opt['verbose'] else 'no',
                                                 minimize='yes' if minimize else 'no',
                                                 output_directory=output_directory,
                                                 timestep=4.0 if opt['hmr'] else 2.0,
                                                 nsteps_per_iteration=opt['nsteps_per_iteration'],
                                                 number_iterations=new_iterations,
                                                 temperature=opt['temperature'],
                                                 pressure=opt['pressure'],
                                                 resume_sim='yes' if resume_sim else 'no',
                                                 resume_setup='yes' if resume_setup else 'no',
                                                 hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                                                 alchemical_pme_treatment=alchemical_pme_treatment,
                                                 checkpoint_interval=checkpoint_interval,
                                                 solvated_pdb_fn=solvated_structure_fn,
                                                 solvated_xml_fn=solvated_omm_sys_serialized_fn,
                                                 solute_pdb_fn=solute_structure_fn,
                                                 solute_xml_fn=solute_omm_sys_serialized_fn,
                                                 solvent_dsl=solvent_str_names)

                opt['yank_template'] = yank_template

                yankutils.run_yank(opt)

                if new_iterations == opt['iterations']:

                    exp_dir = os.path.join(output_directory, "experiments")

                    experiment_to_analyze = ExperimentAnalyzer(exp_dir)
                    analysis = experiment_to_analyze.auto_analyze()

                    # Calculate solvation free energy and its error
                    DeltaG_solvation = analysis['free_energy']['free_energy_diff_unit'].\
                                           in_units_of(unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole
                    dDeltaG_solvation = analysis['free_energy']['free_energy_diff_error_unit'].\
                                            in_units_of(unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole

                    # Create OE Field to save the Solvation Free Energy in kcal/mol
                    DG_Field = OEField('Solvation FE', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                    dG_Field = OEField('Solvation FE Error', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                    record.set_value(DG_Field, DeltaG_solvation)
                    record.set_value(dG_Field, dDeltaG_solvation)

                    opt_1 = '--store={}'.format(exp_dir)

                    result_fn = os.path.join(output_directory, 'results.html')
                    opt_2 = '--output={}'.format(result_fn)

                    opt_3 = '--format=html'

                    try:
                        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

                        subprocess.check_call(['yank', 'analyze', 'report', opt_1, opt_2, opt_3])

                        with open(result_fn, 'r') as f:
                            result_str = f.read()

                        record.set_value(Fields.yank_analysis, result_str)

                        # Upload Floe Report
                        if in_orion():
                            session = OrionSession()

                            file_upload = File.upload(session,
                                                      "{} Yank Report".format(opt['system_title']),
                                                      result_fn)

                            session.tag_resource(file_upload, "floe_report")

                            job_id = environ.get('ORION_JOB_ID')

                            if job_id:
                                session.tag_resource(file_upload, "Job {}".format(job_id))

                    except subprocess.SubprocessError:
                        opt['Logger'].warn("The result file have not been generated")

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload(tar_fn)

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n'+flog.read()

                md_stage_record = MDRecords.MDStageRecord(MDStageNames.FEC,
                                                          MDRecords.MDSystemRecord(system, mdData.structure),
                                                          log=str_logger,
                                                          trajectory=lf)
                # md_stages.append(md_stage_record)
                md_stages[-1] = md_stage_record

                record.set_value(Fields.md_stages, md_stages)

                record.set_value(Fields.primary_molecule, system)

                record.set_value(current_iteration_field, new_iterations)

                self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class SyncBindingFECube(OERecordComputeCube):
    version = "0.0.0"
    title = "SyncSolvationFECube"
    description = """
    This cube is used to synchronize the solvated ligands and the related
    solvated complexes
    """
    classification = [["Synchronization Cube"]]
    tags = [tag for lists in classification for tag in lists]

    solvated_ligand_in_port = RecordInputPort("solvated_ligand_in_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.solvated_ligand_list = []
        self.solvated_complex_list = []

    def process(self, record, port):

        try:

            if port == 'solvated_ligand_in_port':
                self.solvated_ligand_list.append(record)
            else:
                self.solvated_complex_list.append(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

    def end(self):
        try:

            solvated_complex_lig_list = [(i, j) for i, j in
                                         itertools.product(self.solvated_ligand_list, self.solvated_complex_list)
                                         if i.get_value(Fields.id) == j.get_value(Fields.id)]

            for pair in solvated_complex_lig_list:

                new_record = OERecord()

                self.log.info("SYNC = ({} - {} , {} - {})".format(pair[0].get_value(Fields.title),
                                                                  pair[0].get_value(Fields.id),
                                                                  pair[1].get_value(Fields.title),
                                                                  pair[1].get_value(Fields.id)))

                complex_solvated = pair[1].get_value(Fields.primary_molecule)
                new_record.set_value(Fields.primary_molecule, complex_solvated)
                new_record.set_value(Fields.id, pair[1].get_value(Fields.id))
                new_record.set_value(Fields.title, pair[1].get_value(Fields.title))

                ligand_solvated_field = OEField("ligand_solvated", Types.Record)
                complex_solvated_field = OEField("complex_solvated", Types.Record)
                new_record.set_value(ligand_solvated_field, pair[0].get_value(Fields.md_stages)[-1].get_value(Fields.md_system))
                new_record.set_value(complex_solvated_field, pair[1].get_value(Fields.md_stages)[-1].get_value(Fields.md_system))

                self.success.emit(new_record)

        except:
            self.log.error(traceback.format_exc())

        return


class YankBindingFECube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "YankBindingFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = [["Alchemical free energy calculations"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Total Number of Yank iterations for the entire floe. "
                  "A Yank iteration is 500 MD steps")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of MD steps per iteration")

    timestep = parameter.DecimalParameter(
        'timestep',
        default=2.0,
        help_text="Timestep (fs)")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    ligand_resname = parameter.StringParameter(
        'ligand_resname',
        default='LIG',
        help_text='The decoupling ligand residue name')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text="Print verbose YANK logging output")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='On enables Hydrogen Mass Repartitioning')

    sampler = parameter.StringParameter(
        'sampler',
        default='repex',
        choices=['repex', 'sams'],
        help_text='Yank Sampling mode: REPEX Replica Exchange and SAMS Self-Adjusted Mixture Sampling')

    restraints = parameter.StringParameter(
        'restraints',
        default='boresch',
        choices=['harmonic', 'boresch'],
        help_text='Select the restraint types to apply to the ligand during the '
                  'alchemical decoupling. Choices: harmonic, boresch')

    protocol = parameter.StringParameter(
        'protocol',
        default='windows_29',
        choices=['windows_29',
                 'auto_protocol',
                 'windows_sams'],
        help_text='Select the protocol type')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = dict(self.opt)

            current_iteration_field = OEField("current_iterations", Types.Int)

            self.log.warn(">>>>>>> {} verbose {}".format(self.title, self.opt['verbose']))
            self.log.warn(">>>>>>> {} iterations {}".format(self.title, self.opt['iterations']))
            self.log.warn(">>>>>>> {} pressure {}".format(self.title, self.opt['pressure']))
            self.log.warn(">>>>>>> {} temperature {}".format(self.title, self.opt['temperature']))

            # Logger string
            str_logger = '-' * 32 + ' YANK BINDING CUBE PARAMETERS ' + '-' * 32
            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("The Primary Molecule is missing field")
        
            complex = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = complex.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            opt['system_title'] = system_title
            
            if not record.has_value(Fields.id):
                raise ValueError("Missing the ID field")

            opt['system_id'] = record.get_value(Fields.id)

            if opt['iterations'] <= 0:
                raise ValueError("The number of iterations cannot be a non-negative number: {}".format(opt['iterations']))

            if not record.has_value(current_iteration_field):
                raise ValueError("The current number of iterations has not been defined")

            # Current number of iterations
            current_iterations = record.get_value(current_iteration_field)

            if current_iterations == 0:
                minimize = True
                resume_sim = False
                resume_setup = False

            if opt['sampler'] == 'repex':
                self.log.warn('REPEX samplig selected')
                total_time_per_iteration = resources['k80']['w29']['intercept'] + resources['k80']['w29']['slope'] * complex.NumAtoms()
            elif opt['sampler'] == 'sams':
                self.log.warn('SAMS samplig selected')
                total_time_per_iteration = resources['k80']['wsams']['intercept'] + resources['k80']['wsams']['slope'] * complex.NumAtoms()
            else:
                raise ValueError("The selected sampling method is not supported: {}".format(opt['sampler']))

            iterations_per_cube = int(max_cube_running_time/total_time_per_iteration)

            # # TODO DEBUGGING REMOVE NEXT LINE
            # iterations_per_cube = 5

            # Calculate the new number of iterations to run
            if current_iterations + iterations_per_cube > opt['iterations']:
                new_iterations = opt['iterations']
            else:
                new_iterations = current_iterations + iterations_per_cube

            # Checkpoint interval
            checkpoint_interval = new_iterations - current_iterations

            self.log.warn(">>>>>>> {} current iterations {}".format(self.title, current_iterations))
            self.log.warn(">>>>>>> {} iterations per cube {}".format(self.title, iterations_per_cube))
            self.log.warn(">>>>>>> {} new_iterations {}".format(self.title, new_iterations))
            self.log.warn(">>>>>>> {} checkpoint_interval {}".format(self.title, checkpoint_interval))

            # Split the complex in components
            protein_split, ligand_split, water, excipients = oeommutils.split(complex,
                                                                              ligand_res_name=self.opt['ligand_resname'])
            fchg_lig = 0
            for at in ligand_split.GetAtoms():
                fchg_lig += at.GetFormalCharge()

            if fchg_lig != 0:
                alchemical_pme_treatment = 'exact'
            else:
                alchemical_pme_treatment = 'direct-space'

            solvent = water.CreateCopy()

            if not oechem.OEAddMols(solvent, excipients):
                raise ValueError("Solvent merging failure")

            solvent_res_names = set()

            hv = oechem.OEHierView(solvent)
            for chain in hv.GetChains():
                for frag in chain.GetFragments():
                    for hres in frag.GetResidues():
                        solvent_res_names.add(hres.GetOEResidue().GetName())

            solvent_res_names = list(solvent_res_names)

            for i in range(0, len(solvent_res_names)):
                if '+' in solvent_res_names[i]:
                    solvent_res_names[i] = "'"+solvent_res_names[i]+"'"
                if '-' in solvent_res_names[i]:
                    solvent_res_names[i] = "'"+solvent_res_names[i]+"'"

            solvent_str_names = ' '.join(solvent_res_names)

            if current_iterations == 0:
                solvated_ligand_record_field = OEField("ligand_solvated", Types.Record)

                if not record.has_value(solvated_ligand_record_field):
                    raise ValueError("The parametrized solvated ligand is missing")

                solvated_ligand_mdsystem_record = record.get_value(solvated_ligand_record_field)

                solvated_ligand_parmed_structure = solvated_ligand_mdsystem_record.get_value(Fields.structure)

                solvated_complex_record_field = OEField("complex_solvated", Types.Record)

                if not record.has_value(solvated_complex_record_field):
                    raise ValueError("The parametrized solvated complex is missing")

                solvated_complex_mdsystem_record = record.get_value(solvated_complex_record_field)

                solvated_complex_parmed_structure = solvated_complex_mdsystem_record.get_value(Fields.structure)

            else:
                # Extract the MDStageRecord list
                md_stages = record.get_value(Fields.md_stages)

                # Extract the most recent MDStageRecord
                md_stage_record = md_stages[-1]

                complex_mdsystem_record = md_stage_record.get_value(Fields.md_system)

                solvated_complex_parmed_structure = complex_mdsystem_record.get_value(Fields.structure)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                solvated_complex_structure_fn = os.path.join(output_directory, "complex.pdb")
                solvated_ligand_structure_fn = os.path.join(output_directory, "solvent.pdb")
                solvated_complex_omm_serialized_fn = os.path.join(output_directory, "complex.xml")
                solvated_ligand_omm_serialized_fn = os.path.join(output_directory, "solvent.xml")

                if current_iterations != 0:
                    yank_files = md_stage_record.get_value(Fields.trajectory)
                    filename = omm_utils.download(yank_files, delete=True)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)

                    # Remove File after extraction
                    os.remove(filename)

                    # Disable minimization if restart is enabled
                    minimize = False

                    # Enable Yank Restarting
                    resume_sim = True
                    resume_setup = True
                else:

                    with open(solvated_complex_structure_fn, 'w') as f:
                        app.PDBFile.writeFile(solvated_complex_parmed_structure.topology,
                                              solvated_complex_parmed_structure.positions,
                                              file=f)

                    with open(solvated_ligand_structure_fn, 'w') as f:
                        app.PDBFile.writeFile(solvated_ligand_parmed_structure.topology,
                                              solvated_ligand_parmed_structure.positions,
                                              file=f)

                    # Create the solvated OpenMM systems
                    solvated_complex_omm_sys = solvated_complex_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                              ewaldErrorTolerance=1.0e-4,
                                                                                              nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                                              constraints=app.HBonds,
                                                                                              removeCMMotion=False)

                    solvated_ligand_omm_sys = solvated_ligand_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                            ewaldErrorTolerance=1.0e-4,
                                                                                            nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                                            constraints=app.HBonds,
                                                                                            removeCMMotion=False)

                    solvated_complex_omm_serialized = XmlSerializer.serialize(solvated_complex_omm_sys)

                    with open(solvated_complex_omm_serialized_fn, 'w') as solvated_complex_f:
                        solvated_complex_f.write(solvated_complex_omm_serialized)

                    solvated_ligand_omm_serialized = XmlSerializer.serialize(solvated_ligand_omm_sys)

                    with open(solvated_ligand_omm_serialized_fn, 'w') as solvated_ligand_f:
                        solvated_ligand_f.write(solvated_ligand_omm_serialized)

                yank_template = yank_binding_template.format(
                    verbose='yes' if opt['verbose'] else 'no',
                    minimize='yes' if minimize else 'no',
                    output_directory=output_directory,
                    timestep=opt['timestep'],
                    nsteps_per_iteration=opt['nsteps_per_iteration'],
                    number_iterations=new_iterations,
                    temperature=opt['temperature'],
                    pressure=opt['pressure'],
                    resume_sim='yes' if resume_sim else 'no',
                    resume_setup='yes' if resume_setup else 'no',
                    hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                    alchemical_pme_treatment=alchemical_pme_treatment,
                    checkpoint_interval=checkpoint_interval,
                    complex_pdb_fn=solvated_complex_structure_fn,
                    complex_xml_fn=solvated_complex_omm_serialized_fn,
                    solvent_pdb_fn=solvated_ligand_structure_fn,
                    solvent_xml_fn=solvated_ligand_omm_serialized_fn,
                    ligand_resname=opt['ligand_resname'],
                    solvent_dsl=solvent_str_names,
                    sampler=opt['sampler'],
                    restraints=opt['restraints'],
                    protocol=opt['protocol']
                )

                opt['yank_template'] = yank_template

                yankutils.run_yank(opt)

                if new_iterations == opt['iterations']:
                    exp_dir = os.path.join(output_directory, "experiments")

                    experiment_to_analyze = ExperimentAnalyzer(exp_dir)
                    analysis = experiment_to_analyze.auto_analyze()

                    # Calculate binding free energy and its error in kcal/mol
                    DeltaG_binding = analysis['free_energy']['free_energy_diff_unit'].\
                                         in_units_of(unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole
                    dDeltaG_binding = analysis['free_energy']['free_energy_diff_error_unit'].\
                                          in_units_of(unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole

                    # Create OE Field to save the Solvation Free Energy in kcal/mol
                    DG_Field = OEField('Binding Affinity', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                    dG_Field = OEField('Binding Affinity Error', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                    record.set_value(DG_Field, DeltaG_binding)
                    record.set_value(dG_Field, dDeltaG_binding)

                    opt_1 = '--store={}'.format(exp_dir)

                    result_fn = os.path.join(output_directory, 'results.html')
                    opt_2 = '--output={}'.format(result_fn)

                    opt_3 = '--format=html'

                    try:

                        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

                        subprocess.check_call(['yank', 'analyze', 'report', opt_1, opt_2, opt_3])

                        with open(result_fn, 'r') as f:
                            result_str = f.read()

                        record.set_value(Fields.yank_analysis, result_str)

                        # Upload Floe Report
                        if in_orion():
                            session = OrionSession()

                            file_upload = File.upload(session,
                                                      "{} Yank Report".format(opt['system_title']),
                                                      result_fn)

                            session.tag_resource(file_upload, "floe_report")

                            job_id = environ.get('ORION_JOB_ID')

                            if job_id:
                                session.tag_resource(file_upload, "Job {}".format(job_id))

                    except subprocess.SubprocessError:
                        opt['Logger'].warn("The result file have not been generated")

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory+"_"+opt['system_title']) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload(tar_fn)

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n' + flog.read()

                md_stage_record = MDRecords.MDStageRecord(MDStageNames.FEC,
                                                          MDRecords.MDSystemRecord(complex,
                                                                                   solvated_complex_parmed_structure),
                                                          log=str_logger,
                                                          trajectory=lf)
                if current_iterations != 0:
                    # md_stages.append(md_stage_record)
                    md_stages[-1] = md_stage_record
                    record.set_value(Fields.md_stages, md_stages)
                else:
                    record = OERecord()
                    record.set_value(Fields.md_stages, [md_stage_record])
                    record.set_value(Fields.id, opt['system_id'])
                    record.set_value(Fields.title, opt['system_title'])

            record.set_value(current_iteration_field, new_iterations)
            record.set_value(Fields.primary_molecule, complex)

            self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class YankProxyCube(OERecordComputeCube):
    version = "0.0.0"
    title = "YankProxyCube"
    description = """
    This cube is used to implement a cycle with the Yank Solvation FE and
    Yank Binding Cubes.
    """
    classification = [["Alchemical free energy calculations"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    iterations = parameter.IntegerParameter(
        'iterations',
        default=13,
        help_text="Total number of iterations")

    cycle_out_port = RecordOutputPort("cycle_out_port")
    cycle_in_port = RecordInputPort("cycle_in_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        current_iteration_field = OEField("current_iterations", Types.Int)

        try:
            if port == 'intake':
                record.set_value(current_iteration_field, 0)
                current_iteration = record.get_value(current_iteration_field)
                self.opt['Logger'].warn(">>>>>>>>>>> PROXY current iterations {}".format(current_iteration))
                self.opt['Logger'].warn(">>>>>>>>>>> PROXY max iterations {}".format(self.opt['iterations']))
                self.opt['Logger'].warn(">>>>>>>>>>> Forwarding to Cycle...")
                self.cycle_out_port.emit(record)
            else:  # Cycle Port
                current_iteration = record.get_value(current_iteration_field)

                self.opt['Logger'].warn(">>>>>>>>>>> PROXY current iterations {}".format(current_iteration))
                self.opt['Logger'].warn(">>>>>>>>>>> PROXY max iterations {}".format(self.opt['iterations']))

                if current_iteration == self.opt['iterations']:
                    self.opt['Logger'].warn(">>>>>>>>>>> Finishing...")
                    self.success.emit(record)
                else:
                    self.opt['Logger'].warn(">>>>>>>>>>> Forwarding to Cycle...")
                    self.cycle_out_port.emit(record)
        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return
