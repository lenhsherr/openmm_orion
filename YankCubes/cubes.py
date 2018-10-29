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

import MDCubes.utils as omm_utils

from openeye import oechem

from simtk.openmm import (app,
                          unit,
                          XmlSerializer)

from tempfile import TemporaryDirectory

import tarfile

import os

import itertools

from YankCubes import utils as yankutils


from Standards import (MDStageTypes,
                       Fields,
                       MDRecords)

import copy

import textwrap

from oeommtools import utils as oeommutils

from Standards.utils import ParmedData

from MDCubes.utils import MDState


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

    lig_res_name = parameter.StringParameter(
        'lig_res_name',
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
            self.log.info("[{}] verbose {}".format(self.title, self.opt['verbose']))
            self.log.info("[{}] iterations {}".format(self.title, self.opt['iterations']))
            self.log.info("[{}] pressure {}".format(self.title, self.opt['pressure']))
            self.log.info("[{}] temperature {}".format(self.title, self.opt['temperature']))

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
                raise ValueError(
                    "The number of iterations cannot be a non-negative number: {}".format(opt['iterations']))

            if not record.has_value(current_iteration_field):
                raise ValueError("The current number of iterations has not been defined")

            # Current number of iterations
            current_iterations = record.get_value(current_iteration_field)

            if current_iterations == 0:
                opt['minimize'] = True
                opt['resume_sim'] = False
                opt['resume_setup'] = False
                iterations_per_cube = 10
            else:

                iterations_per_cube_field = OEField("iterations_per_cube", Types.Int)

                if not record.has_value(iterations_per_cube_field):
                    raise ValueError("Missing the number of iterations per cube field")

                iterations_per_cube = record.get_value(iterations_per_cube_field)
                self.log.info("{} iterations per cube {}".format(self.title, iterations_per_cube))

            # Calculate the new number of iterations to run
            if current_iterations + iterations_per_cube > opt['iterations']:
                opt['new_iterations'] = opt['iterations']
            else:
                opt['new_iterations'] = current_iterations + iterations_per_cube

            # Checkpoint interval
            opt['checkpoint_interval'] = opt['new_iterations'] - current_iterations

            self.log.info("[{}] current iterations {}".format(self.title, current_iterations))
            self.log.info("[{}] new_iterations {}".format(self.title, opt['new_iterations']))
            self.log.info("[{}] checkpoint_interval {}".format(self.title, opt['checkpoint_interval']))

            prot_split, lig_split, water, excipients = oeommutils.split(system, ligand_res_name=opt['lig_res_name'])

            fchg_lig = 0
            for at in lig_split.GetAtoms():
                fchg_lig += at.GetFormalCharge()

            if fchg_lig != 0:
                opt['alchemical_pme_treatment'] = 'exact'
            else:
                opt['alchemical_pme_treatment'] = 'direct-space'

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

            # Extract from the record the Parmed structure
            if not record.has_value(Fields.pmd_structure):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            solvated_structure = record.get_value(Fields.pmd_structure)

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract the MD State
            mdstate = md_system_record.get_value(Fields.md_state)

            # Update the Parmed structure with the MD State
            solvated_structure.positions = mdstate.get_positions()
            solvated_structure.velocities = mdstate.get_velocities()
            solvated_structure.box_vectors = mdstate.get_box_vectors()

            # Extract the ligand Parmed structure
            solute_structure = solvated_structure.split()[0][0]
            solute_structure.box = None

            solvent_res_names = set()
            for res in solvated_structure.residues:
                solvent_res_names.add(res.name)
            solvent_res_names.remove(solute_structure.residues[0].name)

            solvent_res_names = list(solvent_res_names)

            for i in range(0, len(solvent_res_names)):
                if '+' in solvent_res_names[i]:
                    solvent_res_names[i] = "'" + solvent_res_names[i] + "'"
                if '-' in solvent_res_names[i]:
                    solvent_res_names[i] = "'" + solvent_res_names[i] + "'"

            opt['solvent_str_names'] = ' '.join(solvent_res_names)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                opt['output_directory'] = output_directory
                opt['solvated_structure_fn'] = os.path.join(output_directory, "solvated.pdb")
                opt['solute_structure_fn'] = os.path.join(output_directory, "solute.pdb")

                opt['solvated_omm_sys_serialized_fn'] = os.path.join(output_directory, "solvated.xml")
                opt['solute_omm_sys_serialized_fn'] = os.path.join(output_directory, "solute.xml")

                # Restarting
                if current_iterations != 0:

                    if md_stage_record.has_value(Fields.trajectory):
                        yank_files = md_stage_record.get_value(Fields.trajectory)

                    elif md_stage_record.has_value(Fields.orion_local_trj_field):
                        yank_files = md_stage_record.get_value(Fields.orion_local_trj_field)

                    else:
                        print("No Yank trajectory file have been found in the selected stage record {}".format(
                            md_stage_record.get_value(Fields.stage_name)))

                    filename = os.path.join(output_directory, "yank_files.tar")

                    fn_local = omm_utils.download_file(yank_files, filename, delete=True)

                    with tarfile.open(fn_local) as tar:
                        tar.extractall(path=output_directory)

                    # remove the tar file
                    os.remove(fn_local)

                    # Disable minimization if restart is enabled
                    opt['minimize'] = False
                    # Enable Yank Restarting
                    opt['resume_sim'] = True
                    opt['resume_setup'] = True
                else:
                    with open(opt['solvated_structure_fn'], 'w') as f:
                        app.PDBFile.writeFile(solvated_structure.topology, solvated_structure.positions, file=f)

                    with open(opt['solute_structure_fn'], 'w') as f:
                        app.PDBFile.writeFile(solute_structure.topology, solute_structure.positions, file=f)

                    # Create the solvated and vacuum system
                    solvated_omm_sys = solvated_structure.createSystem(nonbondedMethod=app.PME,
                                                                       nonbondedCutoff=opt[
                                                                                           'nonbondedCutoff'] * unit.angstroms,
                                                                       constraints=app.HBonds,
                                                                       removeCMMotion=False)

                    solute_omm_sys = solute_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                                   constraints=app.HBonds,
                                                                   removeCMMotion=False)

                    solvated_omm_sys_serialized = XmlSerializer.serialize(solvated_omm_sys)
                    with open(opt['solvated_omm_sys_serialized_fn'], 'w') as solvated_f:
                        solvated_f.write(solvated_omm_sys_serialized)

                    solute_omm_sys_serialized = XmlSerializer.serialize(solute_omm_sys)
                    with open(opt['solute_omm_sys_serialized_fn'], 'w') as solute_f:
                        solute_f.write(solute_omm_sys_serialized)

                # Run Yank
                yankutils.run_yank_solvation(opt)

                if current_iterations == 0:

                    iterations_per_cube = yankutils.calculate_iteration_time(output_directory, opt['new_iterations'])

                    if iterations_per_cube == 0:
                        raise ValueError("Total running time per cube > max Orion running time per cube")

                    iterations_per_cube_field = OEField("iterations_per_cube", Types.Int)

                    record.set_value(iterations_per_cube_field, iterations_per_cube)

                    self.log.info("{} iterations per cube saved on the record: {}".format(self.title,
                                                                                          iterations_per_cube))

                # Run Yank analysis
                if opt['new_iterations'] == opt['iterations']:

                    DeltaG_solvation, dDeltaG_solvation = yankutils.run_yank_analysis(opt)

                    # Create OE Field to save the Solvation Free Energy in kcal/mol
                    DG_Field = OEField('Solvation FE', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                    dG_Field = OEField('Solvation FE Error', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                    record.set_value(DG_Field, DeltaG_solvation)
                    record.set_value(dG_Field, dDeltaG_solvation)

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload_file(tar_fn, opt['system_title']+'.tar.gz')

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n' + flog.read()

                md_stage_record = MDRecords.MDStageRecord(self.title, MDStageTypes.FEC,
                                                          MDRecords.MDSystemRecord(system, mdstate),
                                                          log=str_logger,
                                                          trajectory=lf)
                # md_stages.append(md_stage_record)
                md_stages[-1] = md_stage_record

                record.set_value(Fields.md_stages, md_stages)

                record.set_value(Fields.primary_molecule, system)

                record.set_value(current_iteration_field, opt['new_iterations'])

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

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

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

                ligand_solvated_field = OEField("ligand_pmd_solvated", ParmedData)
                complex_solvated_field = OEField("complex_pmd_solvated", ParmedData)

                # Extract and update ligand parmed structure
                ligand_md_state = pair[0].get_value(Fields.md_stages)[-1].get_value(Fields.md_system).get_value(Fields.md_state)
                ligand_pmd_structure = pair[0].get_value(Fields.pmd_structure)

                ligand_pmd_structure.positions = ligand_md_state.get_positions()
                ligand_pmd_structure.velocities = ligand_md_state.get_velocities()
                ligand_pmd_structure.box_vectors = ligand_md_state.get_box_vectors()

                # Extract and update complex parmed structure
                complex_md_state = pair[1].get_value(Fields.md_stages)[-1].get_value(Fields.md_system).get_value(Fields.md_state)
                complex_pmd_structure = pair[1].get_value(Fields.pmd_structure)

                complex_pmd_structure.positions = complex_md_state.get_positions()
                complex_pmd_structure.velocities = complex_md_state.get_velocities()
                complex_pmd_structure.box_vectors = complex_md_state.get_box_vectors()

                new_record.set_value(ligand_solvated_field, ligand_pmd_structure)
                new_record.set_value(complex_solvated_field, complex_pmd_structure)

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

    lig_res_name = parameter.StringParameter(
        'lig_res_name',
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

    protocol_repex = parameter.StringParameter(
        'protocol_repex',
        default='windows_29',
        choices=['auto_protocol',
                 'windows_29'],
        help_text='Select the repex protocol type')

    protocol_sams = parameter.StringParameter(
        'protocol_sams',
        default='windows_sams',
        choices=['auto_protocol',
                 'windows_sams'],
        help_text='Select the sams protocol type')

    user_yank_yaml_file = parameter.FileInputParameter(
        'user_yank_yaml_file',
        default=None,
        help_text="Binding Affinity Yank yaml file. If a file is provided the "
                  "Yank yaml setting will be generated by using part of the Yank "
                  "yaml template parameters and others will be overwritten by the "
                  "provided yaml file")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = dict(self.opt)

            current_iteration_field = OEField("current_iterations", Types.Int)

            self.log.info("[{}] verbose {}".format(self.title, self.opt['verbose']))
            self.log.info("[{}] iterations {}".format(self.title, self.opt['iterations']))
            self.log.info("[{}] pressure {}".format(self.title, self.opt['pressure']))
            self.log.info("[{}] temperature {}".format(self.title, self.opt['temperature']))

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
                raise ValueError(
                    "The number of iterations cannot be a non-negative number: {}".format(opt['iterations']))

            if not record.has_value(current_iteration_field):
                raise ValueError("The current number of iterations has not been defined")

            # Current number of iterations
            current_iterations = record.get_value(current_iteration_field)

            if current_iterations == 0:
                opt['minimize'] = True
                opt['resume_sim'] = False
                opt['resume_setup'] = False
                iterations_per_cube = 10
            else:
                iterations_per_cube_field = OEField("iterations_per_cube", Types.Int)

                if not record.has_value(iterations_per_cube_field):
                    raise ValueError("Missing the number of iterations per cube field")

                iterations_per_cube = record.get_value(iterations_per_cube_field)
                self.log.info("{} iterations per cube {}".format(self.title, iterations_per_cube))

            # Calculate the new number of iterations to run
            if current_iterations + iterations_per_cube > opt['iterations']:
                opt['new_iterations'] = opt['iterations']
            else:
                opt['new_iterations'] = current_iterations + iterations_per_cube

            # Checkpoint interval
            opt['checkpoint_interval'] = opt['new_iterations'] - current_iterations

            self.log.info("[{}] current iterations {}".format(self.title, current_iterations))
            self.log.info("[{}] new_iterations {}".format(self.title, opt['new_iterations']))
            self.log.info("[{}] checkpoint_interval {}".format(self.title, opt['checkpoint_interval']))

            # Split the complex in components
            protein_split, ligand_split, water, excipients = oeommutils.split(complex,
                                                                              ligand_res_name=self.opt['lig_res_name'])

            fchg_lig = 0
            for at in ligand_split.GetAtoms():
                fchg_lig += at.GetFormalCharge()

            if fchg_lig != 0:
                opt['alchemical_pme_treatment'] = 'exact'
            else:
                opt['alchemical_pme_treatment'] = 'direct-space'

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
                    solvent_res_names[i] = "'" + solvent_res_names[i] + "'"
                if '-' in solvent_res_names[i]:
                    solvent_res_names[i] = "'" + solvent_res_names[i] + "'"

            opt['solvent_str_names'] = ' '.join(solvent_res_names)

            if current_iterations == 0:

                solvated_ligand_pmd_field = OEField("ligand_pmd_solvated", ParmedData)
                solvated_complex_pmd_field = OEField("complex_pmd_solvated", ParmedData)

                if not record.has_value(solvated_ligand_pmd_field):
                    raise ValueError("The Parmed solvated ligand structure is missing")

                solvated_ligand_parmed_structure = record.get_value(solvated_ligand_pmd_field)

                if not record.has_value(solvated_complex_pmd_field):
                    raise ValueError("The Parmed solvated complex structure is missing")

                solvated_complex_parmed_structure = record.get_value(solvated_complex_pmd_field)

                mdstate = MDState(solvated_complex_parmed_structure)

            else:
                # Extract the MDStageRecord list
                md_stages = record.get_value(Fields.md_stages)

                # Extract the most recent MDStageRecord
                md_stage_record = md_stages[-1]

                complex_mdsystem_record = md_stage_record.get_value(Fields.md_system)

                mdstate = complex_mdsystem_record.get_value(Fields.md_state)

                solvated_complex_parmed_structure = record.get_value(Fields.pmd_structure)

                solvated_complex_parmed_structure.positions = mdstate.get_positions()
                solvated_complex_parmed_structure.velocities = mdstate.get_velocities()
                solvated_complex_parmed_structure.box_vectors = mdstate.get_box_vectors()

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                opt['output_directory'] = output_directory

                opt['solvated_complex_structure_fn'] = os.path.join(output_directory, "bonded_state.pdb")
                opt['solvated_ligand_structure_fn'] = os.path.join(output_directory, "unbonded_state.pdb")
                opt['solvated_complex_omm_serialized_fn'] = os.path.join(output_directory, "bonded_state.xml")
                opt['solvated_ligand_omm_serialized_fn'] = os.path.join(output_directory, "unbonded_state.xml")

                if current_iterations != 0:

                    if md_stage_record.has_value(Fields.trajectory):
                        yank_files = md_stage_record.get_value(Fields.trajectory)

                    elif md_stage_record.has_value(Fields.orion_local_trj_field):
                        yank_files = md_stage_record.get_value(Fields.orion_local_trj_field)

                    else:
                        print("No Yank trajectory file have been found in the selected stage record {}".format(
                            md_stage_record.get_value(Fields.stage_name)))

                    filename = os.path.join(output_directory, "yank_files.tar")

                    fn_local = omm_utils.download_file(yank_files, filename, delete=True)

                    with tarfile.open(fn_local) as tar:
                        tar.extractall(path=output_directory)

                    # remove the file
                    os.remove(fn_local)

                    # Disable minimization if restart is enabled
                    opt['minimize'] = False

                    # Enable Yank Restarting
                    opt['resume_sim'] = True
                    opt['resume_setup'] = True
                else:

                    with open(opt['solvated_complex_structure_fn'], 'w') as f:
                        app.PDBFile.writeFile(solvated_complex_parmed_structure.topology,
                                              solvated_complex_parmed_structure.positions,
                                              file=f)

                    with open(opt['solvated_ligand_structure_fn'], 'w') as f:
                        app.PDBFile.writeFile(solvated_ligand_parmed_structure.topology,
                                              solvated_ligand_parmed_structure.positions,
                                              file=f)

                    # Create the solvated OpenMM systems
                    solvated_complex_omm_sys = solvated_complex_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                              ewaldErrorTolerance=1.0e-4,
                                                                                              nonbondedCutoff=opt[
                                                                                                                  'nonbondedCutoff'] * unit.angstroms,
                                                                                              constraints=app.HBonds,
                                                                                              removeCMMotion=False)

                    solvated_ligand_omm_sys = solvated_ligand_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                            ewaldErrorTolerance=1.0e-4,
                                                                                            nonbondedCutoff=opt[
                                                                                                                'nonbondedCutoff'] * unit.angstroms,
                                                                                            constraints=app.HBonds,
                                                                                            removeCMMotion=False)

                    solvated_complex_omm_serialized = XmlSerializer.serialize(solvated_complex_omm_sys)

                    with open(opt['solvated_complex_omm_serialized_fn'], 'w') as solvated_complex_f:
                        solvated_complex_f.write(solvated_complex_omm_serialized)

                    solvated_ligand_omm_serialized = XmlSerializer.serialize(solvated_ligand_omm_sys)

                    with open(opt['solvated_ligand_omm_serialized_fn'], 'w') as solvated_ligand_f:
                        solvated_ligand_f.write(solvated_ligand_omm_serialized)

                if opt['sampler'] == 'repex':
                    opt['protocol'] = opt['protocol_repex']
                elif opt['sampler'] == 'sams':
                    opt['protocol'] = opt['protocol_sams']
                else:
                    raise ValueError("The selected sampler method is not currently supported: {}".format(opt['sampler']))

                yankutils.run_yank_binding(opt)

                if current_iterations == 0:

                    iterations_per_cube = yankutils.calculate_iteration_time(output_directory, opt['new_iterations'])

                    if iterations_per_cube == 0:
                        raise ValueError("Total running time per cube > max Orion running time per cube")

                # Run the analysis
                if opt['new_iterations'] == opt['iterations']:

                    DeltaG_binding, dDeltaG_binding = yankutils.run_yank_analysis(opt)

                    # Create OE Field to save the solvation Free Energy in kcal/mol
                    DG_Field = OEField('Solvation FE', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                    dG_Field = OEField('Solvation FE Error', Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                    record.set_value(DG_Field, DeltaG_binding)
                    record.set_value(dG_Field, dDeltaG_binding)

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory + "_" + opt['system_title']) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload_file(tar_fn, opt['system_title']+'.tar.gz')

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n' + flog.read()

                md_stage_record = MDRecords.MDStageRecord(self.title, MDStageTypes.FEC,
                                                          MDRecords.MDSystemRecord(complex,
                                                                                   mdstate),
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
                    record.set_value(Fields.pmd_structure, solvated_complex_parmed_structure)

                    iterations_per_cube_field = OEField("iterations_per_cube", Types.Int)
                    record.set_value(iterations_per_cube_field, iterations_per_cube)
                    self.log.info("[{}] iterations per cube saved on the record: {}".format(self.title,
                                                                                            iterations_per_cube))

            record.set_value(current_iteration_field, opt['new_iterations'])
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
        default=1000,
        help_text="Total number of iterations")

    cycle_out_port = RecordOutputPort("cycle_out_port")
    cycle_in_port = RecordInputPort("cycle_in_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        current_iteration_field = OEField("current_iterations", Types.Int)

        try:
            if not record.has_value(Fields.title):
                self.opt['Logger'].warn("Missing record Title field")

                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("The Primary Molecule is missing field")

                complex = record.get_value(Fields.primary_molecule)
                system_title = complex.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            self.opt['Logger'].info("{} System Title {}".format(self.title, system_title))

            if port == 'intake':
                record.set_value(current_iteration_field, 0)
                current_iteration = record.get_value(current_iteration_field)
                self.opt['Logger'].info("{} current iterations {}".format(self.title, current_iteration))
                self.opt['Logger'].info("{}  max iterations {}".format(self.title, self.opt['iterations']))
                self.opt['Logger'].info("{} Forwarding to Cycle...".format(self.title))
                self.cycle_out_port.emit(record)
            else:  # Cycle Port
                current_iteration = record.get_value(current_iteration_field)

                self.opt['Logger'].info("{} PROXY current iterations {}".format(self.title, current_iteration))
                self.opt['Logger'].info("{} PROXY max iterations {}".format(self.title, self.opt['iterations']))

                if current_iteration == self.opt['iterations']:
                    self.opt['Logger'].info("{} Finishing...".format(self.title))
                    self.success.emit(record)
                else:
                    self.opt['Logger'].warn("{} Forwarding to Cycle...".format(self.title))
                    self.cycle_out_port.emit(record)
        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return