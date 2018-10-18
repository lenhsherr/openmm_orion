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

import MDCubes.utils as utils

from openeye import oechem

import MDCubes.OpenMMCubes.simtools as simtools

from Standards import (Fields,
                       MDRecords,
                       MDStageNames)

from MDCubes.utils import (MDState,
                           md_simulation)

import copy

import textwrap

import os

from oeommtools import utils as oeommutils


class OpenMMminimizeCube(ParallelMixin, OERecordComputeCube):
    title = 'Minimization Cube'

    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "Minimization"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """
    Minimize the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and minimize it.

    Input parameters:
    steps (integer): the number of steps of minimization to apply. If 0
    the minimization will proceed until convergence is reached
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    steps = parameter.IntegerParameter(
        'steps',
        default=0,
        help_text="""Number of minimization steps.
                  If 0 the minimization will continue
                  until convergence""")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=5.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around""")

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if the non-bonded method is NoCutoff""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    center = parameter.BooleanParameter(
        'center',
        default=True,
        description='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity')

    platform = parameter.StringParameter(
        'platform',
        default='Auto',
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='On enables Hydrogen Mass Repartitioning')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=True,
        help_text="""Save the MD simulation stage. If True the MD,
        simulation data will be appended to the md simulation stages 
        otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'min'
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title

            # Logger string
            str_logger = '-'*32 + ' MIN CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

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

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing the Primary Molecule field")

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                raise ValueError("Missing ID Field")

            opt['system_title'] = system_title
            opt['system_id'] = record.get_value(Fields.id)

            # Extract from the record the Parmed structure
            if not record.has_value(Fields.pmd_structure):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            parmed_structure = record.get_value(Fields.pmd_structure)

            # Extract from the record the md stages
            if not record.has_value(Fields.md_stages):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the topology and the MD State
            system = md_system_record.get_value(Fields.topology)
            mdstate = md_system_record.get_value(Fields.md_state)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(system) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                opt['Logger'].info("Updating parameters for molecule: {}\n{}".format(system.GetTitle(), new_args))
                opt.update(new_args)

            opt['molecule'] = system
            opt['str_logger'] = str_logger
            opt['Logger'].info('[{}] MINIMIZING System: {}'.format(opt['CubeTitle'], system_title))

            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the Parmed objects and the primary molecule
            parmed_structure.positions = new_mdstate.get_positions()

            if new_mdstate.get_box_vectors() is not None:
                parmed_structure.box_vectors = new_mdstate.get_box_vectors()

            if opt['SimType'] in ['nvt', 'npt']:
                # numpy array in units of angstrom/picoseconds
                parmed_structure.velocities = new_mdstate.get_velocities()

            # Update the OEMol complex positions to match the new
            # Parmed structure after the simulation
            new_temp_mol = oeommutils.openmmTop_to_oemol(parmed_structure.topology,
                                                         parmed_structure.positions,
                                                         verbose=False)
            new_pos = new_temp_mol.GetCoords()
            system.SetCoords(new_pos)

            record.set_value(Fields.primary_molecule, system)

            md_stage_record = MDRecords.MDStageRecord(MDStageNames.MINIMIZATION,
                                                      MDRecords.MDSystemRecord(system, new_mdstate),
                                                      log=opt['str_logger'])

            record.set_value(Fields.pmd_structure, parmed_structure)

            if opt['save_md_stage']:
                opt['Logger'].info("[{}] Saving MD stage: {}".format(opt['CubeTitle'], opt['SimType']))

                md_stages.append(md_stage_record)
            else:
                md_stages[-1] = md_stage_record

            record.set_value(Fields.md_stages, md_stages)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class OpenMMNvtCube(ParallelMixin, OERecordComputeCube):
    title = 'NVT Cube'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NVT"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """NVT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simulation
    at constant temperature and volume

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to warm up the complex.
      temperature (decimal): target temperature
    """

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

    time = parameter.DecimalParameter(
        'time',
        default=0.01,
        help_text="NVT simulation time in nanoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around""")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    trajectory_interval = parameter.DecimalParameter(
        'trajectory_interval',
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""")

    reporter_interval = parameter.DecimalParameter(
        'reporter_interval',
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""")

    suffix = parameter.StringParameter(
        'suffix',
        default='nvt',
        help_text='Filename suffix for output simulation files')

    center = parameter.BooleanParameter(
        'center',
        default=False,
        help_text='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        help_text='Increase log file verbosity')

    platform = parameter.StringParameter(
        'platform',
        default='Auto',
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        help_text='On enables Hydrogen Mass Repartitioning')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
        simulation data will be appended to the md simulation stages 
        otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'nvt'

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title

            # Logger string
            # Logger string
            str_logger = '-'*32 + ' NVT CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

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

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing the Primary Molecule field")

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                raise ValueError("Missing the ID field")

            opt['system_title'] = system_title
            opt['system_id'] = record.get_value(Fields.id)

            # Extract from the record the Parmed structure
            if not record.has_value(Fields.pmd_structure):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            parmed_structure = record.get_value(Fields.pmd_structure)

            # Extract from the record the md stages
            if not record.has_value(Fields.md_stages):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the topology and the MDState
            system = md_system_record.get_value(Fields.topology)
            mdstate = md_system_record.get_value(Fields.md_state)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(system) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                opt['Logger'].info("Updating parameters for molecule: {}\n{}".format(system.GetTitle(), new_args))
                opt.update(new_args)

            opt['outfname'] = '{}-{}'.format(system_title + '_'+str(opt['system_id']), opt['suffix'])
            opt['molecule'] = system
            opt['str_logger'] = str_logger

            opt['Logger'].info('[{}] START NVT SIMULATION: {}'.format(opt['CubeTitle'], system_title))
            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the Parmed objects and the primary molecule
            parmed_structure.positions = new_mdstate.get_positions()

            if new_mdstate.get_box_vectors() is not None:
                parmed_structure.box_vectors = new_mdstate.get_box_vectors()

            if opt['SimType'] in ['nvt', 'npt']:
                # numpy array in units of angstrom/picoseconds
                parmed_structure.velocities = new_mdstate.get_velocities()

            # Update the OEMol complex positions to match the new
            # Parmed structure after the simulation
            new_temp_mol = oeommutils.openmmTop_to_oemol(parmed_structure.topology,
                                                         parmed_structure.positions,
                                                         verbose=False)
            new_pos = new_temp_mol.GetCoords()
            system.SetCoords(new_pos)
            record.set_value(Fields.primary_molecule, system)

            # Trajectory
            if opt['trajectory_interval']:
                full_path = os.getcwd()
                filename = os.path.join(full_path, opt['outfname']+'.h5')
                lf = utils.upload_file(filename, opt['outfname']+'.h5')
            else:  # Empty Trajectory
                lf = None

            # Read in logging file if any
            if opt['reporter_interval']:

                with open(opt['outfname']+'.log', 'r') as flog:
                    str_logger = flog.read()

                opt['str_logger'] += '\n'+str_logger

            md_stage_record = MDRecords.MDStageRecord(MDStageNames.NVT,
                                                      MDRecords.MDSystemRecord(system, new_mdstate),
                                                      log=opt['str_logger'], trajectory=lf)

            if opt['save_md_stage']:
                opt['Logger'].info("[{}] Saving MD stage: {}".format(opt['CubeTitle'], opt['SimType']))

                md_stages.append(md_stage_record)
            else:
                md_stages[-1] = md_stage_record

            record.set_value(Fields.md_stages, md_stages)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class OpenMMNptCube(ParallelMixin, OERecordComputeCube):
    title = 'NPT Cube'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NPT"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """NPT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simulation at
    constant temperature and pressure.

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to perform the complex simulation.
      temperature (decimal): target temperature
      pressure (decimal): target pressure
    """

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

    time = parameter.DecimalParameter(
        'time',
        default=0.01,
        help_text="NPT simulation time in nanoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD simulation.
        Possible keywords are: ligand, protein, water, ions, ca_protein,
        cofactors. The selection can be refined by using logical tokens:
        not, noh, and, or, diff, around""")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    trajectory_interval = parameter.DecimalParameter(
        'trajectory_interval',
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""")

    reporter_interval = parameter.DecimalParameter(
        'reporter_interval',
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""")

    suffix = parameter.StringParameter(
        'suffix',
        default='npt',
        help_text='Filename suffix for output simulation files')

    center = parameter.BooleanParameter(
        'center',
        default=False,
        help_text='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        help_text='Increase log file verbosity.')

    platform = parameter.StringParameter(
        'platform',
        default='Auto',
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        help_text='On enables Hydrogen Mass Repartitioning')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
        simulation data will be appended to the md simulation stages 
        otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'npt'

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title
            # Logger string
            str_logger = '-'*32 + ' NPT CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

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

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing the Primary Molecule field")

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                raise ValueError("Missing the ID field")

            opt['system_title'] = system_title
            opt['system_id'] = record.get_value(Fields.id)

            # Extract from the record the Parmed structure
            if not record.has_value(Fields.pmd_structure):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            parmed_structure = record.get_value(Fields.pmd_structure)

            # Extract from the record the md stages
            if not record.has_value(Fields.md_stages):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the topology and the Parmed structure
            system = md_system_record.get_value(Fields.topology)
            mdstate = md_system_record.get_value(Fields.md_state)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(system) if dp.GetTag() in
                        ["temperature", "pressure"]}

            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                opt['Logger'].info("Updating parameters for molecule: {}\n{}".format(system.GetTitle(), new_args))

                opt.update(new_args)

            opt['outfname'] = '{}-{}'.format(system_title + '_' + str(opt['system_id']), opt['suffix'])
            opt['molecule'] = system
            opt['str_logger'] = str_logger
            opt['Logger'].info('[{}] START NPT SIMULATION: {}'.format(opt['CubeTitle'], system_title))

            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the Parmed objects and the primary molecule
            parmed_structure.positions = new_mdstate.get_positions()

            if new_mdstate.get_box_vectors() is not None:
                parmed_structure.box_vectors = new_mdstate.get_box_vectors()

            if opt['SimType'] in ['nvt', 'npt']:
                # numpy array in units of angstrom/picoseconds
                parmed_structure.velocities = new_mdstate.get_velocities()

            # Update the OEMol complex positions to match the new
            # Parmed structure after the simulation
            new_temp_mol = oeommutils.openmmTop_to_oemol(parmed_structure.topology,
                                                         parmed_structure.positions,
                                                         verbose=False)
            new_pos = new_temp_mol.GetCoords()
            system.SetCoords(new_pos)
            record.set_value(Fields.primary_molecule, system)

            # Trajectory
            if opt['trajectory_interval']:
                full_path = os.getcwd()
                filename = os.path.join(full_path, opt['outfname']+'.h5')
                lf = utils.upload_file(filename, opt['outfname']+'.h5')
            else:  # Empty Trajectory
                lf = None

            # Read in logging file if any
            if opt['reporter_interval']:
                with open(opt['outfname'] + '.log', 'r') as flog:
                    str_logger = flog.read()

                opt['str_logger'] += '\n' + str_logger

            md_stage_record = MDRecords.MDStageRecord(MDStageNames.NPT,
                                                      MDRecords.MDSystemRecord(system, new_mdstate),
                                                      log=opt['str_logger'], trajectory=lf)

            if opt['save_md_stage']:
                opt['Logger'].info("[{}] Saving MD stage: {}".format(opt['CubeTitle'], opt['SimType']))

                md_stages.append(md_stage_record)
            else:
                md_stages[-1] = md_stage_record

            record.set_value(Fields.md_stages, md_stages)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return