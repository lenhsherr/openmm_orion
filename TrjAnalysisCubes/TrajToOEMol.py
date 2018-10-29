import os
import traceback

from floe.api import (ParallelMixin,
                      parameter)

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from cuberecord import OERecordComputeCube

from Standards import Fields

import MDCubes.utils as omm_utils

import TrjAnalysisCubes.utils as utl

import oetrajanalysis.OETrajBasicAnalysis_utils as oetrjutl

import ensemble2img

from tempfile import TemporaryDirectory


class TrajToOEMolCube(ParallelMixin, OERecordComputeCube):
    title = 'Traj to OEMol Cube'

    version = "0.1.0"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Converting MD Traj into multiconf OEMols for Ligand and Protein

    This cube will take in the MD traj file containing
    the solvated protein:ligand complex and extract
    multiconf OEMols for Ligand and Protein.

    Input parameters:
    -------
    oechem.OEDataRecord - Streamed-in MD results for input

    Output:
    -------
    oechem.OEDataRecord - Stream of output data with trajectory OEMols
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process_failed(self, record, port, last_error):
        print("Failed to process record", record, flush=True)
        self.failure.emit(record)

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Logger string
            opt['Logger'].info(' ')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{}: Attempting MD Traj conversion into OEMols'.format(system_title))

            # Extract the MDStageRecord list
            md_stages = utl.RequestOEFieldType(record, Fields.md_stages)
            if len(md_stages) < 2:
                raise ValueError('{} does not have at least 2 MD Stages'.format(system_title))

            # Extract and verify the traj filename for the last MD stage
            md_stageLast_record = md_stages[-1]

            lastName = utl.RequestOEFieldType(md_stageLast_record, Fields.stage_type)

            if lastName != 'NPT':
                raise ValueError('Cannot find the NPT stage')

            with TemporaryDirectory() as output_directory:

                opt['Logger'].info('Temp Directory: {}'.format(output_directory))

                if md_stageLast_record.has_value(Fields.trajectory):
                    trajID = md_stageLast_record.get_value(Fields.trajectory)

                elif md_stageLast_record.has_value(Fields.orion_local_trj_field):
                    opt['Logger'].info("Orion S3 Trajectory Field Detected")
                    trajID = md_stageLast_record.get_value(Fields.orion_local_trj_field)

                else:
                    raise ValueError("No trajectory have been found in the selected stage record {}".format(
                        md_stageLast_record.get_value(Fields.stage_name)))

                trj_selected_filename = os.path.join(output_directory, "trajectory.h5")

                trajName = omm_utils.download_file(trajID, trj_selected_filename, delete=False)

                opt['Logger'].info('{} Trajectory filename: {}'.format(system_title, trajName))

                # Extract the Setup Topology
                md_stage0_record = md_stages[0]

                setupType = utl.RequestOEFieldType(md_stage0_record, Fields.stage_type)

                if setupType != 'SETUP':
                    raise ValueError('Cannot find the SETUP stage')

                md_system = utl.RequestOEFieldType(md_stage0_record, Fields.md_system)

                setupOEMol = utl.RequestOEFieldType(md_system, Fields.topology)

                opt['Logger'].info('Setup topology has {} atoms'.format(setupOEMol.NumAtoms()))

                # Generate multiconformer protein and ligand OEMols from the trajectory
                opt['Logger'].info('{} Generating protein and ligand trajectory OEMols'.format(system_title))

                ptraj, ltraj = utl.ExtractAlignedProtLigTraj_hdf5(setupOEMol, trajName)

                opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'.format(
                    system_title, ptraj.NumAtoms(), ptraj.NumConfs()))
                opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'.format(
                    system_title, ltraj.NumAtoms(), ltraj.NumConfs()))

            # Generate average and median protein and ligand OEMols from ptraj, ltraj
            opt['Logger'].info('{} Generating protein and ligand median and average OEMols'.format(system_title))

            ligMedian, protMedian, ligAverage, protAverage = oetrjutl.AnalyseProteinLigandTrajectoryOEMols(ltraj, ptraj)

            # Generate interactive trajectory SVG
            opt['Logger'].info('{} Generating interactive trajectory SVG'.format(system_title))
            trajSVG = ensemble2img.run_ensemble2img(ligMedian, protMedian, ltraj, ptraj)

            # Overwrite MDStages with only first (setup) and last (production) stages
            newMDStages = [md_stage0_record, md_stageLast_record]

            record.set_value(Fields.md_stages, newMDStages)

            # Create new record with OETraj results
            oetrajRecord = OERecord()

            oetrajRecord.set_value(OEField('ProtTraj', Types.Chem.Mol), ptraj)

            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ltraj)

            oetrajRecord.set_value(OEField('LigMedian', Types.Chem.Mol), ligMedian)

            oetrajRecord.set_value(OEField('ProtMedian', Types.Chem.Mol), protMedian)

            oetrajRecord.set_value(OEField('LigAverage', Types.Chem.Mol), ligAverage)

            oetrajRecord.set_value(OEField('ProtAverage', Types.Chem.Mol), protAverage)

            TrajSVG_field = OEField('TrajSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

            oetrajRecord.set_value(TrajSVG_field, trajSVG)

            record.set_value(OEField('OETraj', Types.Record), oetrajRecord)

            analysesDone = None

            try:
                analysesDone = utl.RequestOEField(record, 'AnalysesDone', Types.StringVec)
                opt['Logger'].info('{}: found AnalysesDone list'.format(system_title))
                analysesDone.append('OETraj')
            except:
                analysesDone = ['OETraj']
                opt['Logger'].info('{}: created AnalysesDone list'.format(system_title))
            if analysesDone is None:
                raise ValueError('{} AnalysesDone list does not exist'.format(system_title))
            record.set_value(OEField('AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{}: saved protein and ligand traj OEMols'.format(system_title))

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.log.info('Exception in TrajToOEMolCube on {}'.format(system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

