
from floe.api.orion import in_orion

from Standards.utils import ParmedData

from datarecord import OEPrimaryMolField

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

if in_orion():
    from cuberecord import OELargeFileDataType


# ------------ Stage Standard Names ------------- #

class MDStageNames:
    SETUP = 'SETUP'
    MINIMIZATION = 'MINIMIZATION'
    NVT = 'NVT'
    NPT = 'NPT'
    FEC = 'FEC'

# ---------------- Field Standards -------------- #


class Fields:
    # The Title field is used to set the system name
    title = OEField("Title", Types.String, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ID field should be used as identification number for ligands, proteins or complexes
    id = OEField("ID", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Ligand))

    # The Protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Protein))

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure Field
    structure = OEField('Structure', ParmedData)

    # The Stage Name
    stage_name = OEField('Stage_name', Types.String)

    # Topology Field
    topology = OEField('Topology', Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.PrimaryMol))

    # Log Info
    log_data = OEField('Log_data', Types.String)

    # MD System Field
    md_system = OEField("MDSystem", Types.Record)

    # Trajectory
    if in_orion():
        trajectory = OEField("Trajectory", OELargeFileDataType)
    else:
        trajectory = OEField("Trajectory", Types.String)

    # Stage list Field
    md_stages = OEField("MDStages", Types.RecordVec)

    # Stage Field
    md_stage = OEField("MDStages", Types.Record)


# ---------------- Record Standards -------------- #

# The MDSystemRecord class holds the system topology as an OEMol and the system
# parametrization by using a Parmed Structure object

class MDRecords:
    class MDSystemRecord(OERecord):

        def __init__(self, molecule, structure):
            super().__init__()
            self.set_value(Fields.topology, molecule)
            self.set_value(Fields.structure, structure)

    class MDStageRecord(OERecord):

        def __init__(self, name, log, system_record, trajectory=None):
            super().__init__()
            self.set_value(Fields.stage_name, name)
            self.set_value(Fields.log_data, log)
            self.set_value(Fields.md_system, system_record)
            if not trajectory:
                self.set_value(OEField("Trajectory", Types.String), '')
            else:
                self.set_value(Fields.trajectory, trajectory)
