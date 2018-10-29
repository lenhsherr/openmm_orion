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
from datarecord import OERecord
from cuberecord import OERecordComputeCube
from cuberecord.ports import RecordInputPort

from floe.api import (ParallelMixin,
                      parameter)
from openeye import oechem
from oeommtools import (utils as oeommutils,
                        packmol)


from Standards import Fields

from floe.constants import ADVANCED


class SolvationCube(ParallelMixin, OERecordComputeCube):
    title = "Solvation Cube Packmol"
    version = "0.0.0"
    classification = [["Preparation", "OEChem"]]
    tags = ['OEChem', 'PackMol']
    description = """
    This cube solvate a molecular system

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEDataRecord - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=8.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        default='[H]O[H]',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text='Output Packmol log')

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")

    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = self.opt

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing the Primary Molecule Field")

            solute = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Title field")
                solute_title = solute.GetTitle()[0:12]
            else:
                solute_title = record.get_value(Fields.title)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                            ["solvents", "molar_fractions", "density"]}
            if new_args:
                for k in new_args:
                    if k == 'molar_fractions':
                        continue
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute_title, new_args))
                opt.update(new_args)

            # Solvate the system
            sol_system = packmol.oesolvate(solute, **opt)
            self.log.info("[{}] Solvated System atom number: {}".format(self.title,
                                                                        sol_system.NumAtoms()))
            sol_system.SetTitle(solute.GetTitle())

            record.set_value(Fields.primary_molecule, sol_system)
            record.set_value(Fields.title, solute_title)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class ComplexPrepCube(OERecordComputeCube):
    title = "Complex Preparation Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem']
    description = """
        This cube assembles the complex made of the protein and the
        ligands. If a ligand presents multiple conformers, then each conformer
        is bonded to the protein to form the solvated complex. For example if a
        ligand has 3 conformers then 3 complexes are generated.

        Input:
        -------
        oechem.OEDataRecord - Streamed-in of the protein and ligands

        Output:
        -------
        oechem.OEDataRecord - Emits the complex molecules
        """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    # Ligand Residue Name
    lig_res_name = parameter.StringParameter('lig_res_name',
                                             default='LIG',
                                             help_text='The ligand residue name',
                                             level=ADVANCED)

    protein_port = RecordInputPort("protein_port", initializer=True)

    def begin(self):
        for record in self.protein_port:
            self.opt = vars(self.args)
            self.opt['Logger'] = self.log

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing Protein field")
                self.failure.emit(record)
                return

            protein = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Protein Title field")
                self.protein_title = protein.GetTitle()[0:12]
            else:
                self.protein_title = record.get_value(Fields.title)

            self.protein = protein
            return

    def process(self, record, port):
        try:
            if port == 'intake':

                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Missing the Primary Molecule field")

                ligand = record.get_value(Fields.primary_molecule)

                if ligand.NumConfs() > 1:
                    raise ValueError("The ligand {} has multiple conformers: {}".format(ligand.GetTitle(),
                                                                                        ligand.GetNumConfs()))

                if not record.has_value(Fields.title):
                    self.log.warn("Missing Ligand record '{}' field".format(Fields.title.get_name()))
                    ligand_title = ligand.GetTitle()[0:12]
                else:
                    ligand_title = record.get_value(Fields.title)

                if not record.has_value(Fields.id):
                    raise ValueError("Missing Ligand record '{}' field".format(Fields.id.get_name()))

                ligand_id = record.get_value(Fields.id)

                complx = self.protein.CreateCopy()
                oechem.OEAddMols(complx, ligand)

                # Split the complex in components
                protein_split, ligand_split, water, excipients = oeommutils.split(complx,
                                                                                  ligand_res_name=self.opt['lig_res_name'])

                # If the protein does not contain any atom emit a failure
                if not protein_split.NumAtoms():  # Error: protein molecule is empty
                    raise ValueError("The protein molecule does not contains atoms")

                # If the ligand does not contain any atom emit a failure
                if not ligand_split.NumAtoms():  # Error: ligand molecule is empty
                    raise ValueError("The Ligand molecule does not contains atoms")

                # Check if the ligand is inside the binding site. Cutoff distance 3A
                if not oeommutils.check_shell(ligand_split, protein_split, 3):
                    raise ValueError("The ligand is probably outside the protein binding site")

                # Removing possible clashes between the ligand and water or excipients
                if water.NumAtoms():
                    water_del = oeommutils.delete_shell(ligand, water, 1.5, in_out='in')

                if excipients.NumAtoms():
                    excipient_del = oeommutils.delete_shell(ligand, excipients, 1.5, in_out='in')

                # Reassemble the complex
                new_complex = protein_split.CreateCopy()
                oechem.OEAddMols(new_complex, ligand_split)
                if excipients.NumAtoms():
                    oechem.OEAddMols(new_complex, excipient_del)
                if water.NumAtoms():
                    oechem.OEAddMols(new_complex, water_del)

                complex_title = 'p' + self.protein_title + '_' + ligand_title

                new_complex.SetTitle(ligand.GetTitle())

                new_record = OERecord()

                new_record.set_value(Fields.primary_molecule, new_complex)
                new_record.set_value(Fields.title, complex_title)
                new_record.set_value(Fields.ligand, ligand)
                new_record.set_value(Fields.protein, self.protein)
                new_record.set_value(Fields.id, ligand_id)

                self.success.emit(new_record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return