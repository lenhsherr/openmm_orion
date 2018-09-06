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

import os
from orionclient.session import OrionSession
from artemis.wrappers import WorkFloeWrapper, DatasetWrapper, OutputDatasetWrapper
from artemis.test import FloeTestCase
from artemis.decorators import package

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

session = OrionSession()


@package(PACKAGE_DIR)
class TestYankBindingFloes(FloeTestCase):

    def test_yank_binding_repex_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "Binding_free_energy_repex.py"),
            run_timeout=8000,
            queue_timeout=1200
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "toluene.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "lysozyme.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "iterations": 5,
                    "out": output_file.identifier
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

    def test_yank_binding_sams_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "Binding_free_energy_sams.py"),
            run_timeout=8000,
            queue_timeout=1200
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "toluene.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "lysozyme.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "iterations": 5,
                    "out": output_file.identifier
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)