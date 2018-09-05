#!/usr/bin/env python

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

from floe.api import WorkFloe
from cuberecord import (DatasetWriterCube,
                        DatasetReaderCube)
from ComplexPrepCubes.cubes import SolvationCube
from ForceFieldCubes.cubes import ForceFieldCube
from LigPrepCubes.cubes import (LigandChargeCube,
                                LigandSetting)

from YankCubes.cubes import YankSolvationFECube
from MDCubes.OpenMMCubes.cubes import (OpenMMminimizeCube,
                                       OpenMMNvtCube,
                                       OpenMMNptCube)

from YankCubes.yank_templates import number_cubes_solvation

job = WorkFloe("Solvation Free Energy")

job.description = """
    The Solvation Free Energy protocol performs Solvation Free Energy Calculations (SFEC) on 
    a set of provided ligands by using YANK ( http://getyank.org/latest/ ). The ligands need 
    to have coordinates, correct chemistry and must be neutral. Each ligand can have multiple conformers, 
    but each conformer will be treated as a different ligand and prepared to run SFEC. The ligands 
    are solvated in the selected mixture (default water) and parametrized accordingly to the 
    provided force field. A minimization stage is performed on the system followed by a warm up 
    (NVT ensemble) and an equilibration stage (NPT ensemble). In the minimization, warm up 
    and equilibration stage positional harmonic restraints are applied to the ligands with different 
    force constants. At the end of the equilibration stage the SFEC calculation is run by YANK 
    with the selected parameters. Solvation Free Energy values and floe reports are output.    
    
    Required Input Parameters:
    -----------
    ligands: Dataset of the prepared ligands
    
    Outputs:
    --------
    out : Dataset of the solvated systems with the calculated solvation free energies and 
          floe reports
"""

cube_list = []

job.classification = [['Solvation Free Energy']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands",
                        title="Ligand Input File",
                        description="Ligand file name")
job.add_cube(iligs)
cube_list.append(iligs)

chargelig = LigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)

job.add_cube(chargelig)
cube_list.append(chargelig)

ligset = LigandSetting("LigandSetting")
job.add_cube(ligset)
cube_list.append(ligset)

solvate = SolvationCube("Solvation", title="System Solvation")
solvate.promote_parameter("density", promoted_name="density", title="Solution density in g/ml", default=1.0,
                          description="Solution Density in g/ml")
solvate.promote_parameter("solvents", promoted_name="solvents", title="Solvent components",
                          default='[H]O[H]',
                          description="Comma separated smiles strings of solvent components")
solvate.promote_parameter("molar_fractions", promoted_name="molar_fractions",
                          title="Molar fractions",
                          default='1.0',
                          description="Comma separated strings of solvent molar fractions")
solvate.set_parameters(distance_between_atoms=2.5)
solvate.set_parameters(padding_distance=11.0)

job.add_cube(solvate)
cube_list.append(solvate)

ff = ForceFieldCube("ForceField", title="System Parametrization")
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_forcefield', default='GAFF2')
job.add_cube(ff)
cube_list.append(ff)


# First Yank Cube used to build the UI interface
solvationfe0 = YankSolvationFECube("SovationFE0", title="Solvation FE0")
solvationfe0.promote_parameter('iterations', promoted_name='iterations', default=1000,
                               description="Total Number of Yank iterations for the entire floe. "
                                           "A Yank iteration is 500 MD steps")
solvationfe0.promote_parameter('verbose', promoted_name='verbose', default=True)
solvationfe0.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                               description='Temperature (Kelvin)')
solvationfe0.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                               description='Pressure (atm)')
solvationfe0.promote_parameter('hmr', promoted_name='hmr', default=False,
                               description='Hydrogen Mass Repartitioning')

job.add_cube(solvationfe0)

# Minimization
minimize = OpenMMminimizeCube("Minimize", title="System Minimization")
minimize.set_parameters(restraints='noh ligand')
minimize.set_parameters(restraintWt=5.0)
minimize.set_parameters(center=True)
minimize.promote_parameter("hmr", promoted_name="hmr")

job.add_cube(minimize)
cube_list.append(minimize)

# NVT Warm-up
warmup = OpenMMNvtCube('warmup', title='System Warm Up')
warmup.set_parameters(time=0.02)
warmup.promote_parameter("temperature", promoted_name="temperature")
warmup.set_parameters(restraints="noh ligand")
warmup.set_parameters(restraintWt=2.0)
warmup.set_parameters(trajectory_interval=0.0)
warmup.set_parameters(reporter_interval=0.0)
warmup.set_parameters(suffix='warmup')
warmup.promote_parameter("hmr", promoted_name="hmr")

job.add_cube(warmup)
cube_list.append(warmup)

# NPT Equilibration stage
equil = OpenMMNptCube('equil', title='System Equilibration')
equil.set_parameters(time=0.02)
equil.promote_parameter("temperature", promoted_name="temperature")
equil.promote_parameter("pressure", promoted_name="pressure")
equil.promote_parameter("hmr", promoted_name="hmr")
equil.set_parameters(restraints="noh ligand")
equil.set_parameters(restraintWt=0.1)
equil.set_parameters(trajectory_interval=0.0)
equil.set_parameters(reporter_interval=0.0)
equil.set_parameters(suffix='equil')

job.add_cube(equil)
cube_list.append(equil)

# Add YANK first Cube
cube_list.append(solvationfe0)

for i in range(1, number_cubes_solvation):
    solvationfe = YankSolvationFECube("SovationFE"+str(i), title="Solvation FE"+str(i))
    solvationfe.promote_parameter('iterations', promoted_name='iterations')
    solvationfe.promote_parameter("verbose", promoted_name="verbose")
    solvationfe.promote_parameter("temperature", promoted_name="temperature")
    solvationfe.promote_parameter("pressure", promoted_name="pressure")
    solvationfe.promote_parameter("hmr", promoted_name="hmr")

    job.add_cube(solvationfe)
    cube_list.append(solvationfe)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out", description="Data Set Out Name")
job.add_cube(ofs)
cube_list.append(ofs)

fail = DatasetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')
job.add_cube(fail)
cube_list.append(fail)

# Connections
for i in range(0, len(cube_list)-2):
    cube_list[i].success.connect(cube_list[i + 1].intake)
    if i == len(cube_list) - 3:
        cube_list[i].failure.connect(cube_list[i+2].intake)

if __name__ == "__main__":
    job.run()
