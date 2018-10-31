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

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

from abc import ABC, abstractmethod

import parmed

import simtk


class MDState(object):

    def __init__(self, parmed_structure):

        if not parmed_structure.positions:
            raise RuntimeError('Atom positions are not defined')
        else:
            # The returned object is an OpenMM Quantity with units
            self.__positions__ = parmed_structure.positions

        if parmed_structure.velocities is None:
            self.__velocities__ = None
        else:
            # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
            self.__velocities__ = parmed_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond
            # The returned object is an OpenMM Quantity with units

        if parmed_structure.box_vectors is None:
            self.__box_vectors__ = None
        else:
            self.__box_vectors__ = parmed_structure.box_vectors

    def get_positions(self):
        return self.__positions__

    def get_velocities(self):
        return self.__velocities__

    def get_box_vectors(self):
        return self.__box_vectors__

    def set_positions(self, positions):
        if isinstance(positions, simtk.unit.quantity.Quantity):
            self.__positions__ = positions
        else:
            raise ValueError("It was not possible to set the positions")

        return

    def set_velocities(self, velocities):
        if isinstance(velocities, simtk.unit.quantity.Quantity):
            self.__velocities__ = velocities
        else:
            raise ValueError("It was not possible to set the velocities")

        return

    def set_box_vectors(self, box_vectors):
        if isinstance(box_vectors, simtk.unit.quantity.Quantity):
            self.__box_vectors__ = box_vectors
        else:
            raise ValueError("It was not possible to set the box vectors")

        return


class MDSimulations(ABC):
    @abstractmethod
    def __init__(self, mdstate, ff_parameters, opt):

        if not isinstance(mdstate, MDState):
            raise ValueError("{} is not a MDState Object".format(type(mdstate)))

        if not isinstance(ff_parameters, parmed.Structure):
            raise ValueError("{} is not a Parmed Structure Object".format(type(mdstate)))

        if not isinstance(opt, dict):
            raise ValueError("{} is not a dictionary".format(type(opt)))

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def update_state(self):
        pass


def md_simulation(mdstate, ff_parameters, opt):

    if opt['md_engine'] == 'OpenMM':

        from MDCubes.OpenMMCubes.simtools import OpenMMSimulations

        MDSim = OpenMMSimulations(mdstate, ff_parameters, opt)

        MDSim.run()

        new_mdstate = MDSim.update_state()

        return new_mdstate
    else:
        raise ValueError("The selected MD engine is not currently supported: {}".format(opt['md_engine']))


def upload_file(filename, orion_name):

    if in_orion():

        session = OrionSession()

        file_upload = File.upload(session,
                                  orion_name,
                                  filename)

        session.tag_resource(file_upload, "Trajectory")

        job_id = environ.get('ORION_JOB_ID')

        if job_id:
            session.tag_resource(file_upload, "Job {}".format(job_id))

        file_id = file_upload.id

    else:
        file_id = filename

    return file_id


def download_file(file_id, filename, delete=False):

    if in_orion() or isinstance(file_id, int):

        session = OrionSession()

        resource = session.get_resource(File, file_id)

        resource.download_to_file(filename)

        fn_local = filename

        if delete:
            session.delete_resource(resource)
    else:
        fn_local = file_id

    return fn_local


# class MDStructure(object):
#     def __init__(self, parmed_structure, mdengine='OpenMM'):
#
#         self.mdengine = mdengine
#
#         if not parmed_structure.positions:
#             raise RuntimeError('Atom positions are not defined')
#         else:
#             # The returned object is an OpenMM Quantity with units
#             self.__positions__ = parmed_structure.positions
#
#         if parmed_structure.velocities is None:
#             self.__velocities__ = None
#         else:
#             # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
#             self.__velocities__ = parmed_structure.velocities * unit.angstrom/unit.picosecond
#             # The returned object is an OpenMM Quantity with units
#
#         if parmed_structure.box_vectors is None:
#             self.__box_vectors__ = None
#         else:
#             self.__box_vectors__ = parmed_structure.box_vectors
#
#         if self.__box_vectors__:
#             omm_system = parmed_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
#                                                        nonbondedCutoff=10.0 * unit.angstroms,
#                                                        constraints=app.HBonds,
#                                                        removeCMMotion=False)
#         else:
#             omm_system = parmed_structure.createSystem(nonbondedMethod=app.NoCutoff,
#                                                        constraints=app.HBonds,
#                                                        removeCMMotion=False)
#         with TemporaryDirectory() as output_directory:
#
#             # Serialize ethe system bu using OpenMM
#             omm_sys_serialized = XmlSerializer.serialize(omm_system)
#             omm_sys_serialized_fn = os.path.join(output_directory, "system.xml")
#
#             # Write out the Serialized file
#             with open(omm_sys_serialized_fn, 'w') as system_f:
#                 system_f.write(omm_sys_serialized)
#
#             with open(omm_sys_serialized_fn, 'r') as system_f:
#                 self.forcefieldstring = system_f.read()
#
#     def get_state(self):
#         if self.mdengine == 'OpenMM':
#             return MDState(self.__positions__, self.__velocities__, self.__box_vectors__)
#         else:
#             raise ValueError("The selected MD Engine is not currently supported: {}".format(self.mdengine))
#
#     def set_state(self, positions, velocities, box, mdengine='OpenMM'):
#         if mdengine == 'OpenMM':
#             self.__positions__ = positions
#             self.__velocities__ = velocities
#             self.__box_vectors__ = box
#         else:
#             raise ValueError("The selected MD Engine is not currently supported: {}".format(mdengine))

