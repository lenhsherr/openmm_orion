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


from datarecord.types import CustomHandler
import pickle
import parmed

from MDCubes.utils import MDState
import copy


class ParmedData(CustomHandler):

    @staticmethod
    def get_name():
        return 'Parmed'

    @classmethod
    def validate(cls, value):
        return isinstance(value, parmed.structure.Structure)

    @classmethod
    def copy(cls, value):
        return parmed.structure.copy(value)

    @staticmethod
    def serialize(structure):
        struct_dict = structure.__getstate__()
        pkl_obj = pickle.dumps(struct_dict)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_structure = parmed.structure.Structure()
        new_structure.__setstate__(pickle.loads(bytes(data)))
        return new_structure


class MDStateData(CustomHandler):

    @staticmethod
    def get_name():
        return 'MDState'

    @classmethod
    def validate(cls, value):
        return isinstance(value, MDState)

    @classmethod
    def copy(cls, value):
        return copy.deepcopy(value)

    @staticmethod
    def serialize(state):
        pkl_obj = pickle.dumps(state)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_state = pickle.loads(bytes(data))
        return new_state