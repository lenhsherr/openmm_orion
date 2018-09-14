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

from yank.experiment import ExperimentBuilder
import os
import fcntl
import time
from floe.api.orion import in_orion
from simtk import openmm


def local_cluster(sim):
    def wrapper(*args):
        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')
            opt = args[0]
            opt['Logger'].info("LOCAL FLOE CLUSTER OPTION IN USE")
            while True:

                gpu_id = gpus_available_indexes[opt['system_id'] % len(gpus_available_indexes)]

                try:
                    with open(gpu_id + '.txt', 'a') as file:
                        fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        file.write("YANK - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(opt['system_title'],
                                                                                                    opt['system_id'],
                                                                                                    gpus_available_indexes,
                                                                                                    gpu_id))

                        openmm.Platform.getPlatformByName('CUDA').setPropertyDefaultValue('DeviceIndex', gpu_id)

                        sim(opt)
                        time.sleep(1.0)
                        fcntl.flock(file, fcntl.LOCK_UN)
                        break
                except IOError:
                    time.sleep(0.01)

        else:
            sim(*args)

    return wrapper


@local_cluster
def run_yank(opt):

    # Print Yank Template
    opt['Logger'].warn(opt['yank_template'])

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return