import os
import numpy as np
from simtk import unit
import yaml
from yank.analyze import get_analyzer
from openeye import oechem
import tarfile
from floe.api.orion import in_orion, stream_file, upload_file


def analyze_directory(source_directory):
    """
    This Function has been copied and adapted from the Yank ver 0.17.0 source code
    (yank.analyse.analyze_directory)

    Analyze contents of store files to compute free energy differences.

    This function is needed to preserve the old auto-analysis style of YANK. What it exactly does can be refined
    when more analyzers and simulations are made available. For now this function exposes the API.

    Parameters
    ----------
    source_directory : string
        The location of the simulation storage files.

    """
    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    if not os.path.isfile(analysis_script_path):
        err_msg = 'Cannot find analysis.yaml script in {}'.format(source_directory)
        raise RuntimeError(err_msg)
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)

    data = dict()
    for phase_name, sign in analysis:
        phase_path = os.path.join(source_directory, phase_name + '.nc')
        phase = get_analyzer(phase_path)
        data[phase_name] = phase.analyze_phase()
        kT = phase.kT

    # Compute free energy and enthalpy
    DeltaF = 0.0
    dDeltaF = 0.0
    DeltaH = 0.0
    dDeltaH = 0.0
    for phase_name, sign in analysis:
        DeltaF -= sign * (data[phase_name]['DeltaF'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaF += data[phase_name]['dDeltaF'] ** 2
        DeltaH -= sign * (data[phase_name]['DeltaH'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaH += data[phase_name]['dDeltaH'] ** 2
    dDeltaF = np.sqrt(dDeltaF)
    dDeltaH = np.sqrt(dDeltaH)

    DeltaF = DeltaF * kT / unit.kilocalories_per_mole
    dDeltaF = dDeltaF * kT / unit.kilocalories_per_mole
    DeltaH = DeltaH * kT / unit.kilocalories_per_mole
    dDeltaH = dDeltaH * kT / unit.kilocalories_per_mole

    return DeltaF, dDeltaF, DeltaH, dDeltaH


def download(molecule, path):
    file_id = oechem.OEGetSDData(molecule, 'file_id')

    if in_orion():
        tar_file = os.path.basename(path) + '.tar.gz'
        with open(tar_file, mode='w') as archive:
            for chunk in stream_file(file_id):
                archive.write(chunk)
    else:
        tar_file = file_id

    tar = tarfile.open(tar_file)
    tar.extractall(path=path)
    tar.close()

    return


def upload(molecule, path):
    file_id = os.path.basename(path) + '.tar.gz'

    with tarfile.open(file_id, mode='w:gz') as archive:
        archive.add(path, arcname='.', recursive=True)

    if in_orion():  # TODO in Orion the restarting is going to fail
        # file_json = upload_file(file_id, file_id, tags="YANK_TMP")
        # file_id = file_json["id"]
        pass
    else:
        oechem.OESetSDData(molecule, 'file_id', file_id)

    return 