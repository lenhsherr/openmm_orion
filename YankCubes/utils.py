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

from YankCubes.yank_templates import (max_cube_running_time,
                                      yank_binding_template,
                                      yank_solvation_template)

from yank.analyze import ExperimentAnalyzer

from simtk.openmm import unit

from datetime import timedelta

import os

import subprocess

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

import yaml


def yank_solvation_initialize(sim):
    def wrapper(*args):

        opt = args[0]

        yank_template = yank_solvation_template.format(
            verbose='yes' if opt['verbose'] else 'no',
            minimize='yes' if opt['minimize'] else 'no',
            output_directory=opt['output_directory'],
            timestep=4.0 if opt['hmr'] else 2.0,
            nsteps_per_iteration=opt['nsteps_per_iteration'],
            number_iterations=opt['new_iterations'],
            temperature=opt['temperature'],
            pressure=opt['pressure'],
            resume_sim='yes' if opt['resume_sim'] else 'no',
            resume_setup='yes' if opt['resume_setup'] else 'no',
            hydrogen_mass=4.0 if opt['hmr'] else 1.0,
            alchemical_pme_treatment=opt['alchemical_pme_treatment'],
            checkpoint_interval=opt['checkpoint_interval'],
            solvated_pdb_fn=opt['solvated_structure_fn'],
            solvated_xml_fn=opt['solvated_omm_sys_serialized_fn'],
            solute_pdb_fn=opt['solute_structure_fn'],
            solute_xml_fn=opt['solute_omm_sys_serialized_fn'],
            solvent_dsl=opt['solvent_str_names'])

        opt['yank_template'] = yank_template

        # Print Yank Template
        opt['Logger'].info(opt['yank_template'])

        sim(opt)

    return wrapper


def yank_binding_initialize(sim):
    def wrapper(*args):

        opt = args[0]

        yank_template = yank_binding_template.format(
            verbose='yes' if opt['verbose'] else 'no',
            minimize='yes' if opt['minimize'] else 'no',
            output_directory=opt['output_directory'],
            timestep=opt['timestep'],
            nsteps_per_iteration=opt['nsteps_per_iteration'],
            number_iterations=opt['new_iterations'],
            temperature=opt['temperature'],
            pressure=opt['pressure'],
            resume_sim='yes' if opt['resume_sim'] else 'no',
            resume_setup='yes' if opt['resume_setup'] else 'no',
            hydrogen_mass=4.0 if opt['hmr'] else 1.0,
            alchemical_pme_treatment=opt['alchemical_pme_treatment'],
            checkpoint_interval=opt['checkpoint_interval'],
            complex_pdb_fn=opt['solvated_complex_structure_fn'],
            complex_xml_fn=opt['solvated_complex_omm_serialized_fn'],
            solvent_pdb_fn=opt['solvated_ligand_structure_fn'],
            solvent_xml_fn=opt['solvated_ligand_omm_serialized_fn'],
            ligand_resname=opt['lig_res_name'],
            solvent_dsl=opt['solvent_str_names'],
            sampler=opt['sampler'],
            restraints=opt['restraints'],
            protocol=opt['protocol'])

        if opt['user_yank_yaml_file'] is not None:

            opt['Logger'].info("Yank User template file in use ")

            # Loading the yank_template as yaml
            yank_yaml_template = yaml.load(yank_template)

            # Loading the user Yaml file
            try:
                fn = opt['user_yank_yaml_file']
                with open(fn, 'r') as yaml_file:
                    yank_yaml_user = yaml.load(yaml_file)
            except:
                raise IOError("It was not possible to load the provided yaml file")

            yank_yaml_merge = dict(yank_yaml_template)

            # Cleaning
            del yank_yaml_merge['harmonic']
            del yank_yaml_merge['boresch']

            if 'experiments' in yank_yaml_user:
                if isinstance(yank_yaml_user['experiments'], list):

                    for exp in yank_yaml_user['experiments']:
                        yank_yaml_merge[exp] = yank_yaml_user[exp]

                    yank_yaml_merge['experiments'] = yank_yaml_user['experiments']

                    # yank_template = yank_yaml_merge
                else:
                    yank_yaml_merge['experiments'] = ["exp"]

                    yank_yaml_merge['exp'] = yank_yaml_user['experiments']

                if len(yank_yaml_merge['experiments']) > 1:
                    raise ValueError("Currently just one Yank experiment is supported. Provided: {}".format(yank_yaml_merge['experiments']))

                for exp in yank_yaml_merge['experiments']:

                    # System section
                    yank_yaml_merge[exp]['system'] = 'system'

                    # Protocol Section
                    if 'protocol' in yank_yaml_merge[exp]:
                        yank_yaml_merge['protocols'][yank_yaml_merge[exp]['protocol']] = yank_yaml_user['protocols'][yank_yaml_merge[exp]['protocol']]

                    else:  # If the protocol section is not present use the selected one
                        yank_yaml_merge[exp]['protocol'] = yank_yaml_template[yank_yaml_template['experiments'][0]]['protocol']

                    if 'sampler' in yank_yaml_merge[exp]:
                        yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']] = yank_yaml_user['samplers'][yank_yaml_merge[exp]['sampler']]
                        if 'mcmc_moves' in yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]:
                            yank_yaml_merge['mcmc_moves'][yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves']] = yank_yaml_user['mcmc_moves'][yank_yaml_user['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves']]
                        else:  # mcmc move is set to the default one
                            yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves'] = yank_yaml_merge['samplers'][yank_yaml_template[yank_yaml_template['experiments'][0]]['sampler']]['mcmc_moves']
                    else:  # The sampler is set to the default one
                        yank_yaml_merge[exp]['sampler'] = yank_yaml_template[yank_yaml_template['experiments'][0]]['sampler']

                    # Sampler section number of iterations
                    yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['number_of_iterations'] = opt['new_iterations']

                yank_template = yank_yaml_merge
            else:
                raise ValueError("Missing Experiments section in the provided Yank Yaml file")

        opt['yank_template'] = yank_template

        # Print Yank Template
        if isinstance(opt['yank_template'], dict):
            opt['Logger'].info(yaml.dump(opt['yank_template']))
        else:
            opt['Logger'].info(opt['yank_template'])

        sim(opt)

    return wrapper


@yank_solvation_initialize
def run_yank_solvation(opt):

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


@yank_binding_initialize
def run_yank_binding(opt):

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


def run_yank_analysis(opt):

    exp_dir = os.path.join(opt['output_directory'], "experiments")

    experiment_to_analyze = ExperimentAnalyzer(exp_dir)
    analysis = experiment_to_analyze.auto_analyze()

    # Calculate free energy and its error
    DeltaG = analysis['free_energy']['free_energy_diff_unit'].\
                 in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
    dDeltaG = analysis['free_energy']['free_energy_diff_error_unit'].\
                  in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole

    opt_1 = '--store={}'.format(exp_dir)

    result_fn = os.path.join(opt['output_directory'], 'results.html')
    opt_2 = '--output={}'.format(result_fn)

    opt_3 = '--format=html'

    try:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        subprocess.check_call(['yank', 'analyze', 'report', opt_1, opt_2, opt_3])

        # Upload Floe Report
        if in_orion():
            session = OrionSession()

            file_upload = File.upload(session,
                                      "{}.html".format(opt['system_title']),
                                      result_fn)

            session.tag_resource(file_upload, "floe_report")

            job_id = environ.get('ORION_JOB_ID')

            if job_id:
                session.tag_resource(file_upload, "Job {}".format(job_id))

    except subprocess.SubprocessError:
        opt['Logger'].warn("The result file have not been generated")

    return DeltaG, dDeltaG


def calculate_iteration_time(output_directory, num_iterations):

    UNITS = {"s": "seconds", "m": "minutes", "h": "hours", "d": "days", "w": "weeks"}

    def convert_to_seconds(s):
        count = int(float(s[:-1]))
        user_unit = UNITS[s[-1]]
        td = timedelta(**{user_unit: count})
        return td.seconds + 60 * 60 * 24 * td.days

    log_fn = os.path.join(output_directory, "experiments/experiments.log")

    f = open(log_fn, 'r')

    log = f.read()

    times = []
    for line in log.splitlines():
        if "Iteration took" in line:
            times.append(line)

    clean_times = [i.split('-')[-1].split(" ")[-1][:-1] for i in times]

    # times in seconds
    tms_list = []

    for tm in clean_times:
        tms = convert_to_seconds(tm)
        tms_list.append(tms)

    # average time per iterations in hrs
    avg_time_per_iteration_phase_1 = (sum(tms_list[:num_iterations]) / num_iterations) / 3600.0
    avg_time_per_iteration_phase_2 = (sum(tms_list[num_iterations:]) / num_iterations) / 3600.0

    avg_time_per_iteration = avg_time_per_iteration_phase_1 + avg_time_per_iteration_phase_2

    iterations_per_cube = int(max_cube_running_time / avg_time_per_iteration)

    return iterations_per_cube
