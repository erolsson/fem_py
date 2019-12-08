import os
import pickle
from subprocess import Popen

import numpy as np

from one_element_input_file import write_input_file
from common import package_directory
from materials.transformation_materials import test_material


def one_element_abaqus(simulation_directory, material, boundary_conditions, simulation_name='oneElement',
                       element_size=1., element_type='C3D8R', user_subroutine=None,
                       time_period=1., max_increment=1., martensite_fraction=None, temperature=None, **_):
    run_directory = os.getcwd()
    file_directory = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)

    os.chdir(simulation_directory)

    write_input_file(filename=simulation_name + '.inp', material=material, boundary_conditions=boundary_conditions,
                     element_size=element_size, element_type=element_type, umat_file=user_subroutine,
                     time_period=time_period, max_increment=max_increment, martensite_fraction=martensite_fraction,
                     temperature=temperature)

    with open('abaqus_v6.env', 'w') as env_file:
        env_file.write('ask_delete = OFF\n')

    # Running the abaqus simulation in an external script
    job_string = 'abaqus j=' + simulation_name + ' interactive'
    if user_subroutine:
        job_string += ' user=' + user_subroutine
    abaqus_job = Popen(job_string, shell=True)
    abaqus_job.wait()
    os.chdir(run_directory)
    odb_path = os.path.abspath(simulation_directory + '/' + simulation_name + '.odb')
    odb_abs_dir = os.path.dirname(odb_path)
    os.chdir(file_directory)
    abaqus_post_processing_job = Popen('abaqus python one_element_post_processing.py ' +
                                       odb_abs_dir + ' ' + simulation_name, shell=True)

    abaqus_post_processing_job.wait()

    os.chdir(run_directory)
    with open(simulation_directory + '/stressStrain' + simulation_name + '.pkl', 'rb') as pickle_handle:
        data = pickle.load(pickle_handle)
    stresses = data[:, 0:6]
    strains = data[:, 6:12]
    epl = data[:, 12]
    fM = data[:, 13]
    return strains, stresses, epl, fM
