import os
try:
    import distro
except ImportError:
    import platform as distro
import pickle
import subprocess

import numpy as np

from fem_py.abaqus_material_test.one_element_input_file import write_input_file
from fem_py.common import package_directory
from fem_py.materials.transformation_materials import test_material


def one_element_abaqus(simulation_directory, material, boundary_conditions, simulation_name='oneElement',
                       element_size=1., element_type='C3D8', user_subroutine=None,
                       time_period=1., max_increment=1., martensite_fraction=None, temperature=None, output=True,
                       **_):
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
    job_string = 'j=' + simulation_name + ' interactive'
    if user_subroutine:
        job_string += ' user=' + user_subroutine

    if distro.linux_distribution()[0] == 'Ubuntu':
        abq = '\"singularity exec --nv ' + os.path.expanduser('~/imgs/sing/abaqus-2018-centos-7.img') + \
                   ' vglrun /opt/abaqus/2018/Commands/abq2018\"'
    else:
        abq = '/scratch/users/erik/SIMULIA/CAE/2018/linux_a64/code/bin/ABQLauncher'

    file_lines = ['export LD_PRELOAD=""',
                  'abq=' + abq,
                  '${abq} ' + job_string,
                  '${abq} ' + ' python ' + file_directory + '/one_element_post_processing.py ' + simulation_name]

    with open(simulation_name + '.sh', 'w') as run_file:
        for line in file_lines:
            run_file.write(line + '\n')
    subprocess.Popen('chmod u+x ' + simulation_name + '.sh', shell=True)
    if output is True:
        abaqus_job = subprocess.Popen('./' + simulation_name + '.sh', shell=True)
    else:
        FNULL = open(os.devnull, 'w')
        abaqus_job = subprocess.Popen('./' + simulation_name + '.sh', shell=True, stdout=FNULL,
                                      stderr=subprocess.STDOUT)
    abaqus_job.wait()

    with open('stressStrain' + simulation_name + '.pkl', 'rb') as pickle_handle:
        data = pickle.load(pickle_handle, encoding='latin1')
    stresses = data[:, 0:6]
    strains = data[:, 6:12]
    epl = data[:, 12]
    fM = data[:, 13]
    os.chdir(run_directory)
    return strains, stresses, epl, fM
