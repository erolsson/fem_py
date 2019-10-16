import os
import pickle
from subprocess import Popen

import numpy as np

from common import package_directory
from materials.SS2506 import test_material


def one_element_abaqus(material_parameters, exx=None, eyy=None, ezz=None, pxx=None, pyy=None, pzz=None, time_period=1.,
                       umat_file=None, simulation_name='oneElement', **_):
    run_directory = os.getcwd()
    file_directory = os.path.dirname(os.path.abspath(__file__))
    simulation_directory = file_directory + '/one_element/'
    material = material_parameters

    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)

    os.chdir(simulation_directory)

    with open('parameter_pickle.pkl', 'w') as parameter_pickle:
        pickle.dump((exx, eyy, ezz, pxx, pyy, pzz), parameter_pickle)
        pickle.dump(material, parameter_pickle)
        pickle.dump(time_period, parameter_pickle)
        pickle.dump(umat_file, parameter_pickle)
        pickle.dump(simulation_name)

    # Running the abaqus simulation in an external script
    abaqus_job = Popen('abaqus cae noGUI=' + file_directory + '/one_element.py -- ' + package_directory, shell=True)
    abaqus_job.wait()
    os.chdir(run_directory)

    with open(simulation_directory + '/stressStrain.pkl', 'rb') as pickle_handle:
        data = pickle.load(pickle_handle)
    stresses = data[:, 0:6]
    strains = data[:, 6:12]

    return strains, stresses


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    time = np.linspace(0., 1., 1000)
    strain_z = np.zeros((1000, 2))
    strain_z[:, 0] = time
    strain_z[:, 1] = 0.02*np.sin(2*np.pi*time)
    e, s = one_element_abaqus(material_parameters=test_material, ezz=strain_z)
    plt.plot(e[:, 2], s[:, 2])
    plt.figure()
    plt.plot(e[:, 2], e[:, 0])
    plt.xlabel(r'Strain, $\varepsilon_{xx}$', fontsize=20)
    plt.ylabel(r'Stress, $\sigma_{xx}$', fontsize=20)
    plt.tight_layout()
    plt.show()
