from collections import namedtuple
from copy import deepcopy

import os

import numpy as np

import matplotlib.pyplot as plt

from abaqus_material_test.material_test import one_element_abaqus
from one_element_test import one_element_simulation
from abaqus_material_test.one_element_input_file import BC

from materials.SS2506 import test_material
from materials.SS2506 import neu_sehitoglu

Simulation = namedtuple('Simulation', ['model', 'label', 'color', 'umat_file', 'name', 'fm'])
simulations = [  # Simulation(model=one_element_abaqus, label='New', color='b', umat_file=None, name='oneElementAbaqus',
                 #            fm=None),
               Simulation(model=one_element_abaqus, label='Abaqus', color='r', name='oneElementUmat',
                          umat_file=os.path.expanduser('~/fem_py/transformation_subroutines/'
                                                       'transformation_subroutine.o'), fm=0.8)]

increments = 100
neu_sehitoglu2 = deepcopy(neu_sehitoglu)
neu_sehitoglu2.beta = 0
neu_sehitoglu.a1 = 0
neu_sehitoglu.a2 = 0

time = np.linspace(0., 1., increments)
strain_z = np.zeros((increments, 2))
strain_z[:, 0] = time
strain_z[:, 1] = 1.*np.sin(np.pi*time/2)

pressure_z = np.copy(strain_z)
pressure_z[:, 1] = 3000*np.sin(np.pi*time/2)
temperature = np.array([[0, 22.], [time[-1], 22.]])

boundary_conditions = [BC(amplitude=strain_z, direction='z', mode='strain')]
args = {'boundary_conditions': boundary_conditions, 'simulation_directory': 'abaqus_material_test/one_element',
        'temperature': temperature, 'max_increment': 1.}

for mat, lw in zip([neu_sehitoglu, neu_sehitoglu2], [1, 2]):
    for simulation in simulations:
        args['user_subroutine'] = simulation.umat_file
        args['simulation_name'] = simulation.name
        args['max_time_increment'] = 1.
        args['martensite_fraction'] = simulation.fm
        args['material'] = mat
        e, s, epl, fM = simulation.model(**args)
        plt.figure(1)
        plt.plot(e[:, 2], s[:, 2], '-x' + simulation.color, lw=lw, label=simulation.label)
        plt.figure(2)
        plt.plot(e[:, 2], epl, '-x' + simulation.color, lw=lw, label=simulation.label)
        plt.figure(3)
        plt.plot(e[:, 2], fM, '-x' + simulation.color, lw=lw, label=simulation.label)
plt.show()
