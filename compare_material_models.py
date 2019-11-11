from collections import namedtuple

import os

import numpy as np

import matplotlib.pyplot as plt

from abaqus_material_test.material_test import one_element_abaqus
from one_element_test import one_element_simulation
from abaqus_material_test.one_element_input_file import BC

from materials.SS2506 import test_material
from materials.SS2506 import neu_sehitoglu

Simulation = namedtuple('Simulation', ['model', 'label', 'color', 'umat_file', 'name', 'fm'])
simulations = [Simulation(model=one_element_abaqus, label='New', color='b', umat_file=None, name='oneElementAbaqus',
                          fm=None),
               Simulation(model=one_element_abaqus, label='Abaqus', color='r', name='oneElementUmat',
                          umat_file=os.path.expanduser('~/fem_py/transformation_subroutines/'
                                                       'transformation_subroutine.o'), fm=0.8)]

increments = 100

time = np.linspace(0., 1., increments)
strain_z = np.zeros((increments, 2))
strain_z[:, 0] = time
strain_z[:, 1] = 0.02*np.sin(np.pi*time/2)

pressure_z = np.copy(strain_z)
pressure_z[:, 1] = 3000*np.sin(np.pi*time/2)
temperature = np.array([[0, 22.], [time[-1], 22.]])
boundary_conditions = [BC(amplitude=pressure_z, direction='z', mode='stress')]
args = {'boundary_conditions': boundary_conditions, 'simulation_directory': 'abaqus_material_test/one_element',
        'material': neu_sehitoglu, 'temperature': temperature}

for inc, lw in zip([1.], [1]):
    for simulation in simulations:
        args['user_subroutine'] = simulation.umat_file
        args['simulation_name'] = simulation.name
        args['max_time_increment'] = inc
        args['martensite_fraction'] = simulation.fm
        e, s = simulation.model(**args)
        plt.figure(1)
        plt.plot(e[:, 2], s[:, 2], '-x' + simulation.color, lw=lw, label=simulation.label)
        plt.figure(2)
        plt.plot(e[:, 2], e[:, 1], '-x' + simulation.color, lw=lw, label=simulation.label)
plt.show()
