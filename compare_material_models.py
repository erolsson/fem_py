from collections import namedtuple

import os

import numpy as np

import matplotlib.pyplot as plt

from abaqus_material_test.material_test import one_element_abaqus
from one_element_test import one_element_simulation

from materials.SS2506 import test_material

Simulation = namedtuple('Simulation', ['model', 'label', 'color', 'umat_file', 'name'])
simulations = [Simulation(model=one_element_abaqus, label='New', color='b', umat_file=None, name='oneElementAbaqus'),
               Simulation(model=one_element_abaqus, label='Abaqus', color='r', name='oneElementUmat',
                          umat_file=os.path.expanduser('~/fem_py/transformation_subroutines/'
                                                       'transformation_subroutine.o'))]

increments = 1000

time = np.linspace(0., 1., increments)
strain_z = np.zeros((increments, 2))
strain_z[:, 0] = time
strain_z[:, 1] = 0.02*np.sin(np.pi*time)

pressure_z = np.copy(strain_z)
pressure_z[:, 1] = 3000*np.sin(np.pi*time)

args = {'pzz': pressure_z, 'increments': increments, 'material_parameters': test_material}

# for inc, lw in zip([1., 1e-1, 1e-2, 1e-3], [1, 2, 3, 4]):
for inc, lw in zip([1.], [1]):
    for simulation in simulations:
        args['umat_file'] = simulation.umat_file
        args['simulation_name'] = simulation.name
        args['max_time_increment'] = inc
        e, s = simulation.model(**args)
        plt.figure(1)
        plt.plot(e[:, 2], s[:, 2], '-x' + simulation.color, lw=lw, label=simulation.label)
        plt.figure(2)
        plt.plot(e[:, 2], e[:, 1], '-x' + simulation.color, lw=lw, label=simulation.label)
        plt.figure(3)
        plt.plot(np.diff(s[:, 2]), '-x' + simulation.color, lw=lw, label=simulation.label)
plt.show()
