from collections import namedtuple

import numpy as np

import matplotlib.pyplot as plt

from abaqus_material_test.material_test import one_element_abaqus
from one_element_test import one_element_simulation

from materials.SS2506 import test_material

Simulation = namedtuple('Simulation', ['model', 'label', 'color', 'umat_file'])
simulations = [Simulation(model=one_element_simulation, label='New', color='b', umat_file=None),
               Simulation(model=one_element_abaqus, label='Abaqus', color='r',
                          umat_file='transformation_subroutines/transformation_umat.o')]

increments = 1000

time = np.linspace(0., 1., increments)
strain_z = np.zeros((increments, 2))
strain_z[:, 0] = time
strain_z[:, 1] = 0.02*np.sin(np.pi*time)

pressure_z = np.copy(strain_z)
pressure_z[:, 1] = 3000*np.sin(np.pi*time)

args = {'pzz': pressure_z, 'increments': increments, 'material_parameters': test_material}

for simulation in simulations:
    e, s = simulation.model(**args)
    plt.figure(1)
    plt.plot(e[:, 2], s[:, 2], '-' + simulation.color, lw=2, label=simulation.label)
    plt.figure(2)
    plt.plot(e[:, 2], e[:, 1], '-' + simulation.color, lw=2, label=simulation.label)
plt.show()
