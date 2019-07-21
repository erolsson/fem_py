from collections import namedtuple

import numpy as np

import matplotlib.pyplot as plt

from abaqus_material_test.material_test import one_element_abaqus
from one_element_test import one_element_simulation

from materials.SS2506 import test_material

Simulation = namedtuple('Simulation', ['model', 'label', 'color'])
simulations = [Simulation(model=one_element_simulation, label='New', color='b'),
               Simulation(model=one_element_abaqus, label='Abaqus', color='r')]

time = np.linspace(0., 1., 1000)
strain_z = np.zeros((1000, 2))
strain_z[:, 0] = time
strain_z[:, 1] = 0.02*np.sin(2*np.pi*time)
args = {'ezz': strain_z, 'increments': 1000, 'material_parameters': test_material}
for simulation in simulations:
    e, s = simulation.model(**args)
    plt.plot(e[:, 2], s[:, 2], simulation.color, lw=2, label=simulation.label)
plt.show()
