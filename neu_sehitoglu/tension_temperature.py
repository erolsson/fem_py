import os

import numpy as np

import matplotlib.pyplot as plt

from materials.SS2506 import neu_sehitoglu
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC


def run_sim_from_experiment(name, temperature, stress_strain_data):
    stress_bc = np.array([[0., 0], [1., np.max(stress_strain_data[:, 1])]])
    strain_bc = np.array([[0., 0], [1., np.max(stress_strain_data[:, 0])]])
    e, s, _, _ = one_element_abaqus(simulation_dir, material=neu_sehitoglu,
                                    boundary_conditions=[BC(stress_bc, 'z', 'stress')],
                                    simulation_name=name,
                                    temperature=np.array([[0, temperature], [1, temperature]]),
                                    user_subroutine=umat_file)
    return e, s


umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/fem_py/abaqus_material_test/neu_sehitoglu')

neu_sehitoglu.sde = 0
neu_sehitoglu.sy0M = 500
for temp, color in zip([150., 22., 22.], ['r', 'b', 'g']):
    data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig_4_' + str(int(temp)) + 'C'),
                         delimiter=',')
    plt.plot(data[:, 0], data[:, 1], color + '*')
    sim_name = 'stress_' + str(int(temp)) + '_tens'
    strain, stress = run_sim_from_experiment(sim_name, temp, data)
    plt.plot(strain[:, 2], stress[:, 2], color, lw=2)

data_comp = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig_2_compression'),
                          delimiter=',')
strain, stress = run_sim_from_experiment('stress_22_comp', 22., data_comp)
plt.plot(strain[:, 2], stress[:, 2], 'g', lw=2)
plt.show()
