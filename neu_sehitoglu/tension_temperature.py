import os

import numpy as np

import matplotlib.pyplot as plt

from materials.SS2506 import neu_sehitoglu
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC


def run_sim_from_experiment(name, temperature, stress_strain_data, sign=1):
    stress_bc = np.array([[0., 0], [1., np.max(stress_strain_data[:, 1])*sign]])
    strain_bc = np.array([[0., 0], [1., np.max(stress_strain_data[:, 0])*sign]])
    e, s, _, _ = one_element_abaqus(simulation_dir, material=neu_sehitoglu,
                                    boundary_conditions=[BC(strain_bc, 'z', 'strain')],
                                    simulation_name=name,
                                    temperature=np.array([[0, temperature], [1, temperature]]),
                                    user_subroutine=umat_file,
                                    max_increment=0.01)
    return e, s


umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/fem_py/abaqus_material_test/neu_sehitoglu')

neu_sehitoglu.sde = 0
temp = 22
for direction, color in zip(['compression', 'tension'], ['r', 'b']):
    data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig2_' + direction),
                         delimiter=',')
    plt.figure(1)
    plt.plot(data[:, 0], data[:, 1], color + '*')
    sim_name = 'stress_22_' + direction
    sign = 1
    if direction == 'compression':
        sign = -1
    strain, stress = run_sim_from_experiment(sim_name, temp, data, sign=sign)

    plt.plot(np.abs(strain[:, 2]), np.abs(stress[:, 2]), color, lw=2)

    plt.figure(2)
    dv = np.sum(strain[:, 0:3], 1) - stress[:, 2]/neu_sehitoglu.E*(1-2*neu_sehitoglu.v)
    print dv
    plt.plot(np.abs(strain[:, 2]), dv, color, lw=2)

plt.show()
