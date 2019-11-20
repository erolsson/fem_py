import os

import numpy as np

import matplotlib.pyplot as plt

from materials.SS2506 import neu_sehitoglu
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC

umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')

for temp, color in zip([150., 22.], ['r', 'b']):
    data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig_4_' + str(int(temp)) + 'C'),
                         delimiter=',')
    plt.plot(data[:, 0], data[:, 1], color + '*')
    stress_bc = np.array([[0., 0], [1., data[-1, 1]]])
    strain_bc = np.array([[0., 0], [1., data[0, 1]]])

    stress, strain, _, _ = one_element_abaqus('abaqus_material_test/neu_sehitoglu/', material=neu_sehitoglu,
                                              boundary_conditions=[BC(stress_bc, 'z', 'stress')],
                                              simulation_name='stress_temp',
                                              temperature=np.array([[0, temp], [1, temp]]), user_subroutine=umat_file)
    plt.plot(strain[:, 2], stress[:, 2], color, lw=2)
plt.show()
