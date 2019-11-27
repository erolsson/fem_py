import os

import numpy as np

import matplotlib.pyplot as plt

from uniaxial_experiments import experiments
from materials.transformation_materials import hazar_et_al
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC

umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/fem_py/abaqus_material_test/hazar_et_al')

for experiment in experiments:
    experiment.plot_stress_strain()
    stress_bc = BC(amplitude=[[0, 0], [1., 1.2*experiment.stress_strain[-1, 1]]], direction='z', mode='stress')
    if experiment.stress_strain[-1, 1] > 0:
        name = 'tension'
    else:
        name = 'compression'
    e, s, _, _ = one_element_abaqus(simulation_dir, material=hazar_et_al,
                                    boundary_conditions=[stress_bc],
                                    simulation_name=name + '_' + str(int(experiment.temperature)),
                                    temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                    user_subroutine=umat_file,
                                    max_increment=0.1)
    print e
    plt.plot(e, s, '--' + experiment.color)

plt.show()
