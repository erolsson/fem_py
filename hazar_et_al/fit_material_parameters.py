from __future__ import print_function
import os

import numpy as np

from uniaxial_experiments import experiments
from materials.transformation_materials import hazar_et_al
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC

umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/fem_py/abaqus_material_test/hazar_et_al')


def run_fe_simulation(parameter_values, experiment, parameter_names):
    material = hazar_et_al
    for par_value, par_name in zip(parameter_values, parameter_names):
        setattr(material, par_name, par_value)
    strain_bc = BC(amplitude=[[0, 0], [1., 1.*experiment.stress_strain[-1, 0]]], direction='z', mode='strain')
    if experiment.stress_strain[-1, 1] > 0:
        name = 'tension'
    else:
        name = 'compression'
    e, s, _, fm = one_element_abaqus(simulation_dir, material=hazar_et_al,
                                     boundary_conditions=[strain_bc],
                                     simulation_name=name + '_' + str(int(experiment.temperature)),
                                     temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                     user_subroutine=umat_file,
                                     max_increment=0.01)
    return e, s, fm


def residual(par, data):
    parameter_names, experiment_list = data


if __name__ == '__main__':
    run_fe_simulation([0], experiments[1], ['R2'])
