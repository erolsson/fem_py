from __future__ import print_function
import os

import numpy as np

from scipy.optimize import fmin

from uniaxial_experiments import experiments
from materials.transformation_materials import hazar_et_al
from abaqus_material_test.material_test import one_element_abaqus
from abaqus_material_test.one_element_input_file import BC

umat_file = os.path.expanduser('~/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/fem_py/abaqus_material_test/hazar_et_al/')

signed_parameters = ['Mss']


def write_initial_file(fsb0, initial_austenite=0.22):
    with open(simulation_dir + '/austenite.dat', 'w') as initial_file:
        for gp in (range(1, 9)):
            initial_file.write('1, ' + str(gp) + ', ' + str(initial_austenite) + ', ' + str(fsb0) + '\n')


def run_fe_simulation(parameter_values, experiment, parameter_names):
    material = hazar_et_al
    try:
        fsb0 = parameter_values[parameter_names.index('fsb0')]
    except ValueError:
        fsb0 = material.fsb0
    write_initial_file(fsb0)
    for par_value, par_name in zip(parameter_values, parameter_names):
        if par_name not in signed_parameters:
            par_value = abs(par_value)
        if par_name != 'fsb0':
            setattr(material, par_name, par_value)
    if experiment.volume_expansion is not None:
        e_max = max(experiment.stress_strain[-1, 0], np.max(experiment.transformation_data[:, 0]),
                    1.1*np.max(experiment.volume_expansion[:, 0]))
    else:
        e_max = max(experiment.stress_strain[-1, 0], np.max(experiment.transformation_data[:, 0]))
    strain_bc = BC(amplitude=[[0, 0], [1., 1.1*e_max]], direction='z', mode='strain')
    if experiment.stress_strain[-1, 1] > 0:
        name = 'tension'
    else:
        name = 'compression'

    e, s, _, fm = one_element_abaqus(simulation_dir, material=hazar_et_al,
                                     boundary_conditions=[strain_bc],
                                     simulation_name=name + '_' + str(int(experiment.temperature)),
                                     temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                     user_subroutine=umat_file,
                                     max_increment=0.01, output=False)
    return e, s, fm


def residual(par, *data):
    parameter_names, experiment_list = data
    res = 0
    for experiment in experiment_list:
        e_fem, s_fem, fm_fem = run_fe_simulation(par, experiment, parameter_names)
        e_exp = experiment.stress_strain[:, 0]
        s_exp = experiment.stress_strain[:, 1]
        s_intep = np.interp(e_exp, e_fem[:, 2], s_fem[:, 2])
        stress_residual = np.sum((1 - s_intep[s_exp > 0]/s_exp[s_exp > 0])**2)/s_exp.shape[0]
        print('=======================================================================================================')
        print(' *** *** *** Temperature ' + str(experiment.temperature) + ' *** *** ***')
        print('=======================================================================================================')
        print("Stress at end of test: " + str(np.interp(e_exp[-1], e_fem[:, 2], s_fem[:, 2])) +
              " Exp. is " + str(s_exp[-1]))
        fm_exp = experiment.transformation_data[:, 1]
        fm_interp = np.interp(experiment.transformation_data[:, 0], e_fem[:, 2], fm_fem)
        print('Martensite fractions is ' + str(fm_interp) + ' Exp. is '
              + str(fm_exp))
        martensite_residual = np.sum((fm_exp - fm_interp)**2)/fm_exp.shape[0]/(np.max(fm_fem) - 0.78)**2
        volume_residual = 0
        if experiment.volume_expansion is not None:
            inelastic_strain = e_fem[:, 2] - s_fem[:, 2]/hazar_et_al.E
            vol_fem = np.sum(e_fem[:, 0:3], 1) - s_fem[:, 2]/hazar_et_al.E*(1 - 2*hazar_et_al.v)
            vol_exp = experiment.volume_expansion[:, 1]
            dv_interp = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, vol_fem)
            fm_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, fm_fem)
            e_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, e_fem[:, 2])
            volume_residual = np.sum((dv_interp - vol_exp)**2)/np.max(vol_exp)**2
            print('Volume expansion is ' + str(dv_interp) + ' Exp. is '
                  + str(vol_exp) + ' with a martensite fraction of ' + str(fm_vol) + ' at strain' + str(e_vol))
        res += stress_residual + martensite_residual + volume_residual

    parameter_str = ''
    for name, val in zip(parameter_names, par):
        parameter_str += name + '=' + str(val) + ', '

    print('=======================================================================================================')
    print(parameter_str)
    print('Residual:' + str(res))
    print('')

    return res


if __name__ == '__main__':
    parameters = {'beta': 4.92067433e+02, 'alpha': 1.06625032e+02, 'a1': 0.05, 'Mss': -86., 'fsb0': 6.40671990e-02,
                  'R1': 2.46866223e-02, 'R2': 2.83373527e-03, 'dV': 2.35190122e-02, 'g1': 1.78036231e+02,
                  'g_std': 1.78036231e+02}
    print(fmin(residual, list(parameters.values()), args=(list(parameters.keys()), experiments),
               maxfun=1e6, maxiter=1e6))
