from __future__ import print_function
import os

import numpy as np

from scipy.optimize import fmin

import matplotlib
import matplotlib.pyplot as plt

from fem_py.materials.transformation_materials import SS2506
from fem_py.abaqus_material_test.material_test import one_element_abaqus
from fem_py.abaqus_material_test.one_element_input_file import BC

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

umat_file = os.path.expanduser('~/python_projects/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/SS2506/')

signed_parameters = ['Mss']
fit_lines = []


def write_initial_file(fsb0, sim_dir, initial_austenite=0.20):
    with open(sim_dir + '/austenite.dat', 'w') as initial_file:
        for gp in (range(1, 9)):
            initial_file.write('1, ' + str(gp) + ', ' + str(initial_austenite) + ', ' + str(fsb0) + '\n')


def run_fe_simulation(parameter_values, experiment, parameter_names):
    material = SS2506
    try:
        fsb0 = parameter_values[parameter_names.index('fsb0')]
    except ValueError:
        fsb0 = material.fsb0
    write_initial_file(fsb0, simulation_dir)
    for par_value, par_name in zip(parameter_values, parameter_names):
        if par_name not in signed_parameters:
            par_value = abs(par_value)
        if par_name != 'fsb0':
            setattr(material, par_name, par_value)

    max_exp_strain = np.max(experiment[:, 0])
    strain_bc = BC(amplitude=[[0, 0], [1., 1.1*max_exp_strain]], direction='z', mode='strain')
    if experiment[-1, 1] > 0:
        name = 'tension'
    else:
        name = 'compression'

    e, s, _, fm = one_element_abaqus(simulation_dir, material=material,
                                     boundary_conditions=[strain_bc],
                                     simulation_name=name,
                                     temperature=np.array([[0, 22.], [1, 22.]]),
                                     user_subroutine=umat_file,
                                     max_increment=0.01, output=False)
    return e, s, fm


def residual(par, *data):
    parameter_names, experiment = data
    res = 0
    if len(fit_lines) > 0:
        fit_lines[0].remove()
        fit_lines[1].remove()
        fit_lines[:] = []
    e_fem, s_fem, fm_fem = run_fe_simulation(par, experiment, parameter_names)
    plt.figure(0)
    fit_lines.append(plt.plot(e_fem[:, 2], s_fem[:, 2], '--b', lw=2)[0])
    plt.figure(0)
    fit_lines.append(plt.plot(e_fem[:, 2], s_fem[:, 0], '--b', lw=2)[0])
    plt.draw()
    plt.pause(0.001)
    e_exp = experiment[:, 0]
    s_exp = experiment[:, 1]
    et_exp = experiment[:, 2]
    s_intep = np.interp(e_exp, e_fem[:, 2], s_fem[:, 2])
    et_intep = np.interp(e_exp, e_fem[:, 2], e_fem[:, 1])
    stress_residual = np.sum((1 - s_intep[s_exp > 300]/s_exp[s_exp > 300])**2)/s_exp.shape[0]
    print("Stress at end of test: " + str(np.interp(e_exp[-1], e_fem[:, 2], s_fem[:, 2])) +
          " Exp. is " + str(s_exp[-1]))
    print("Stress at end of test: " + str(np.interp(e_exp[-1], e_fem[:, 2], e_fem[:, 1])) +
          " Exp. is " + str(et_exp[-1]))
    strain_residual = np.sum((1 - et_intep[s_exp > 300]/et_exp[s_exp > 300])**2)/s_exp.shape[0]
    res += stress_residual + strain_residual

    parameter_str = ''
    for name, val in zip(parameter_names, par):
        parameter_str += name + '=' + str(val) + ', '

    print('=======================================================================================================')
    print(parameter_str)
    print('Residual: ' + str(res))
    print('')

    return res


if __name__ == '__main__':
    stress_strain = np.genfromtxt('experimental_data/tension_data_case_EN.dat', delimiter=',')
    strains = np.genfromtxt('experimental_data/transversal_strain_tension_EN.dat', delimiter=',')
    experimental_data = np.zeros((stress_strain.shape[0], 3))
    experimental_data[:, 0:2] = stress_strain
    experimental_data[:, 2] = strains[:, 1]
    plt.figure(0)
    plt.ion()
    plt.show()
    plt.plot(experimental_data[:, 0], experimental_data[:, 1], 'b', lw=3)
    plt.xlabel(r'$\varepsilon_{zz}$ [-]', fontsize=24)
    plt.ylabel(r'$\sigma_{zz}$ [MPa]', fontsize=24)
    plt.tight_layout()
    plt.figure(1)
    plt.plot(experimental_data[:, 0], experimental_data[:, 2], 'b', lw=3)
    plt.xlabel(r'$\varepsilon_{zz}$ [-]', fontsize=24)
    plt.ylabel(r'$\varepsilon_{xx}$ [-]', fontsize=24)
    plt.tight_layout()

    plt.draw()
    plt.pause(0.001)
    parameters = {'a1': 0.056, 'dV': 0.037,
                  'R1': 0.02, 'R2': 0.02, 'Mss': -10}
    print(fmin(residual, list(parameters.values()), args=(list(parameters.keys()), experimental_data),
               maxfun=1e6, maxiter=1e6))
