from __future__ import print_function
import os

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from scipy.optimize import fmin

from fem_py.hazar_et_al.uniaxial_experiments import experiments
from fem_py.materials.transformation_materials import hazar_et_al
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
simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/hazar_et_al/')

signed_parameters = ['Mss', 'g0']

fit_lines = []


def write_initial_file(fsb0, sim_dir, initial_austenite=0.22):
    with open(sim_dir + '/austenite.dat', 'w') as initial_file:
        for gp in (range(1, 9)):
            initial_file.write('1, ' + str(gp) + ', ' + str(initial_austenite) + ', ' + str(fsb0) + '\n')


def run_fe_simulation(parameter_values, experiment, parameter_names):

    material = hazar_et_al
    try:
        fsb0 = parameter_values[parameter_names.index('fsb0')]
    except ValueError:
        fsb0 = material.fsb0

    write_initial_file(fsb0, simulation_dir)

    for par_value, par_name in zip(parameter_values, parameter_names):
        if par_name not in signed_parameters:
            par_value = abs(par_value)
        if par_name not in ['fsb0']:
            setattr(material, par_name, par_value)

    max_exp_strain = np.max(abs(experiment.stress_strain[:, 0]))
    if experiment.volume_expansion is not None:
        e_max = max(abs(max_exp_strain), np.max(abs(experiment.transformation_data[:, 0])),
                    1.1*np.max(abs(experiment.volume_expansion[:, 0])))
    else:
        e_max = max(max_exp_strain, np.max(abs(experiment.transformation_data[:, 0])))
    if experiment.mode == 'compression':
        e_max *= -1
    strain_bc = BC(amplitude=[[0, 0], [1., 1.1*e_max]], direction='z', mode='strain')

    e, s, _, fm = one_element_abaqus(simulation_dir, material=hazar_et_al,
                                     boundary_conditions=[strain_bc],
                                     simulation_name=experiment.mode + '_' + str(int(experiment.temperature)),
                                     temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                     user_subroutine=umat_file,
                                     max_increment=0.01, output=False)
    return e, s, fm


def residual(par, *data):
    if len(fit_lines) > 0:
        for line in fit_lines:
            line.remove()
    fit_lines[:] = []
    parameter_names, experiment_list, parameter_bounds = data
    res = 0
    for bounded_par_name, bound in parameter_bounds.items():
        try:
            idx = parameter_names.index(bounded_par_name)
            if bound[0] is not None and par[idx] < bound[0]:
                par[idx] = bound[0]
            if bound[1] is not None and par[idx] > bound[1]:
                par[idx] = bound[1]
        except ValueError:
            pass

    for experiment in experiment_list:
        if experiment.mode == 'compression':
            fig = 3
        else:
            fig = 0
        e_fem, s_fem, fm_fem = run_fe_simulation(par, experiment, parameter_names)
        e_exp = experiment.stress_strain[:, 0]
        s_exp = experiment.stress_strain[:, 1]
        if e_exp[-1] > 0:
            s_interp = np.interp(e_exp, e_fem[:, 2], s_fem[:, 2])
        else:
            print(e_exp[::-1], e_fem[::-1, 2], s_fem[::-1, 2])
            s_interp = np.interp(e_exp[::-1], e_fem[::-1, 2], s_fem[::-1, 2])
            s_interp = s_interp[::-1]
        stress_residual = np.sum((1 - s_interp[abs(s_exp) > 500]/s_exp[abs(s_exp) > 500])**2) / \
            s_exp[abs(s_exp) > 500].shape[0]
        plt.figure(fig)
        fit_lines.append(plt.plot(abs(e_exp), abs(s_interp), '--x' + experiment.color, lw=2)[0])

        print('=======================================================================================================')
        print(' *** *** *** Temperature ' + str(experiment.temperature) + ' *** *** ***')
        print('=======================================================================================================')
        print("Stress at end of test: " + str(np.interp(e_exp[-1], e_fem[:, 2], s_fem[:, 2])) +
              " Exp. is " + str(s_exp[-1]) + " R = " + str(stress_residual))

        fm_exp = experiment.transformation_data[:, 1]
        fm_interp = np.interp(experiment.transformation_data[:, 0], e_fem[:, 2], fm_fem)
        plt.figure(fig+1)
        fit_lines.append(plt.plot(abs(e_fem[:, 2]), fm_fem, '--x' + experiment.color, lw=2)[0])
        martensite_residual = np.sum((fm_exp - fm_interp)**2)/fm_exp.shape[0]
        print('Martensite fractions is ' + str(fm_interp) + ' Exp. is '
              + str(fm_exp) + " R = " + str(martensite_residual))

        volume_residual = 0
        inelastic_strain = e_fem[:, 2] - s_fem[:, 2]/hazar_et_al.E
        vol_fem = np.sum(e_fem[:, 0:3], 1) - s_fem[:, 2]/hazar_et_al.E*(1 - 2*hazar_et_al.v)
        plt.figure(fig+2)
        fit_lines.append(plt.plot(abs(inelastic_strain), vol_fem, '--' + experiment.color, lw=2)[0])
        if experiment.volume_expansion is not None and experiment.mode == 'tension':
            vol_exp = experiment.volume_expansion[:, 1]
            dv_interp = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, vol_fem)
            fm_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, fm_fem)
            e_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, e_fem[:, 2])
            volume_residual = np.sum((1 - vol_exp/dv_interp)**2)
            print('Volume expansion is ' + str(dv_interp) + ' Exp. is '
                  + str(vol_exp) + ' with a martensite fraction of ' + str(fm_vol) + ' at strain' + str(e_vol)
                  + " R = " + str(volume_residual))

        res += stress_residual + martensite_residual + volume_residual
        plt.draw()
        plt.pause(0.001)
    parameter_str = ''
    for name, val in zip(parameter_names, par):
        parameter_str += name + '=' + str(val) + ', '

    print('=======================================================================================================')
    print(parameter_str)
    print('Residual: ' + str(res))
    print('')

    return res


if __name__ == '__main__':
    # parameters = {'beta': 1176, 'fsb0:': 0.0602, 'alpha': 110.435, 'R1': 0.01046,
    #               'R2': 0.008255, 'a1': 0.0258, 'Mss': -128.44}
    # parameters = {'beta': 916, 'g0': 0., 'g1': 5.244, 'M_sigma': 22, 'M_d': 350}
    bounds = {'M_sigma': (22., 75), 'M_d': (160, None)}
    parameters = {'beta': 100}
    experiments = experiments[4:5]
    plt.figure(0)
    plt.ion()
    plt.show()

    for experiment in experiments:
        if experiment.mode == 'compression':
            fig = 3
        else:
            fig = 0

        plt.figure(fig)
        experiment.plot_stress_strain()
        plt.figure(fig+1)
        experiment.plot_transformation()
        plt.figure(fig+2)
        experiment.plot_volume_expansion()
    plt.draw()
    plt.pause(0.01)

    print(fmin(residual, list(parameters.values()), args=(list(parameters.keys()), experiments, bounds),
               maxfun=1e6, maxiter=1e6))
    plt.show()
