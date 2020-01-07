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

signed_parameters = ['Mss']

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
        if par_name != 'fsb0':
            setattr(material, par_name, par_value)
    stabilization_temp = -np.log(1-0.78)/material.k - material.Ms - material.a1*800 + 22
    material.Mss = stabilization_temp
    max_exp_strain = np.max(experiment.stress_strain[:, 0])
    if experiment.volume_expansion is not None:
        e_max = max(max_exp_strain, np.max(experiment.transformation_data[:, 0]),
                    1.1*np.max(experiment.volume_expansion[:, 0]))
    else:
        e_max = max(max_exp_strain, np.max(experiment.transformation_data[:, 0]))
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
    if len(fit_lines) > 0:
        for line in fit_lines:
            line.remove()
    fit_lines[:] = []
    parameter_names, experiment_list = data
    res = 0
    for experiment in experiment_list:
        if experiment.mode == 'compression':
            fig = 3
        else:
            fig = 0
        e_fem, s_fem, fm_fem = run_fe_simulation(par, experiment, parameter_names)
        e_exp = experiment.stress_strain[:, 0]
        s_exp = experiment.stress_strain[:, 1]
        s_intep = np.interp(e_exp, e_fem[:, 2], s_fem[:, 2])
        stress_residual = np.sum((1 - s_intep/s_exp)**2)/s_exp.shape[0]
        plt.figure(fig)
        fit_lines.append(plt.plot(e_exp, s_intep, '--x' + experiment.color, lw=2)[0])

        print('=======================================================================================================')
        print(' *** *** *** Temperature ' + str(experiment.temperature) + ' *** *** ***')
        print('=======================================================================================================')
        print("Stress at end of test: " + str(np.interp(e_exp[-1], e_fem[:, 2], s_fem[:, 2])) +
              " Exp. is " + str(s_exp[-1]))
        martensite_residual = 0

        fm_exp = experiment.transformation_data[:, 1]
        fm_interp = np.interp(experiment.transformation_data[:, 0], e_fem[:, 2], fm_fem)
        plt.figure(fig+1)
        fit_lines.append(plt.plot(e_exp, fm_interp, '--x' + experiment.color, lw=2)[0])
        if experiment.temperature < 120:
            print('Martensite fractions is ' + str(fm_interp) + ' Exp. is '
                  + str(fm_exp))
            martensite_residual = np.sum((fm_exp - fm_interp)**2)/fm_exp.shape[0]/(np.max(fm_fem) - 0.78)**2
        volume_residual = 0
        inelastic_strain = e_fem[:, 2] - s_fem[:, 2]/hazar_et_al.E
        vol_fem = np.sum(e_fem[:, 0:3], 1) - s_fem[:, 2]/hazar_et_al.E*(1 - 2*hazar_et_al.v)
        plt.figure(fig+2)
        fit_lines.append(plt.plot(inelastic_strain, vol_fem, '--' + experiment.color, lw=2)[0])
        if experiment.volume_expansion is not None:
            vol_exp = experiment.volume_expansion[:, 1]
            dv_interp = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, vol_fem)
            fm_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, fm_fem)
            e_vol = np.interp(experiment.volume_expansion[:, 0], inelastic_strain, e_fem[:, 2])
            volume_residual = np.sum((1 - vol_exp/dv_interp)**2)
            print('Volume expansion is ' + str(dv_interp) + ' Exp. is '
                  + str(vol_exp) + ' with a martensite fraction of ' + str(fm_vol) + ' at strain' + str(e_vol))
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
    parameters = {'fsb0': 0.11879158373330523,
                  'R1': 0.01599183565870769, 'R2': 0.00793495890945268,
                  'g1': 60.62406811305142, 'g2': 5.,
                  'g_std': 18.056870261724143, 'g0': 14.361769920692034}
    experiments = experiments[1:]
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

    plt.ioff()
    print(fmin(residual, list(parameters.values()), args=(list(parameters.keys()), experiments[1:]),
               maxfun=1e6, maxiter=1e6))
    plt.show()
