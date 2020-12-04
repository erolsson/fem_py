from __future__ import print_function
from collections import namedtuple
import os

import numpy as np

from scipy.optimize import fmin

import matplotlib
import matplotlib.pyplot as plt

from fem_py.materials.transformation_materials import neu_sehitoglu
from fem_py.abaqus_material_test.material_test import one_element_abaqus
from fem_py.abaqus_material_test.one_element_input_file import BC
from transformation_fatigue.multiprocesser.multiprocesser import multi_processer
from fem_py.SS2506.fit_transformation_parameters import write_initial_file

matplotlib.style.use('classic')
matplotlib.use("Qt5agg")
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

umat_file = os.path.expanduser('~/python_projects/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/neu_sehitoglu/')

signed_parameters = ['Mss']
fit_lines = []


Experiment = namedtuple('Experiment', ['name', 'temperature', 'boundary_condition', 'stress_data', 'volume_data'])


def mypause(interval):
    backend = plt.rcParams['backend']
    if backend in matplotlib.rcsetup.interactive_bk:
        fig_manager = matplotlib._pylab_helpers.Gcf.get_active()
        if fig_manager is not None:
            canvas = fig_manager.canvas
            if canvas.figure.stale:
                canvas.draw()
            canvas.start_event_loop(interval)
            return


def run_sim_from_experiment(experiment):
    e, s, _, _ = one_element_abaqus(simulation_dir, material=neu_sehitoglu,
                                    boundary_conditions=experiment.boundary_condition,
                                    simulation_name=experiment.name,
                                    temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                    user_subroutine=umat_file,
                                    output=False,
                                    max_increment=0.01)
    return e, s


def residual(parameter_values, parameter_names, experiments):
    res = 0
    if len(fit_lines) > 0:
        for i in range(len(fit_lines)):
            fit_lines[i].remove()
        fit_lines[:] = []
    write_initial_file(0, simulation_dir, initial_austenite=0.35)
    for par_value, par_name in zip(parameter_values, parameter_names):
        if par_name not in signed_parameters:
            par_value = abs(par_value)
        setattr(neu_sehitoglu, par_name, par_value)

    job_list = []
    for experiment in experiments:
        job_list.append((run_sim_from_experiment, [experiment], {}))
    results = multi_processer(job_list, delay=0., timeout=300.)
    plt.figure(0)
    fit_lines.append(plt.plot(results[0][0][:, 2], results[0][1][:, 2], '--b', lw=2)[0])
    fit_lines.append(plt.plot(-results[1][0][:, 2], -results[1][1][:, 2], '--r', lw=2)[0])

    plt.figure(1)
    dv_tension = np.sum(results[0][0][:, 0:3], 1) - results[0][1][:, 2]/neu_sehitoglu.E*(1 - 2*neu_sehitoglu.v)
    fit_lines.append(plt.plot(results[0][0][:, 2], dv_tension, '--b', lw=2)[0])

    dv_comp = np.sum(results[1][0][:, 0:3], 1) - results[1][1][:, 2]/neu_sehitoglu.E*(1 - 2*neu_sehitoglu.v)
    fit_lines.append(plt.plot(-results[1][0][:, 2], dv_comp, '--r', lw=2)[0])

    plt.figure(2)
    gamma = results[2][0][:, 1] - results[2][0][:, 2]
    tau = (results[2][1][:, 1] - results[2][1][:, 2])/2
    fit_lines.append(plt.plot(gamma, tau, '--k')[0])

    plt.figure(3)
    dv = np.sum(results[2][0][:, 0:3], axis=1)
    fit_lines.append(plt.plot(tau, dv, '--k')[0])

    plt.draw()
    mypause(0.001)
    num_points = 100
    min_strain = min(results[0][0][-1, 2], experiments[0].stress_data[-1, 0])
    strain_points = np.linspace(1e-4, min_strain, num_points)
    exp_stress = np.interp(strain_points, experiments[0].stress_data[:, 0], experiments[0].stress_data[:, 1])
    num_stress = np.interp(strain_points, results[0][0][:, 2], results[0][1][:, 2])
    res += np.sum((num_stress - exp_stress)**2)/np.max(exp_stress)**2

    min_strain = min(results[0][0][-1, 2], experiments[0].volume_data[-1, 0])
    strain_points = np.linspace(1e-4, min_strain, num_points)
    exp_volume = np.interp(strain_points, experiments[0].volume_data[:, 0], experiments[0].volume_data[:, 1])
    num_dv = np.interp(strain_points, results[0][0][:, 2], dv_tension)
    res += np.sum((num_dv - exp_volume)**2)/np.max(exp_volume)**2

    min_strain = min(-results[1][0][-1, 2], experiments[1].stress_data[-1, 0])
    strain_points = np.linspace(1e-4, min_strain, num_points)
    exp_stress = np.interp(strain_points, experiments[1].stress_data[:, 0], experiments[1].stress_data[:, 1])
    num_stress = np.interp(strain_points, -results[1][0][:, 2], -results[1][1][:, 2])
    # res += np.sum((num_stress - exp_stress)**2)/np.max(exp_stress)**2

    min_strain = min(-results[1][0][-1, 2], experiments[1].volume_data[-1, 0])
    strain_points = np.linspace(1e-4, min_strain, num_points)
    exp_volume = np.interp(strain_points, experiments[0].volume_data[:, 0], experiments[0].volume_data[:, 1])
    num_dv = np.interp(strain_points, -results[1][0][:, 2], dv_comp)
    # res += np.sum((num_dv - exp_volume)**2)/np.max(exp_volume)**2

    strain_points = np.linspace(1e-4, 0.04, num_points)
    exp_stress = np.interp(strain_points, experiments[2].stress_data[:, 0], experiments[2].stress_data[:, 1])
    num_stress = np.interp(strain_points, gamma, tau)
    # res += np.sum((num_stress - exp_stress)**2)/np.max(exp_stress)**2

    stress_points = np.linspace(0, 2000, num_points)
    exp_volume = np.interp(stress_points, experiments[2].volume_data[:, 0], experiments[2].volume_data[:, 1])
    num_volume = np.interp(stress_points, tau, dv)
    res += np.sum((num_volume - exp_volume)**2)/np.max(exp_volume)**2

    parameter_str = ''
    for name, val in zip(parameter_names, parameter_values):
        parameter_str += name + '=' + str(val) + ', '

    print('=======================================================================================================')
    print('\t' + parameter_str)
    print('\tResidual: ' + str(res))
    print('')

    return res


def main():
    tension_stress_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig2_tension'),
                                        delimiter=',')
    idx = np.argsort(tension_stress_data[:, 0])
    tension_stress_data = tension_stress_data[idx, :]
    plt.figure(0)
    plt.plot(tension_stress_data[:, 0], tension_stress_data[:, 1], '-b*', ms=12)
    tension_volume_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig3_tension'),
                                        delimiter=',')
    idx = np.argsort(tension_volume_data[:, 0])
    tension_volume_data = tension_volume_data[idx, :]
    plt.figure(1)
    plt.plot(tension_volume_data[:, 0], tension_volume_data[:, 1], '-b*', ms=12)

    comp_stress_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig2_compression'),
                                     delimiter=',')
    idx = np.argsort(comp_stress_data[:, 0])
    comp_stress_data = comp_stress_data[idx, :]
    plt.figure(0)
    plt.plot(comp_stress_data[:, 0], comp_stress_data[:, 1], '-r*', ms=12)
    comp_volume_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig3_compression'),
                                     delimiter=',')
    idx = np.argsort(comp_volume_data[:, 0])
    comp_volume_data = comp_volume_data[idx, :]
    plt.figure(1)
    plt.plot(comp_volume_data[:, 0], comp_volume_data[:, 1], '-r*', ms=12)

    plt.figure(2)
    torsion_stress_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig7'), delimiter=',')
    torsion_stress_data[:, 0] *= 2
    plt.plot(torsion_stress_data[:, 0], torsion_stress_data[:, 1], '-b*', ms=12)
    plt.xlim(0, 0.14)
    plt.ylim(0, 2000)

    plt.figure(3)
    torsion_volume_data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig8'), delimiter=',')
    plt.plot(torsion_volume_data[:, 0], torsion_volume_data[:, 1], '-b*', ms=12)
    plt.xlim(0, 2000)
    experiments = [Experiment(name='tension', temperature=22.,
                              boundary_condition=[BC(amplitude=[(0., 0.), (1., 1300)],
                                                     direction='z',
                                                     mode='stress')],
                              stress_data=tension_stress_data,
                              volume_data=tension_volume_data),
                   Experiment(name='compression', temperature=22.,
                              boundary_condition=[BC(amplitude=[(0., 0.), (1., -2400)],
                                                     direction='z',
                                                     mode='stress')],
                              stress_data=comp_stress_data,
                              volume_data=comp_volume_data),
                   Experiment(name='torsion', temperature=22.,
                              boundary_condition=[BC(amplitude=[(0., 0.), (1., 4000)],
                                                     direction='y', mode='stress'),
                                                  BC(amplitude=[(0., 0.), (1., -4000)],
                                                     direction='z', mode='stress')],
                              stress_data=torsion_stress_data,
                              volume_data=torsion_volume_data)
                   ]
    plt.draw()
    plt.pause(0.001)
    plt.ion()
    plt.show()
    neu_sehitoglu.Ms = 185.51653373817538
    parameters = {'a1': 0.028107490537271854,
                  'a2': 1.2959257140209983e-05,
                  'a3': 2.8149370838159175e-07,
                  'R1': -6.014172555589034e-07,
                  'R2': 0.024553598689201442,
                  'Mss': -45}
    print(fmin(residual, list(parameters.values()), args=(list(parameters.keys()), experiments)))
    plt.show()


if __name__ == '__main__':
    main()
