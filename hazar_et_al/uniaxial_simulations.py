from __future__ import print_function
import os

import numpy as np

import matplotlib.pyplot as plt

from uniaxial_experiments import experiments
from fit_material_parameters import write_initial_file

from fem_py.materials.transformation_materials import hazar_et_al
from fem_py.abaqus_material_test.material_test import one_element_abaqus
from fem_py.abaqus_material_test.one_element_input_file import BC

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

umat_file = os.path.expanduser('~/python_projects/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/hazar_et_al')

write_initial_file(hazar_et_al.fsb0, simulation_dir)

for experiment in experiments:
    fig_idx = 0
    if experiment.stress_strain[-1, 1] > 0:
        name = 'tension'
    else:
        name = 'compression'
        fig_idx = 2
    plt.figure(0 + fig_idx)
    experiment.plot_stress_strain()
    plt.figure(1 + fig_idx)
    experiment.plot_transformation()
    stress_bc = BC(amplitude=[[0, 0], [1., 1.2*experiment.stress_strain[-1, 0]]], direction='z', mode='strain')

    e, s, _, fm = one_element_abaqus(simulation_dir, material=hazar_et_al,
                                     boundary_conditions=[stress_bc],
                                     simulation_name=name + '_' + str(int(experiment.temperature)),
                                     temperature=np.array([[0, experiment.temperature], [1, experiment.temperature]]),
                                     user_subroutine=umat_file,
                                     max_increment=0.01)
    
    plt.figure(0 + fig_idx)
    plt.plot(abs(e[:, 2]), abs(s[:, 2]), '--' + experiment.color, lw=3)
    plt.figure(1 + fig_idx)
    plt.plot(abs(e[:, 2]), fm, '--' + experiment.color, lw=3)

for i, fig in enumerate('abcd'):
    plt.figure(i)
    plt.xlabel('Strain [-]')
    ax = plt.gca()
    plt.text(0.05, 0.9, '(' + fig + ')', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    if i % 2 == 0:
        plt.ylabel('Stress [MPa]')
    else:
        plt.ylabel('Fraction martensite [-]')
    if i == 0:
        plt.legend(loc='upper left', bbox_to_anchor=(0.5, 0.5), framealpha=0.9)
    plt.tight_layout()
plt.show()
