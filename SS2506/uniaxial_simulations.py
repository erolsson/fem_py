import os

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from fit_transformation_parameters import run_fe_simulation, write_initial_file
from materials.transformation_materials import SS2506

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

umat_file = os.path.expanduser('~/python_projects/fem_py/transformation_subroutines/transformation_subroutine.o')
simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/SS2506/')


def main():
    stress_strain = np.genfromtxt('experimental_data/tension_data_case_EN.dat', delimiter=',')
    strains = np.genfromtxt('experimental_data/transversal_strain_tension_EN.dat', delimiter=',')
    experimental_data = np.zeros((stress_strain.shape[0], 3))
    experimental_data[:, 0:2] = stress_strain
    experimental_data[:, 2] = strains[:, 1]
    austenite = np.arange(0.15, 0.21, 0.01)
    print(austenite)
    SS2506.a1 = 0.08
    SS2506.a2 = 3e-4
    SS2506.a3 = 0
    SS2506.Mss = -148.38998859637542
    for ra in austenite:
        print(ra)
        e_fem, s_fem, fm_fem = run_fe_simulation([], experimental_data, [], SS2506, initial_austenite=ra)
        plt.figure(0)
        plt.plot(e_fem[:, 2], s_fem[:, 2])

    plt.figure(0)
    plt.plot(experimental_data[:, 0], experimental_data[:, 1], '-bx', lw=3)
    plt.xlabel(r'$\varepsilon_{zz}$ [-]', fontsize=24)
    plt.ylabel(r'$\sigma_{zz}$ [MPa]', fontsize=24)
    plt.tight_layout()
    plt.figure(1)
    plt.plot(experimental_data[:, 0], experimental_data[:, 2], '-bx', lw=3)
    plt.xlabel(r'$\varepsilon_{zz}$ [-]', fontsize=24)
    plt.ylabel(r'$\varepsilon_{xx}$ [-]', fontsize=24)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
