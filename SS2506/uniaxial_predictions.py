from __future__ import print_function
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from fit_transformation_parameters import run_fe_simulation

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

for name, color in zip(['Tension', 'Compression'], ['b', 'r']):
    data = np.genfromtxt('experimental_data/' + name.lower() + '_data_case_EN.dat', delimiter=',')
    plt.figure(0)
    plt.plot(abs(data[:, 0]), abs(data[:, 1]), color, lw=3, label=name)
    plt.figure(1)
    strains = np.genfromtxt('experimental_data/transversal_strain_' + name.lower() + '_EN.dat', delimiter=',')
    plt.plot(abs(strains[:, 0]), abs(strains[:, 1]), color, lw=3, label=name)

    e, s, _ = run_fe_simulation([], data, [])
    print(s, e)
    plt.figure(0)
    plt.plot(abs(e[:, 2]), abs(s[:, 2]), color + '--', lw=3)

    plt.figure(1)
    plt.plot(abs(e[:, 2]), abs(e[:, 0]), color + '--', lw=3)

plt.figure(0)
plt.xlim(0, 0.01)
plt.ylim(0, 2000)
plt.legend(loc='best')
plt.xlabel('Strain [-]')
plt.ylabel('Stress [MPa]')
plt.tight_layout()
plt.savefig('stress_strain2506.png')

plt.figure(1)
plt.xlim(0, 0.01)
plt.ylim(0, 0.0035)
plt.legend(loc='best')
plt.xlabel('Longitudinal strain magnitude  [-]')
plt.ylabel('Transversal strain magnitude [-]')
plt.tight_layout()
plt.savefig('strains.png')

plt.show()
