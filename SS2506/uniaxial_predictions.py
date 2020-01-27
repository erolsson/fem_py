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

for name in ['tension', 'compression']:
    data = np.genfromtxt('experimental_data/' + name + '_data_case_EN.dat', delimiter=',')
    plt.figure(0)
    plt.plot(abs(data[:, 0]), abs(data[:, 1]), lw=3, label=name)
    plt.figure(1)
    strains = np.genfromtxt('experimental_data/transversal_strain_' + name + '_EN.dat', delimiter=',')
    plt.plot(abs(strains[:, 0]), abs(strains[:, 1]), lw=3, label=name)

    e, s, _ = run_fe_simulation([], data, [])
    plt.figure(0)
    plt.plot(abs(e[:, 2]), abs(s[:, 2]), '--', lw=3)

plt.figure(0)
plt.xlim(0, 0.01)
plt.ylim(0, 2000)
plt.legend(loc='best')
plt.xlabel('Strain')
plt.ylabel('Stress [MPa]')
plt.tight_layout()
plt.savefig('stress_strain2506.png')

plt.show()
