import numpy as np

import matplotlib
import matplotlib.pyplot as plt

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

case_data = -np.genfromtxt('experimental_data/compression_data_case_EN.dat', delimiter=',')
plt.plot(case_data[:, 0], case_data[:, 1], lw=3, label='Compression')
case_data = np.genfromtxt('experimental_data/tension_data_case_EN.dat', delimiter=',')
plt.plot(case_data[:, 0], case_data[:, 1], lw=3, label='Tension')

plt.xlim(0, 0.01)
plt.ylim(0, 2000)
plt.legend(loc='best')
plt.xlabel('Strain')
plt.ylabel('Stress [MPa]')
plt.tight_layout()
plt.savefig('stress_strain2506.png')
plt.show()
