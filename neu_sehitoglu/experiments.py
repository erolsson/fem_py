from collections import namedtuple
import os
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


matplotlib.style.use('classic')
matplotlib.use("Qt5agg")
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


Serie = namedtuple('Serie', ['tension', 'compression', 'torsion'])


class Experiment:
    def __init__(self, stress_data_file, volume_data_file, num_points=100, mode=''):
        stress_data = np.genfromtxt(stress_data_file, delimiter=', ')
        volume_data = np.genfromtxt(volume_data_file, delimiter=', ')

        if mode == 'torsion':
            min_stress = np.min([stress_data[-1, 1], volume_data[-1, 0]])
            min_strain = np.interp(min_stress, stress_data[:, 1], stress_data[:, 0])
            self.strain = np.linspace(0, min_strain, num_points)
            self.stress = np.interp(self.strain, stress_data[:, 0], stress_data[:, 1])
            self.dv = np.interp(self.stress, volume_data[:, 0], volume_data[:, 1])
            # self.strain *= 2
        else:
            min_strain = np.min([stress_data[-1, 0], volume_data[-1, 0]])
            self.strain = np.linspace(0, min_strain, num_points)
            self.dv = np.interp(self.strain, volume_data[:, 0], volume_data[:, 1])
            self.stress = np.interp(self.strain, stress_data[:, 0], stress_data[:, 1])


data_directory = os.path.expanduser('~/phase_transformations/neu_sehitoglu/')

experiments = Serie(tension=Experiment(stress_data_file=os.path.join(data_directory, 'fig2_tension'),
                                       volume_data_file=os.path.join(data_directory, 'fig3_tension')),
                    compression=Experiment(stress_data_file=os.path.join(data_directory, 'fig2_compression'),
                                           volume_data_file=os.path.join(data_directory, 'fig3_compression')),
                    torsion=Experiment(stress_data_file=os.path.join(data_directory, 'fig7'),
                                       volume_data_file=os.path.join(data_directory, 'fig8'), mode='torsion'))


def main():
    plt.figure(0)
    plt.plot(experiments.tension.strain, experiments.tension.stress, '-xb', lw=2, ms=12)
    plt.plot(experiments.compression.strain, experiments.compression.stress, '-xr', lw=2, ms=12)

    plt.figure(1)
    plt.plot(experiments.torsion.strain, experiments.torsion.stress, '-xk', lw=2, ms=12)

    plt.figure(2)
    plt.plot(experiments.tension.stress, experiments.tension.dv, '-xb', lw=2, ms=12)
    plt.plot(experiments.compression.stress, experiments.compression.dv, '-xr', lw=2, ms=12)
    plt.plot(experiments.torsion.stress, experiments.torsion.dv, '-xk', lw=2, ms=12)
    plt.show()


if __name__ == '__main__':
    main()
