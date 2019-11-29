from __future__ import print_function
import os

import numpy as np

import matplotlib.pyplot as plt

data_directory = os.path.expanduser('~/phase_transformations/hazar_et_al/')


class Experiment:
    def __init__(self, temperature, color, mode='tension'):
        self.temperature = temperature
        self.color = color
        stress_strain = np.genfromtxt(data_directory + '/fig6a_' + str(self.temperature) + 'C', delimiter=',')
        self.stress_strain = stress_strain[np.argsort(stress_strain[:, 0]), :]

        transformation_data = np.genfromtxt(data_directory + '/fig8a_' + str(self.temperature) + 'C',
                                            delimiter=',')
        if len(transformation_data.shape) < 2:
            transformation_data = np.expand_dims(transformation_data, 0)
        self.transformation_data = transformation_data

        try:
            self.volume_expansion = np.genfromtxt(data_directory + '/fig18_' + mode + '_' + str(self.temperature) + 'C',
                                                  delimiter=',')
            if len(self.volume_expansion.shape) < 2:
                self.volume_expansion = np.expand_dims(self.volume_expansion, 0)
        except IOError:
            self.volume_expansion = None

    def plot_stress_strain(self):
        plt.plot(self.stress_strain[:, 0], self.stress_strain[:, 1], self.color, lw=2)

    def plot_transformation(self):
        plt.plot(self.transformation_data[:, 0], self.transformation_data[:, 1], 'x' + self.color, ms=12, mew=2)

    def plot_volume_expansion(self):
        if self.volume_expansion is not None:
            plt.plot(self.volume_expansion[:, 0], self.volume_expansion[:, 1], 'x' + self.color, ms=12, mew=2)


experiments = [Experiment(temperature=22, color='r'), Experiment(temperature=75, color='b'),
               Experiment(temperature=100, color='g'), Experiment(temperature=150, color='k')]

if __name__ == '__main__':
    for experiment in experiments:
        plt.figure(0)
        experiment.plot_stress_strain()
        plt.figure(1)
        experiment.plot_transformation()

    par = None
    k = 0
    a = 0
    Ms = 220
    Mss = 0

    for experiment in experiments[0:2]:
        temp = experiment.temperature
        print(temp)
        if temp in [22., 75.]:
            trans_stress = np.zeros(experiment.transformation_data.shape[0] + 1)
            y = 0*trans_stress
            trans_stress[1:] = np.interp(experiment.transformation_data[:, 0], experiment.stress_strain[:, 0],
                                         experiment.stress_strain[:, 1])
            y[1:] = -np.log(1 - experiment.transformation_data[:, 1])
            y[0] = -np.log(1 - 0.78)
            if temp == 22.:
                trans_stress[0] = 850
            else:
                trans_stress[0] = 1160

            plt.figure(2)
            print(trans_stress)

            plt.plot(trans_stress, y, 'x' + experiment.color, ms=12, mew=2)
            if temp == 22.:
                par = np.polyfit(trans_stress, y, 1)

    k = 0.01
    a = par[0]/k
    Mss = -(np.log(1-0.78) + 0*3.5*0.78**4)/k - Ms - a*850 + 22
    print(k, a, Mss)
    B = k*(Ms + a*np.array([1150, 1350]) + Mss - np.array([75, 100]))
    print(B)
    print (1 - np.exp(-B))
    print(((-np.log(1-0.78) - B)/3.05)**0.25)
    for experiment in experiments:
        y = k*(Ms + a*experiment.stress_strain[:, 1] + Mss - experiment.temperature)
        plt.figure(2)
        plt.plot(experiment.stress_strain[:, 1], y, experiment.color)
        plt.figure(1)
        fM = 1 - np.exp(-y)
        plt.plot(experiment.stress_strain[:, 0], fM, '--' + experiment.color)

        plt.figure(3)
        plt.plot(experiment.stress_strain[:, 1], experiment.stress_strain[:, 0] - experiment.stress_strain[:, 1]/200.5e3,
                 experiment.color)

    plt.figure(1)
    plt.xlim(0, 0.015)
    plt.ylim(0.78, 1.)
    plt.show()
