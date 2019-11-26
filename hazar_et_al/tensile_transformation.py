from __future__ import print_function
from collections import namedtuple
import os

import numpy as np

import matplotlib.pyplot as plt

data_directory = os.path.expanduser('~/phase_transformations/hazar_et_al/')


class Experiment:
    def __init__(self, temperature, color):
        self.temperature = temperature
        self.color = color
        self.stress_strain = None
        stress_strain = np.genfromtxt(data_directory + '/fig6a_' + str(self.temperature) + 'C', delimiter=',')
        self.stress_strain = stress_strain[np.argsort(stress_strain[:, 0]), :]

        transformation_data = np.genfromtxt(data_directory + '/fig8a_' + str(self.temperature) + 'C',
                                            delimiter=',')
        if len(transformation_data.shape) < 2:
            transformation_data = np.expand_dims(transformation_data, 0)
        self.transformation_data = transformation_data

    def plot_stress_strain(self):
        plt.plot(self.stress_strain[:, 0], self.stress_strain[:, 1], self.color, lw=2)

    def plot_transformation(self):
        plt.plot(self.transformation_data[:, 0], self.transformation_data[:, 1], 'x' + self.color, ms=12, mew=2)


experiments = [Experiment(temperature=22, color='r'), Experiment(temperature=75, color='b'),
               Experiment(temperature=100, color='g'), Experiment(temperature=150, color='k')]

for experiment in experiments:
    plt.figure(0)
    experiment.plot_stress_strain()
    plt.figure(1)
    experiment.plot_transformation()

par = None
k = 0
a = 0
Ms = 0
for experiment in experiments[0:2]:
    trans_stress = np.interp(experiment.transformation_data[:, 0], experiment.stress_strain[:, 0],
                             experiment.stress_strain[:, 1])
    plt.figure(2)
    k = 0.010
    temp = experiment.temperature
    y = -np.log(1 - experiment.transformation_data[:, 1])
    if temp == 22.:
        par = np.polyfit(trans_stress, y, 1)
        print(par)
        print(par[1]/0.01)
    plt.plot(trans_stress, y, 'x' + experiment.color, ms=12, mew=2)
    if temp == 75.:
        new_p = trans_stress[0]*par[0] + par[1]
        k = (new_p - y)/(75-22)
        a = par[0]/k
        Ms = par[1]/k + 22
    print(k, a, Ms)

for experiment in experiments:
    y = k*(Ms + a*experiment.stress_strain[:, 1] - experiment.temperature)
    plt.figure(2)
    plt.plot(experiment.stress_strain[:, 1], y, experiment.color)
    plt.figure(1)
    fM = 1 - np.exp(-y)
    plt.plot(experiment.stress_strain[:, 0], fM, '--' + experiment.color)

plt.figure(1)
plt.xlim(0, 0.015)
plt.ylim(0.78, 1.)
plt.show()
