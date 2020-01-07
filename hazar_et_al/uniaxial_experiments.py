from __future__ import print_function
import os

import numpy as np

from scipy.optimize import fmin
from scipy.stats import norm

import matplotlib.pyplot as plt

from fem_py.materials.transformation_materials import hazar_et_al

data_directory = os.path.expanduser('~/phase_transformations/hazar_et_al/')


class Experiment:
    def __init__(self, temperature, color, mode='tension'):
        self.temperature = temperature
        self.color = color
        self.mode = mode
        fig = 'a'
        if mode == 'compression':
            fig = 'b'
        print(mode)
        stress_strain = np.genfromtxt(data_directory + '/fig6' + fig + '_' + str(self.temperature) + 'C', delimiter=',')
        self.stress_strain = stress_strain[np.argsort(stress_strain[:, 0]), :]

        transformation_data = np.genfromtxt(data_directory + '/fig8' + fig + '_' + str(self.temperature) + 'C',
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
               Experiment(temperature=100, color='g'), Experiment(temperature=150, color='k'),
               Experiment(temperature=75, color='b', mode='compression')]


def uniaxial_stress_strain_curve_plastic(material, epl, fM=0.78):
    data = np.zeros((epl.shape[0]+1, 2))
    s = 0*epl + material.sy0(fM)
    s += material.Q*(1-np.exp(-material.b*epl))
    for C, g in zip(material.Cm, material.gamma_m):
        s += 2./3*C/g*(1-np.exp(-g*epl))
    s /= (1. + material.sde)
    # print(s)
    e = s/material.E + epl
    data[1:, 0] = e
    data[1:, 1] = s
    return data


def strain_transformation(par, epl, temperature):
    par = np.abs(par)
    fsb = 1 - (1 - par[0])*np.exp(-par[1]*epl)
    # print(fsb, epl, par)
    Gamma = par[3] - par[4]*np.array(temperature)/hazar_et_al.Ms
    p = norm.cdf(Gamma, loc=0, scale=abs(par[5]))
    c = (1 - 0.78)/np.exp(-par[2]*par[0]**4.*p)
    return 1 - c*np.exp(-par[2]*fsb**4.*p)


def platic_trans_residual(par, *data):
    epl, fm, temperature = data
    fm_model = strain_transformation(par, epl, temperature)
    return np.sum((fm - fm_model)**2)


ms_start = {22.: 800, 75: 1000, 100: 1290}


if __name__ == '__main__':
    plt.figure(0)
    plastic_strain = np.linspace(0, 0.01, 1000)
    transformation_free_data = uniaxial_stress_strain_curve_plastic(hazar_et_al, plastic_strain)
    plt.plot(transformation_free_data[:, 0], transformation_free_data[:, 1])
    all_epl = []
    all_fMe = []
    temp = []
    for experiment in experiments:
        plt.figure(0)
        experiment.plot_stress_strain()
        plt.figure(1)
        experiment.plot_transformation()
        plt.figure(2)
        experiment.plot_volume_expansion()
        plt.figure(3)
        if experiment.temperature == 75 and experiment.mode == 'tension':
            e_inel = experiment.stress_strain[:, 0] - experiment.stress_strain[:, 1]/205E3
            plt.plot(experiment.stress_strain[:, 0], e_inel)

    plt.show()
