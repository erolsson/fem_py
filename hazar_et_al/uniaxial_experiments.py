from __future__ import print_function
import os

import numpy as np

from scipy.optimize import fmin
from scipy.stats import norm

import matplotlib.pyplot as plt

from materials.transformation_materials import hazar_et_al

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


def uniaxial_stress_strain_curve_plastic(material, epl):
    data = np.zeros((epl.shape[0]+1, 2))
    s = 0*epl + material.sy0(0.78)
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


ms_start = {22.: 850, 75: 1000, 100: 1290}


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
        e = experiment.transformation_data[:, 0]
        s = np.interp(e, experiment.stress_strain[:, 0], experiment.stress_strain[:, 1])
        fM = experiment.transformation_data[:, 1]
        plt.figure(2)
        if experiment.temperature in ms_start:
            s = np.hstack(([ms_start[experiment.temperature]], s))
            fM = np.hstack(([0.78], fM))
        plt.plot(s, fM, 'x' + experiment.color, ms=12, mew=2)

        k = 0.01
        Ms = 220
        a = 0.0677272727272727
        Mss = -104.13818181818179
        fMsigma = 1 - np.exp(-k*(Ms + a*experiment.stress_strain[:, 1] + Mss - experiment.temperature))
        fMsigma[fMsigma <= 0.78] = 0.78
        plt.plot(experiment.stress_strain[:, 1], fMsigma, '--' + experiment.color)

        if experiment.temperature == 22.:
            e_tr_free = np.interp(experiment.stress_strain[:, 1],
                                  transformation_free_data[:, 1],
                                  transformation_free_data[:, 0])
            e_tr = experiment.stress_strain[:, 0] - e_tr_free
            plt.figure(3)
            x = experiment.stress_strain[1:, 1]/hazar_et_al.sy0A
            y = np.diff(e_tr)/np.diff(fMsigma) - hazar_et_al.dV/3
            x = x[~np.isinf(y)]
            y = y[(~np.isinf(y))]

            par_r = np.polyfit(x[2:], y[2:], 1)
            plt.plot(x, par_r[0]*x + par_r[1])
            print(par_r)
            par_r = [0.02/2**0.5, 0.02/2**0.5]
            plt.plot(x, y, '-*')


        fMep = fM[fM > 0.78] - np.interp(s[fM > 0.78], experiment.stress_strain[:, 1], fMsigma)
        dfm = (fM[fM > 0.78] - 0.78)
        s_eq = np.interp(experiment.transformation_data[:, 0], experiment.stress_strain[:, 0],
                         experiment.stress_strain[:, 1])
        e_tr = (par_r[1] + par_r[0]*s_eq/hazar_et_al.sy0A + hazar_et_al.dV/3)*dfm
        print(e_tr)
        fMep[fMep < 0] = 0
        e = experiment.transformation_data[:, 0]
        print(e)
        epl = e - e_tr - s_eq/hazar_et_al.E
        print(s_eq/hazar_et_al.E)
        plt.figure(5)
        plt.plot(epl, fMep + 0.78, 'x' + experiment.color, ms=12, mew=2)
        all_epl += epl.tolist()
        all_fMe += fMep.tolist()
        temp += [experiment.temperature]*len(fMep)
    all_epl = np.array(all_epl)
    all_fMe = np.array(all_fMe) + 0.78
    par_ep = fmin(platic_trans_residual, x0=[0.8, 10, 10, 60, 10, 1], args=(all_epl, all_fMe, temp), maxfun=1e6,
                  maxiter=1e6)
    plt.figure(5)
    epl = np.linspace(0, 0.002, 1000)
    # par = [0.5, 4., 3.]
    print(par_ep)
    plt.plot(epl, strain_transformation(par_ep, epl, 22), 'r')
    plt.plot(epl, strain_transformation(par_ep, epl, 75), 'b')
    plt.plot(epl, strain_transformation(par_ep, epl, 100), 'g')
    plt.plot(epl, strain_transformation(par_ep, epl, 150), 'k')
    plt.plot()
    plt.show()
