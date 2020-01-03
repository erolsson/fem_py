from __future__ import print_function
import os

import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt


def kinematic_hardening(strain, par):
    back_stress = 0*strain
    for i in range(len(par)/2):
        C = abs(par[2*i])
        g = abs(par[2*i + 1])
        back_stress += C/g*(1 - np.exp(-g*strain))
    return back_stress


def kinematic_params_residual(params, *data):
    epl_exp, back_stress_exp = data
    back_stress = kinematic_hardening(epl_exp, params)
    return np.sum((back_stress - back_stress_exp)**2)


if __name__ == '__main__':
    minimum_yield_stress = 300
    yield_stress_offset = 0.0001
    starting_params = [15432, 5., 281622, 236, 470894, 2301]

    uniaxial_data_comp = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig2_compression'),
                                       delimiter=',')
    uniaxial_data_comp = uniaxial_data_comp[uniaxial_data_comp[:, 0].argsort(), :]

    e_exp_comp = uniaxial_data_comp[:, 0]
    s_exp_comp = uniaxial_data_comp[:, 1]

    plt.figure(0)
    plt.plot(e_exp_comp, s_exp_comp, '*')

    E = np.polyfit(e_exp_comp[s_exp_comp <= minimum_yield_stress], s_exp_comp[s_exp_comp <= minimum_yield_stress], 1)[0]
    print(E)
    epl = e_exp_comp - s_exp_comp/E
    sy0 = np.interp(yield_stress_offset, epl, s_exp_comp)
    sy = s_exp_comp[epl > yield_stress_offset]
    epl = epl[epl > yield_stress_offset]
    plt.figure(1)
    plt.plot(epl, sy)

    params = fmin(kinematic_params_residual, starting_params, args=(epl, sy - sy0), maxfun=1e6, maxiter=1e6)
    alpha = kinematic_hardening(epl, params)
    plt.plot(epl, sy0 + alpha)
    print(sy0, params)
    plt.show()
