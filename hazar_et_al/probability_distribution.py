from __future__ import print_function

import numpy as np

from scipy.optimize import fmin
from scipy.stats import norm

import matplotlib
import matplotlib.pyplot as plt

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def beta(temperature, beta_0, loc, scale):
    return norm.cdf(-temperature, loc=loc, scale=scale)*beta_0


def residual(par, *data):
    beta_0, g0, g1 = par
    temperature, val = data
    print(par, beta(temperature, beta_0, g0, g1))
    return np.sum((val - beta(temperature, beta_0, g0, g1))**2)


if __name__ == '__main__':
    temp = np.array([75, 100, 150])
    beta_vals = np.array([211, 125, 67.5])

    # par = fmin(residual, [250, -150, 50], args=(temp, beta_vals), maxfun=1e6, maxiter=1e6)
    # x = par[1] + par[2]*(temp + 275.15)/(220+273.15)
    plt.plot(temp, beta_vals, 'x', lw=2, mew=2, ms=16)
    temp = np.linspace(22, 200, 1000)
    plt.plot(temp, beta(temp, 250, -110, 50))
    plt.show()


