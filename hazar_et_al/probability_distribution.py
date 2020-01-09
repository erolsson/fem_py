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

Md = 354.
Msigma = 22.


def beta(temperature, beta_0, loc, scale):
    x = (temperature-Msigma)/(Md - Msigma)
    return norm.cdf(-x, loc=loc, scale=scale)*beta_0


def beta_g(temperature, triax, beta0, g0, g1, g2):
    x = (temperature-Msigma)/(Md - Msigma)
    gamma = g0 - g1*x + g2*triax
    return norm.cdf(gamma)*beta0


def residual(par, *data):
    beta_0, g0, g1 = par
    temperature, val = data
    r = np.sum((val - beta(temperature, beta_0, g0, g1))**2)
    print(par, beta(temperature, beta_0, g0, g1), r)

    return np.sum(((val - beta(temperature, beta_0, g0, g1))/val)**2)


if __name__ == '__main__':
    temp = np.array([22, 75, 100, 150])
    beta_vals = np.array([1815.3888, 211, 125, 67.5])

    # par = fmin(residual, [1000, -0.65, 0.1], args=(temp, beta_vals), maxfun=1e6, maxiter=1e6)
    std = 0.2
    mean = 0
    x = (temp-Msigma)/(Md - Msigma)
    plt.plot(-x, beta_vals, 'x', lw=2, mew=2, ms=16)
    x_comp = -norm.ppf(0.1161132/1000, mean, std)
    print(x_comp - x[1])
    g2 = (x_comp - x[1])/2/std
    print(g2)
    g0 = -mean/std - g2
    g = g0 - 1./std*x + g2
    temp = np.linspace(-50, 350, 1000)
    x = (temp-Msigma)/(Md - Msigma)
    f = beta(temp, 1000, mean, std)

    plt.plot(-x, f)
    plt.figure(2)
    plt.plot(g, beta_vals, 'x', lw=2, mew=2, ms=16)
    g0 = 0 - g2

    g = g0 - 5.244*x + g2

    plt.plot(g, beta_g(temp, 1, 916, g0, 5.244, g2), lw=2, mew=2, ms=16)
    plt.plot(g0 - 5.244*(75 - Msigma)/(Md - Msigma) - g2, beta_g(75, -1, 916, g0, 5.244, g2), 'rx', lw=2, mew=2, ms=16)
    print(g0+g2, 1./std, g2)

    plt.show()


