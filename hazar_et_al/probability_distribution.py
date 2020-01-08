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

Md = 250


def beta(temperature, beta_0, loc, scale):
    x = (Md-50)/(Md - temperature)
    return norm.cdf(-x, loc=loc, scale=scale)*beta_0


def beta_g(temperature, triax, beta0, g0, g1, g2):
    x = (Md - 50)/(Md - temperature)
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
    beta_vals = np.array([741.70, 211, 125, 67.5])

    # par = fmin(residual, [1000, -0.65, 0.1], args=(temp, beta_vals), maxfun=1e6, maxiter=1e6)
    std = 0.5
    mean = -1
    x = (Md - 50.)/(Md - temp)
    plt.plot(-x, beta_vals, 'x', lw=2, mew=2, ms=16)
    x_comp = -norm.ppf(0.1161132/1000, mean, std)
    print(x_comp - x[1])
    g2 = (x_comp - x[1])/2/std
    print(g2)
    g0 = -mean/std - g2
    g = g0 - 1./std*x + g2
    temp = np.linspace(-50, 150, 1000)
    x = (Md - 50)/(Md - temp)
    f = beta(temp, 1000, mean, std)

    plt.plot(-x, f)
    plt.figure(2)
    plt.plot(g, beta_vals, 'x', lw=2, mew=2, ms=16)
    g = g0 - x/std + g2
    plt.plot(g, beta_g(temp, 1, 1000, g0, 1/std, g2), lw=2, mew=2, ms=16)
    plt.plot(g0 - (Md - 50)/(Md - 75)/std - g2, beta_g(75, -1, 1000, g0, 1/std, g2), 'rx', lw=2, mew=2, ms=16)
    print(g0+g2, 1./std, g2)

    plt.show()


