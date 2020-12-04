import matplotlib
import matplotlib.pyplot as plt

import numpy as np

from scipy.optimize import fmin

from materials.transformation_materials import SS2506

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    ra = 0.18
    st = 600
    stress_strain = np.genfromtxt('experimental_data/tension_data_case_EN.dat', delimiter=',')
    strains = np.genfromtxt('experimental_data/transversal_strain_tension_EN.dat', delimiter=',')

    experimental_data = np.zeros((stress_strain.shape[0], 3))
    experimental_data[:, 0:2] = stress_strain
    experimental_data[:, 2] = strains[:, 1]
    plt.figure(0)
    dv = experimental_data[:, 0] + 2*experimental_data[:, 2]
    s = experimental_data[:, 1]
    austenite = ra - (dv - experimental_data[:, 1]/SS2506.E*(1-2*SS2506.v))/SS2506.dV
    plt.plot(s, austenite)

    def fa(par, x):
        par[0] = 0.02
        par[1] = 2e-4
        par[2] = 0e-7
        # par[3] = -148.38998859637542
        m_stress = par[0]*x + par[1]*x**2/3 + par[2]*2/27*x**3
        return np.exp(-SS2506.k*(SS2506.Ms + m_stress + par[3] - 22))

    def residual(par, x, y):
        return np.sum((fa(par, x)*1000 - y*1000)**2)

    [a1, a2, a3, Mss] = fmin(residual, [0.06, 4e-4, 2e-7, -50], args=(s[s > st], austenite[s > st]), maxfun=1e6)
    print(a1, a2, a3, Mss)
    # Mss = -61.067138726239385
    plt.plot(s, fa([a1, a2, a3, Mss], s))
    plt.figure(3)
    plt.plot(s, a1*s + a2*s**2/3 + a3*2/27*s**3 + Mss)
    plt.plot(s, 4e-4*s**2/3)
    plt.figure(5)
    martensite = 1 - austenite
    et = experimental_data[:, 0] - experimental_data[:, 1]/SS2506.E
    detdfm = (np.diff(et)/np.diff(martensite) - SS2506.dV/3)*3./2
    plt.plot(s[1:][s[1:] > st], detdfm[s[1:] > st])
    y = detdfm[s[1:] > st]
    x = s[1:][s[1:] > st]
    print(np.polyfit(x, y, 1))

    stress_strain = np.genfromtxt('experimental_data/compression_data_case_EN.dat', delimiter=',')
    strains = np.genfromtxt('experimental_data/transversal_strain_compression_EN.dat', delimiter=',')

    experimental_data = np.zeros((stress_strain.shape[0], 3))
    experimental_data[:, 0:2] = stress_strain
    experimental_data[:, 2] = strains[:, 1]
    # SS2506.dV = 0.037
    plt.figure(2)
    dv = experimental_data[:, 0] + 2*experimental_data[:, 2]
    s = experimental_data[:, 1]
    austenite = ra - (dv - experimental_data[:, 1]/205e3*(1-2*SS2506.v))/SS2506.dV
    plt.plot(s, austenite)

    plt.show()


if __name__ == '__main__':
    main()
