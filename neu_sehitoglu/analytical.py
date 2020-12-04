import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from scipy.optimize import fmin

from materials.transformation_materials import neu_sehitoglu
from experiments import experiments

matplotlib.style.use('classic')
matplotlib.use("Qt5agg")
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{gensymb}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    plt.figure(0)
    Mss = -45

    austenite = 0.35 - experiments.torsion.dv/neu_sehitoglu.dV
    plt.figure(0)
    plt.plot(experiments.torsion.stress, austenite)

    plt.figure(1)
    plt.plot(experiments.torsion.stress**2, -np.log(austenite)/neu_sehitoglu.k - Mss + 22)
    [a, Ms] = np.polyfit(experiments.torsion.stress**2, -np.log(austenite)/neu_sehitoglu.k - Mss + 22, 1)
    a2 = a
    print(a, Ms)
    plt.plot(experiments.torsion.stress**2, a2*experiments.torsion.stress**2 + Ms, '--')

    fa = np.exp(-neu_sehitoglu.k*(Ms + a2*experiments.torsion.stress**2 + Mss - 22))
    plt.figure(0)
    plt.plot(experiments.torsion.stress, fa, '--')

    austenite_tension = 0.35 - experiments.tension.dv/neu_sehitoglu.dV
    austenite_compression = 0.35 - experiments.compression.dv/neu_sehitoglu.dV
    plt.figure(2)
    plt.plot(experiments.tension.stress, austenite_tension)
    plt.plot(experiments.compression.stress, austenite_compression)

    plt.figure(3)
    y_tension = -np.log(austenite_tension)/neu_sehitoglu.k - Mss + 22 - Ms - a2*experiments.tension.stress**2/3
    y_comp = -np.log(austenite_compression)/neu_sehitoglu.k - Mss + 22 - Ms - a2*experiments.compression.stress**2/3
    plt.plot(experiments.tension.stress, y_tension)
    plt.plot(experiments.compression.stress, y_comp)

    def fun(par, x):
        return par[0]*x + par[1]*x**3

    def residual(par, x, y):
        return np.sum((y - fun(par, x))**2)
    [a1, a3] = fmin(residual, [0.01, 2e-9], args=(experiments.tension.stress, y_tension))
    a3 *= 27./2.
    s = experiments.tension.stress
    plt.plot(s, a1*s + 2*a3*s**3/27)
    plt.plot(s, -a1*s + -2*a3*s**3/27)
    m_stress = a1*s + a2*s**2 + 2*a3*s**3/27
    fa = np.exp(-neu_sehitoglu.k*(Ms + m_stress + Mss - 22.))
    plt.figure(2)
    plt.plot(s, fa)

    s = experiments.compression.stress
    m_stress = -a1*s + a2*s**2/3 + -2*a3*s**3/27
    fa = np.exp(-neu_sehitoglu.k*(Ms + m_stress + Mss - 22.))
    plt.plot(s, fa)
    print(a1, a2, a3)
    plt.show()


if __name__ == '__main__':
    main()
