import numpy as np
import matplotlib.pyplot as plt


def transformation_function(stress, martensite_fraction):
    a1 = 0.056/3
    a2 = 0
    m_stress = a1*stress + a2*stress
    k = 0.017
    ms = 169
    mss = -94.29909071576407
    temp = 22
    return 1 - np.exp(-k*(ms + m_stress + mss - temp)) - martensite_fraction


if __name__ == '__main__':
    fm = np.linspace(0.6, 1.0, 1000)
    s = 701.1
    plt.plot(fm, transformation_function(s, fm))
    plt.show()
