import numpy as np
import matplotlib.pyplot as plt

K = 200E3/3/(1-2*.03)


def transformation_function(st, martensite_fraction):
    sigma = st - K*0.037/3*(martensite_fraction - 0.65)
    print sigma
    a1 = 0.056/3
    a2 = 0
    m_stress = a1*sigma + a2*sigma
    k = 0.017
    ms = 169
    mss = -94.29909071576407
    temp = 22
    return 1 - np.exp(-k*(ms + m_stress + mss - temp)) - martensite_fraction


if __name__ == '__main__':
    fm = np.linspace(0.6, 1.0, 1000)
    s = 701.1
    plt.plot(fm, transformation_function(s, fm))
    print (transformation_function(s, 0.6501) - transformation_function(s, 0.65))/0.0001
    plt.show()
