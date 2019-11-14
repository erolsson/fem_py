import numpy as np
import matplotlib.pyplot as plt

K = 200E3/3/(1-2*.3)
print K


def transformation_function(st, martensite_fraction):
    sigma = np.outer(0*martensite_fraction + 1, st) - K*0.037/3*np.outer((martensite_fraction - 0.65), np.ones(3))
    print sigma
    a1 = 0.056/3
    a2 = 0
    m_stress = a1*np.sum(sigma, 1)
    print "ms", m_stress
    k = 0.017
    ms = 169
    mss = -94.29909071576407
    temp = 22
    return 1 - np.exp(-k*(ms + m_stress + mss - temp)) - martensite_fraction


def d_transformation_function(st, martensite_fraction):
    sigma = np.outer(0*martensite_fraction + 1, st) - K*0.037/3*np.outer((martensite_fraction - 0.65), np.ones(3))
    print sigma
    a1 = 0.056/3
    a2 = 0
    m_stress = a1*np.sum(sigma, 1)
    print "ms", m_stress
    k = 0.017
    ms = 169
    mss = -94.29909071576407
    temp = 22
    print k*np.exp(-k*(ms + m_stress + mss - temp))*a1
    print -K*0.037/3
    return k*np.exp(-k*(ms + m_stress + mss - temp))*(-K*0.037*a1) - 1


if __name__ == '__main__':
    fm = np.linspace(0.65, 1.0, 1000)
    s = np.array([0, 0, 701.1])
    plt.plot(fm, transformation_function(s, fm), '*')
    print (transformation_function(s, 0.6501) - transformation_function(s, 0.65))/0.0001
    print d_transformation_function(s, 0.65)
    plt.show()
