import numpy as np
import matplotlib.pyplot as plt

K = 200E3/3/(1-2*.3)
G = 200E3/2/(1+0.3)
R1 = 0.02
R2 = 0.02
f0 = 0.65


def transformation_function(st, martensite_fraction, plastic_strain):
    a1 = 0.056/3
    a2 = 0.028

    k = 0.017
    ms = 169
    mss = -94.29909071576407
    temp = 22

    if not isinstance(martensite_fraction, np.ndarray):
        martensite_fraction = np.array([martensite_fraction])
    if not isinstance(plastic_strain, np.ndarray):
        plastic_strain = np.array([plastic_strain])

    func = np.zeros((martensite_fraction.shape[0], plastic_strain.shape[0]))
    for i, f in enumerate(martensite_fraction):
        for j, e in enumerate(plastic_strain):
            s_dev = st - 1./3*sum(st)
            seqp = np.sqrt(1.5*np.sum(s_dev*s_dev))
            se = seqp - 3*G*(e + R1*f)/(1+3*G*R2*f/1000)
            RA = R1+R2*se/1000
            sigma = st - 2*G*(e * RA*f) - K*f*0.037/3
            m_stress = a1*np.sum(sigma, 0) + a2*se
            func[i, j] = 1 - np.exp(-k*(ms + m_stress + mss - temp)) - (e + f0)
    return func


if __name__ == '__main__':
    fm = np.linspace(0, 0.1, 1000)
    s = np.array([0, 0, 701.1])
    print
    h = transformation_function(s, fm, 0)
    print h.shape
    plt.plot(fm, h)
    plt.show()
