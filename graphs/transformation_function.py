import numpy as np
import matplotlib.pyplot as plt

K = 200E3/3/(1-2*.3)
G = 200E3/2/(1+0.3)
R1 = 0.02
R2 = 0.0
f0 = 0.65
b = 100.
Q = 2100.


def yield_function(st, martensite_fraction, plastic_strain):
    if not isinstance(martensite_fraction, np.ndarray):
        martensite_fraction = np.array([martensite_fraction])
    if not isinstance(plastic_strain, np.ndarray):
        plastic_strain = np.array([plastic_strain])
    s_dev = st - 1./3*sum(st)
    seqp = np.sqrt(1.5*np.sum(s_dev*s_dev))
    func = np.zeros((martensite_fraction.shape[0], plastic_strain.shape[0]))
    for i, phase in enumerate(martensite_fraction):
        for j, e in enumerate(plastic_strain):
            seq2 = seqp - 3*G*(e + R1*phase)/(1+3*G*R2*phase/1000)
            sy = 1000 + b*Q*e/(1 + b*e)
            func[i, j] = seq2 - sy

    return func


def transformation_function(st, martensite_fraction, plastic_strain):
    a1 = 0.056/3
    a2 = 0.028

    k = 0.017
    ms = 169
    mss = -107.87909071576405
    temp = 22

    if not isinstance(martensite_fraction, np.ndarray):
        martensite_fraction = np.array([martensite_fraction])
    if not isinstance(plastic_strain, np.ndarray):
        plastic_strain = np.array([plastic_strain])

    func = np.zeros((martensite_fraction.shape[0], plastic_strain.shape[0]))
    for i, phase in enumerate(martensite_fraction):
        for j, e in enumerate(plastic_strain):
            s_dev = st - 1./3*sum(st)
            seqp = np.sqrt(1.5*np.sum(s_dev*s_dev))
            se = seqp - 3*G*(e + R1*phase)/(1+3*G*R2*phase/1000)
            RA = R1+R2*se/1000
            print RA*phase
            nij = 3./2*s_dev/seqp
            sigma = st - 2*G*(e + RA*phase)*nij - K*phase*0.037/3
            print sigma
            m_stress = a1*np.sum(sigma, 0) + a2*se
            func[i, j] = 1 - np.exp(-k*(ms + m_stress + mss - temp)) - (phase + f0)
    return func


if __name__ == '__main__':
    fm = np.linspace(0, 0.1, 1000)
    dl = np.linspace(0, 0.01, 1000)
    # s = np.array([-27.55,     -27.55,       1095])
    s = np.array([2.557e-14,  1.145e-13,      701.1])
    # hdl = transformation_function(s, 0, dl)[0, :]
    hfm = transformation_function(s, fm, 0)
    # print (hdl[1] - hdl[0])/(dl[1] - dl[0])
    # print hdl.shape
    # plt.plot(dl, hdl, '*')
    print (hfm[1] - hfm[0])/(fm[1] - fm[0])

    # plt.figure(2)
    # ffm = yield_function(s, fm, 0)
    # print (ffm[1] - ffm[0])/(fm[1] - fm[0])
    # plt.plot(fm, ffm)
    # fdl = yield_function(s, 0, dl)[0, :]
    # plt.figure(3)
    # plt.plot(dl, fdl)
    # print (fdl[1] - fdl[0])/(dl[1] - dl[0])
    plt.show()
