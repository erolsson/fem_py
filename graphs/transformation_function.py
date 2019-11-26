import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from materials.transformation_materials import neu_sehitoglu

K = neu_sehitoglu.K
G = neu_sehitoglu.G
R1 = neu_sehitoglu.R1
R2 = neu_sehitoglu.R2
f0 = 0.65
b = neu_sehitoglu.b
Q = neu_sehitoglu.Q


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
            seq2 = seqp - 3*G*(e + R1*phase)/(1+3*G*R2*phase/neu_sehitoglu.sy0A)
            sy0 = (f0 + phase)*neu_sehitoglu.sy0M + (1 - (f0 + phase))*neu_sehitoglu.sy0A
            sy = sy0 + b*Q*e/(1 + b*e)
            func[i, j] = seq2 - sy

    return func


def transformation_function(st, martensite_fraction, plastic_strain):
    a1 = neu_sehitoglu.a1
    a2 = neu_sehitoglu.a2

    k = neu_sehitoglu.k
    ms = neu_sehitoglu.Ms
    mss = neu_sehitoglu.Mss
    temp = 22

    if not isinstance(martensite_fraction, np.ndarray):
        martensite_fraction = np.array([martensite_fraction])
    if not isinstance(plastic_strain, np.ndarray):
        plastic_strain = np.array([plastic_strain])

    func = np.zeros((martensite_fraction.shape[0], plastic_strain.shape[0]))
    for i, phase in enumerate(martensite_fraction):
        for j, e in enumerate(plastic_strain):
            st_dev = st - 1./3*sum(st)
            seqp = np.sqrt(1.5*np.sum(st_dev*st_dev))
            se = (seqp - 3*G*(e + R1*phase))/(1+3*G*R2*phase/neu_sehitoglu.sy0A)
            RA = R1+R2*se/neu_sehitoglu.sy0A
            nij = 3./2*st_dev/seqp
            sigma = st - 2*G*(e + RA*phase)*nij - K*phase*0.037/3
            sdev = sigma - 1./3*sum(sigma)
            s2 = sigma + np.array([0, 0, 1e-3])

            sdev2 = s2 - 1./3*sum(s2)
            se2 = np.sqrt(1.5*sum(sdev2*sdev2))
            m_stress = a1*np.sum(sigma, 0) + a2*se
            func[i, j] = 1 - np.exp(-k*(ms + m_stress + mss - temp)) - (phase + f0)
    return func


if __name__ == '__main__':
    DfM = 0.06
    DL = 0.001
    fm = np.linspace(-DfM, DfM, 100)
    dl = np.linspace(-DL, DL, 100)
    # s = np.array([-27.55,     -27.55,       1095])
    s = np.array([  0,     0,      1094.8674])
    # hdl = transformation_function(s, 0, dl)[0, :]
    hfm = transformation_function(s, fm, 0)
    hdl = transformation_function(s, 0, dl)[0]
    h = transformation_function(s, fm, dl)
    print (hdl[1] - hdl[0])/(dl[1] - dl[0])
    # print hdl.shape
    plt.plot(fm, hfm, '*')
    print hfm[0]
    print (hfm[1] - hfm[0])/(fm[1] - fm[0])

    plt.figure(2)
    ffm = yield_function(s, fm, 0)
    print (ffm[1] - ffm[0])/(fm[1] - fm[0])
    plt.plot(fm, ffm)
    fdl = yield_function(s, 0, dl)[0]
    f = yield_function(s, fm, dl)/1000
    idx = np.unravel_index(np.argmin(h**2 + f**2), f.shape)
    print idx
    print fm[idx[0]], dl[idx[1]]
    print np.min(h**2 + f**2)
    x, y = np.meshgrid(fm, dl)
    plt.figure(10)
    ax = plt.gca(projection='3d')
    ax.plot_surface(x, y, h**2 + f**2)
    plt.figure(3)
    plt.plot(dl, fdl)
    print (fdl[1] - fdl[0])/(dl[1] - dl[0])
    plt.show()
