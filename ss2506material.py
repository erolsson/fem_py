import numpy as np


def matrix_double_contract(a, b):
    return np.sum(a*b) + np.sum(a[3:]*b[3:])


def vector_det(a):
    return a[0]*(a[1]*a[2] - a[5]**2) - a[3]*(a[3]*a[2]-a[5]*a[4]) + a[4]*(a[1]*a[4]-a[3]*a[5])


def matrix_contract(a, b):
    return np.array([a[0]*b[0] + a[3]*b[3] + a[4]*b[4],
                     a[1]*b[1] + a[3]*b[3] + a[5]*b[5],
                     a[2]*b[2] + a[4]*b[4] + a[5]*b[5],
                     a[0]*b[3] + a[3]*b[1] + a[4]*b[5],
                     a[0]*b[4] + a[4]*b[2] + a[3]*b[5],
                     a[1]*b[5] + a[3]*b[4] + a[5]*b[2]])


class SS2506:
    # noinspection PyPep8Naming
    J = np.array([[2./3,  -1./3, -1./3, 0,    0,    0],
                  [-1./3,  2./3, -1./3, 0,    0,    0],
                  [-1./3, -1./3,  2./3, 0,    0,    0],
                  [0,      0,     0,    1./2, 0,    0],
                  [0,      0,     0,    0,    1./2, 0],
                  [0,      0,     0,    0,    0,    1./2]])

    E3 = np.zeros((6, 6))
    E3[0:3, 0:3] = np.ones((3, 3))

    I3 = np.array([1, 1, 1, 0, 0, 0])

    # noinspection PyPep8Naming
    def __init__(self, E, v, a, dV, R_params, sy0_au, fM):
        # Elastic constants
        self.G = E/2./(1+v)
        self.K = E/3./(1-2*v)

        # Koistinen - Marburger parameter
        self.k = 0.04

        # Martensite start temperature
        self.Ms = 169

        # Ambient temperature
        self.T = 20

        # Martensite transformation parameters
        self.a1 = a[0]
        self.a2 = a[1]
        self.a3 = a[2]

        s_lim = 485.
        self.Mss = -1./self.k*np.log(1-fM) - self.Ms - s_lim*(self.a1 + self.a2 + 2*self.a3/27) + self.T

        # Martensite expansion
        self.dV = dV

        # Austenite plasticity parameters
        self.R1 = R_params[0]
        self.R2 = R_params[1]

        # Initial yield stress of austenite
        self.sy0_au = sy0_au

        # Elastic stiffness matrix
        self.D_el = 2*self.G*self.J + self.K*self.E3
        self.D_alg = np.copy(self.D_el)
        # state variables
        self._strain = np.zeros(6)
        self._stress = np.zeros(6)

        self._alpha = np.zeros(6)   # For kinematic hardening
        self._fM = fM               # Fraction of martensite

    def stress(self):
        return self._stress

    # noinspection PyPep8Naming
    def update(self, strain):
        de = strain - self._strain

        # Calculating trial stress
        st = self._stress + np.dot(self.D_el, de)
        # Check if phase transformations occur
        phase_transformations = self._h(st, self._fM) > 0

        plastic = False
        elastic = not phase_transformations and not plastic

        if elastic:
            self._stress = st
            self.D_alg = self.D_el
        else:
            st_dev = np.copy(st)
            st_dev[0:3] = st[0:3] - np.sum(st[0:3])/3
            st_eq = np.sqrt(3./2*matrix_double_contract(st_dev - self._alpha, st_dev - self._alpha))
            s2 = st

            if not plastic:
                # Newton - Rahpson algorithm for determining the increase in martensite fraction DfM
                dDfM = 1e99
                DfM = 0
                while abs(dDfM) > 1e-9:
                    # Equivalent stress at time (2)
                    s2_eq = (st_eq - 3*self.G*self.R1*DfM)/(1+3*self.G*self.R2/self.sy0_au*DfM)

                    RA = self.R1 + self.R2*s2_eq/self.sy0_au

                    # Deviatoric stress at time (2)
                    s2_dev = (1-3*self.G*RA*DfM/st_eq)*(st_dev - self._alpha) + self._alpha

                    #  Stress at time 2
                    s2 = np.copy(s2_dev)
                    s2[0:3] += (np.sum(st[0:3]) - self.K*self.dV*DfM)/3

                    # B - tensor derivative of Mstress due to sigma_ij
                    F = self.k*np.exp(-self.k*(self.Ms + self._M_stress(s2) + self.Mss - self.T))
                    B = (self.a1*self.I3
                         + self.a2*3./2*s2_dev/s2_eq
                         + self.a3*(matrix_contract(s2_dev, s2_dev) - 2./9*s2_eq*self.I3))

                    # Derivative of RA with respect to dFM
                    dRAdDfM = (3*self.G*self.R2*(self.R2*st_eq + self.R1*self.sy0_au)
                               / (3*self.G*self.R2*DfM + self.sy0_au)**2)

                    ds2dDfM = -3*self.G*(RA/st_eq + DfM/st_eq*dRAdDfM)*(st_dev - self._alpha) - self.K*self.dV*self.I3/3

                    dhddfM = F*matrix_double_contract(B, ds2dDfM) - 1
                    dDfM = self._h(s2, self._fM + DfM)/dhddfM
                    DfM -= dDfM

                self._fM += DfM
                self._stress = s2

            else:
                print st
                raise NotImplemented

        self._strain = strain
    # noinspection PyPep8Naming

    def _h(self, s, fM):
        return 1. - np.exp(-self.k*(self.Ms + self._M_stress(s) + self.Mss - self.T)) - fM

    # noinspection PyPep8Naming
    def _M_stress(self, s):
        s_dev = np.copy(s)
        s_dev[0:3] = s_dev[0:3] - np.sum(s_dev[0:3])/3
        m_stress = (self.a1*np.sum(s[0:3])
                    + self.a2*np.sqrt(np.sum(s[0:3]**2) - s[0]*s[1] - s[0]*s[2] - s[1]*s[2] + 3*np.sum(s[3:6]**2))
                    + self.a3*vector_det(s_dev))
        return m_stress

    def tangent(self):
        return self.D_alg


if __name__ == '__main__':
    material = SS2506(E=200e3, v=0.3, a=[0.056, 0.028, 0.], dV=0.037, R_params=[0.02, 0.02], sy0_au=400, fM=0.8)
    exx = 0.003
    material.update(np.array([exx, 0, 0, 0, 0, 0]))
