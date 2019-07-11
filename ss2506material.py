import numpy as np


def matrix_double_contract(a, b):
    return np.sum(a*b) + np.sum(a[3:]*b[3:])


def vector_det(a):
    return a[0]*(a[1]*a[2] - a[5]**2) - a[3]*(a[3]*a[2]-a[5]*a[4]) + a[4]*(a[1]*a[4]-a[3]*a[5])


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
        phase_transformations = self._h(st) > 0

        plastic = False
        elastic = not phase_transformations and not plastic

        if elastic:
            self._stress = st
            self.D_alg = self.D_el
        else:
            st_dev = st - np.sum(st[0:3])/3
            st_eq = np.sqrt(3./2*matrix_double_contract(st_dev - self._alpha, st_dev - self._alpha))
            s2 = st

            B =

            if not plastic:
                # Newton - Rahpson algorithm for determining the increase in martensite fraction dfM
                residual = 1e99
                dfM = 0
                while residual > 1e-9:
                    s_eq_2 = (st_eq - 3*self.G*self.R1*dfM)/(1+3*self.G*self.R2/self.sy0_au*dfM)
                    RA = self.R1 + self.R2*s_eq_2/self.sy0_au
                    s2 = (1-3*self.G*RA*dfM/st_eq)*(st_dev - self._alpha) - self._alpha
                    s2[0:3] += self.K*self.dV*dfM

                    dhddfm = self.k*np.exp(-self.k*(self.Ms + self._M_stress(s2) + self.Mss - self.T)) - 1

                    dfM -= self._h(s2)/dhddfm
                    residual = dfM/(dfM + self._fM)

                self._fM += dfM
                self._stress = s2

            else:
                print st
                raise NotImplemented

        self._strain = strain
    # noinspection PyPep8Naming

    def _h(self, s):
        return 1. - np.exp(-self.k*(self.Ms + self._M_stress(s) + self.Mss - self.T)) - self._fM

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