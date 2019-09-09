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


class TransformationMaterial:
    # noinspection PyPep8Naming
    J = np.array([[2./3,  -1./3, -1./3, 0,    0,    0],
                  [-1./3,  2./3, -1./3, 0,    0,    0],
                  [-1./3, -1./3,  2./3, 0,    0,    0],
                  [0,      0,     0,    1./2, 0,    0],
                  [0,      0,     0,    0,    1./2, 0],
                  [0,      0,     0,    0,    0,    1./2]])

    I1 = np.eye(6)
    I1[3:6, 3:6] /= 2

    E3 = np.zeros((6, 6))
    E3[0:3, 0:3] = np.ones((3, 3))

    I3 = np.array([1, 1, 1, 0, 0, 0])

    # noinspection PyPep8Naming
    def __init__(self, parameters, fM):
        # Elastic constants
        self.parameters = parameters

        # Ambient temperature
        self.T = 20.

        self.Mss = (-1./self.parameters.k*np.log(1-fM) - self.parameters.Ms -
                    self.parameters.s_lim*(self.parameters.a1 + self.parameters.a2 +
                                           2*self.parameters.a3/27) + self.T)

        # Elastic stiffness matrix
        self.D_el = 2*self.parameters.G*self.J + self.parameters.K*self.E3
        self.D_alg = np.copy(self.D_el)

        # state variables
        self._strain = np.zeros(6)
        self._stress = np.zeros(6)

        self._alpha_m = np.zeros((self.parameters.Cm.shape[0], 6))   # For kinematic hardening
        self._fM = fM                                 # Fraction of martensite
        self._R = 0                                   # Increase of yield surface
        self._ep_eff = 0.

    def stress(self):
        return self._stress

    # noinspection PyPep8Naming
    def update(self, strain):
        sy0 = self.parameters.sy0

        de = strain - self._strain

        # Calculating trial stress
        st = self._stress + np.dot(self.D_el, de)
        # Check if phase transformations occur
        phase_transformations = self._h(st, self._fM) > 0
        plastic = self._f(st, self.alpha(), sy0 + self._R) > 0
        elastic = not phase_transformations and not plastic
        if elastic:
            self._stress = np.copy(st)
            self.D_alg = np.copy(self.D_el)
        else:
            st_dev = np.copy(st)
            st_dev[0:3] = st[0:3] - np.sum(st[0:3])/3
            st_eq = np.sqrt(3./2*matrix_double_contract(st_dev - self.alpha(), st_dev - self.alpha()))

            s2 = None
            dDfM = 1e99
            DfM = 0

            DL = 0

            bij = np.zeros(6)
            s2_eq = None
            s2_dev = None
            RA = None

            G = self.parameters.G
            K = self.parameters.K
            R1 = self.parameters.R1
            R2 = self.parameters.R1

            Ms = self.parameters.Ms

            a1 = self.parameters.a1
            a2 = self.parameters.a2
            a3 = self.parameters.a3

            sy0_au = self.parameters.sy0_au
            dV = self.parameters.dV
            k = self.parameters.k
            nij = None
            Hinv = 0
            Y = 0
            A = np.copy(self.I1)
            if not plastic:
                # Newton - Rahpson algorithm for determining the increase in martensite fraction DfM
                while abs(dDfM) > 1e-14:
                    # Equivalent stress at time (2)
                    s2_eq = (st_eq - 3*G*R1*DfM)/(1+3*G*R2/sy0_au*DfM)

                    RA = R1 + R2*s2_eq/sy0_au

                    # Deviatoric stress at time (2)
                    s2_dev = (1-3*G*RA*DfM/st_eq)*(st_dev - self.alpha()) + self.alpha()

                    #  Stress at time 2
                    s2 = np.copy(s2_dev)
                    s2[0:3] += (np.sum(st[0:3]) - K*dV*DfM)/3

                    # B - tensor derivative of Mstress due to sigma_ij
                    F = k*np.exp(-k*(Ms + self._M_stress(s2) + self.Mss - self.T))
                    bij = F*(a1*self.I3 + a2*3./2*s2_dev/s2_eq +
                             a3*(matrix_contract(s2_dev, s2_dev) - 2./9*s2_eq**2*self.I3))

                    # Derivative of RA with respect to dFM
                    dRAdDfM = (3*G*R2*(R2*st_eq + R1*sy0_au)/(3*G*R2*DfM + sy0_au)**2)

                    ds2dDfM = -3*G*(RA/st_eq + DfM/st_eq*dRAdDfM)*(st_dev - self.alpha()) - K*dV*self.I3/3

                    dhddfM = F*matrix_double_contract(bij, ds2dDfM) - 1
                    dDfM = self._h(s2, self._fM + DfM)/dhddfM
                    DfM -= dDfM
                    print DfM, dDfM
                # Calculating the tangent
                s2_dev_tilde = s2_dev - self.alpha()
                nij = 1.5*s2_dev_tilde/s2_eq

                # A1 = (np.eye(6)
                #       + 9./2*G/s2_eq*DfM/st_eq*R2/sy0_au*np.outer(st_dev_tilde, s2_dev_tilde)
                #       + F*np.outer(3*G*RA/st_eq*st_dev_tilde + K/3*dV*self.I3, bij))

                # A2 = (2*G*(1-3*G*RA*DfM/st_eq)*self.J
                #       + 9./2*G**2*RA*DfM/st_eq**3*np.dot(np.outer(st_dev_tilde, st_dev_tilde), self.J)
                #       + K*self.E3)

                # self.D_alg = np.dot(np.linalg.inv(A1), A2)
                self._fM += DfM

            else:
                Q = self.parameters.Q
                b = self.parameters.b

                C = self.parameters.Cm
                gamma = self.parameters.gamma_m
                residual = 1e99
                alpha_2 = None
                R = 0
                iterations = 0
                residuals = []
                dsydDL = None
                R = self._R
                while residual > 1e-14:
                    iterations += 1
                    R = (self._R + b*Q*DL)/(1 + b*DL)
                    dsydDL = b*(Q - self._R)

                    sy = sy0 + R
                    RA = R1 + R2*sy/sy0_au
                    theta = 1./(1+gamma*DL)
                    Am = gamma/(1 + gamma*DL)**2
                    D = sy + 3*G*(DL+RA*DfM) + DL*np.sum(theta*C)
                    nij = 1.5*(st_dev - np.sum(np.outer(theta, np.ones(6))*self._alpha_m, 0))/D

                    s2_eq = np.sqrt(matrix_double_contract(nij*sy/1.5, nij*sy))
                    f = 2./3*matrix_double_contract(nij, nij) - 1
                    dCijdDL = np.sum(np.outer(Am, np.ones(6))*self._alpha_m, 0)
                    dDdDL = dsydDL + 3*G*(1+R2/sy0_au*dsydDL*DfM) + np.sum(theta**2*C)
                    dndDL = (1.5*dCijdDL - nij*dDdDL)/D
                    J11 = 4./3*matrix_double_contract(nij, dndDL)
                    if not phase_transformations:
                        dDL = f/J11
                        DL -= dDL
                        residual = np.abs(dDL)
                        residuals.append(residual)

                    else:
                        raise NotImplemented
                # Updating effective plastic strain
                self._ep_eff += DL
                H = 0

                # updating isostropic hardening
                self._R = R

                # Updating back -stress
                alpha_2 = np.copy(self._alpha_m)
                for i in range(self._alpha_m.shape[0]):
                    alpha_2[i, :] /= (1 + gamma[i]*DL)
                    alpha_2[i, :] += 2./3*(C[i]/(1 + gamma[i]*DL))*DL*nij
                self._alpha_m = np.copy(alpha_2)

                Am = gamma/(1 + gamma*DL)**2
                theta = 1./(1 + gamma*DL)

                H = np.sum(C) - matrix_double_contract(nij, np.sum(np.outer(gamma, np.ones(6))*self._alpha_m, 0))
                H += dsydDL
                Hinv = 1./H
                Y = s2_eq + np.sum(theta*C)*DL
                nnt = np.outer(nij, nij)
                A1 = 3*G/Y*(DL + RA*DfM)*self.J
                A2 = 2*G*(Hinv - (DL + RA*DfM)/Y*Hinv*(np.sum(theta*C) + dsydDL))*nnt
                A3 = 3*G*Hinv/Y*(DL + RA*DfM)*np.outer(np.sum(np.outer(Am/theta, np.ones(6))*self._alpha_m, 0), nij)
                # Updating the stress
                s2_dev = st_dev
                if DL > 0:
                    s2_dev = st_dev - 2*G*(DL + RA*DfM)*nij
                s2 = np.copy(s2_dev)
                s2[0:3] += (np.sum(st[0:3]) - K*dV*DfM)/3

                A += (A1 + A2 + A3)

            if DfM > 0:
                A4 = 2*G*RA*np.outer(nij, bij)
                A5 = K*dV/3*np.outer(self.I3, bij)
                A += A4 + A5
            self._stress = st_dev - 2*G*DL*nij - DfM*(RA*nij + 1./3*dV*self.I3)

            # self.D_alg = self.D_el
            # Calculating the tangent
            A *= 2
            A[0:3, 0:3] /= 2
            self.D_alg = np.dot(np.linalg.inv(A), self.D_el)

        # print self.D_alg[0:3, 0:3]
        self._strain = strain

    # noinspection PyPep8Naming
    def _h(self, s, fM):
        return 1. - np.exp(-self.parameters.k*(self.parameters.Ms + self._M_stress(s) + self.Mss - self.T)) - fM

    @staticmethod
    def _f(s, a, sy):
        sdev = np.copy(s)
        sdev[0:3] -= sum(sdev[0:3])/3
        return np.sqrt(3./2*matrix_double_contract(sdev - a, sdev - a)) - sy

    # noinspection PyPep8Naming
    def _M_stress(self, s):
        a1, a2, a3 = self.parameters.a1, self.parameters.a2, self.parameters.a3
        s_dev = np.copy(s)
        s_dev[0:3] = s_dev[0:3] - np.sum(s_dev[0:3])/3
        m_stress = (a1*np.sum(s[0:3])
                    + a2*np.sqrt(np.sum(s[0:3]**2) - s[0]*s[1] - s[0]*s[2] - s[1]*s[2] + 3*np.sum(s[3:6]**2))
                    + a3*vector_det(s_dev))
        return m_stress

    def tangent(self):
        return self.D_alg

    def alpha(self):
        return np.sum(self._alpha_m, 0)

    # noinspection PyPep8Naming
    def lambda_function(self, DL, st, DfM=0, ep_eff=0):
        gamma_i = self.parameters.gamma_i
        Ci = self.parameters.Ci
        if gamma_i != 0.:
            R = Ci/gamma_i*(1 - np.exp(-gamma_i*(ep_eff + DL)))
        else:
            R = Ci*(ep_eff + DL)
        sy = self.parameters.sy0 + R
        RA = self.parameters.R1 + self.parameters.R2*sy/self.parameters.sy0_au
        theta = 1./(1+self.parameters.gamma*DL)
        Am = self.parameters.gamma/(1 + self.parameters.gamma*DL)**2
        D = sy + 3*self.parameters.G*(DL+RA*DfM) + DL*np.sum(theta*self.parameters.C)
        st_dev = np.copy(st)
        st_dev -= np.sum(st_dev[0:3])/3
        nij = 1.5*(st_dev - np.sum(np.outer(theta, np.ones(6))*self._alpha_m, 0))/D
        return 2./3*matrix_double_contract(nij, nij) - 1


if __name__ == '__main__':
    pass
