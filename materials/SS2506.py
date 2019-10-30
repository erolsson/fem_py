import numpy as np


# noinspection PyPep8Naming
class ElasticPlasticTransformMaterial:
    def __init__(self, E, v, sy0M, sy0A, Q, b, Cm, gamma_m, a):
        # Elastic parameters
        self.E = float(E)
        self.v = float(v)

        self.G = self.E/2/(1+self.v)
        self.K = self.E/3/(1-2*self.v)

        # Initial yield stress of Martensite and Austenite
        self.sy0M = sy0M
        self.sy0A = sy0A
        self.fM = 0.8
        # Parameters for isostropic hardening
        self.Q = Q
        self.b = b

        # parameters for Kinematic hardening
        self.Cm = Cm
        self.gamma_m = gamma_m

        # Martensite start temperature
        self.Ms = 169

        # parameters for phase transformations
        self.a1 = a[0]
        self.a2 = a[1]
        self.a3 = a[2]

        self.R1 = 0.0
        self.R2 = 0.0

        self.dV = 0.037

        self.k = 0.04
        self.s_lim = 485.

    def umat_depvar(self):
        if self.gamma_m.shape[0] > 0:
            return 3 + (self.gamma_m.shape[0]+1)*6
        return 3

    def umat_parameters(self):
        kinematic_hardening_params = []
        for C, g in zip(self.Cm, self.gamma_m):
            kinematic_hardening_params += [C, g]
        return [self.E, self.v, self.sy0M, self.sy0M,  self.Q, self.b, self.gamma_m.shape[0]] + \
            kinematic_hardening_params + [self.R1, self.R2]


test_material = ElasticPlasticTransformMaterial(E=200e3, v=0.3, sy0M=1000., sy0A=485, Q=0*180., b=100.,
                                                Cm=np.array([135e3, 700e3, 50e3]),
                                                gamma_m=0*np.array([950., 500., 50.]), a=[0.056, 0.028, 0.])

if __name__ == '__main__':
    print test_material.umat_parameters()
