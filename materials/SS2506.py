import numpy as np


# noinspection PyPep8Naming
class ElasticPlasticTransformMaterial:
    def __init__(self, E, v, sy0, Q, b, Cm, gamma_m, a):
        # Elastic parameters
        self.E = float(E)
        self.v = float(v)

        self.G = self.E/2/(1+self.v)
        self.K = self.E/3/(1-2*self.v)

        # Initial yield stress
        self.sy0 = sy0

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

        self.R1 = 0.02
        self.R2 = 0.02

        self.dV = 0.037
        self.sy0_au = 400

        self.k = 0.04
        self.s_lim = 485.


test_material = ElasticPlasticTransformMaterial(E=200e3, v=0.3, sy0=600, Q=200, b=100, Cm=np.array([10000.]),
                                                gamma_m=np.array([1.]), a=[0., 0., 0.])
