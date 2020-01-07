from __future__ import print_function
import numbers

import numpy as np


# noinspection PyPep8Naming
class ElasticPlasticTransformMaterial:
    def __init__(self, E, v, sy0M, sy0A, Q, b, Cm, gamma_m, a, Ms, name, Mss, fM, beta, alpha, n, sde,
                 g0, g1, g2, g_mean, g_std, fsb0=0., yield_multi=1.):
        # Elastic parameters
        self.E = float(E)
        self.v = float(v)

        self.G = self.E/2/(1+self.v)
        self.K = self.E/3/(1-2*self.v)

        # Initial yield stress of Martensite and Austenite
        self.sy0M = sy0M
        self.sy0A = sy0A
        self.fM = fM
        # Parameters for isostropic hardening
        self.Q = Q
        self.b = b

        # parameters for Kinematic hardening
        self.Cm = Cm
        self.gamma_m = gamma_m

        # Martensite start temperature
        self.Ms = Ms

        # parameters for phase transformations
        self.a1 = a[0]
        self.a2 = a[1]
        self.a3 = a[2]

        # Parameters for the plastic strain transformations
        self.beta = beta
        self.alpha = alpha
        self.n = n

        self.R1 = 0.02
        self.R2 = 0.012

        self.dV = 0.037

        self.k = 0.017

        self.sde = sde

        self.g0 = g0
        self.g1 = g1
        self.g2 = g2

        self.g_mean = g_mean
        self.g_std = g_std

        self.fsb0 = fsb0

        self.back_stresses = Cm.shape[0]
        self.name = name
        self.yield_multi = 1.
        if isinstance(Mss, numbers.Real):
            self.Mss = Mss
        else:
            self.Mss = (-1./self.k*np.log(1 - Mss[0]) - self.Ms - Mss[2]*(self.a1 + self.a2) +
                        Mss[1])

    def abaqus_material_string(self):
        material_string = ['\t*Elastic',
                           '\t\t' + str(self.E) + ', ' + str(self.v)]
        sy0 = self.fM*self.sy0M + (1-self.fM)*self.sy0A
        if self.back_stresses > 0:
            material_string.append('\t*Plastic, hardening=COMBINED, datatype=PARAMETERS, number backstresses='
                                   + str(self.back_stresses))
            back_stress_string = '\t\t' + str(sy0)
            for i in range(self.back_stresses):
                back_stress_string += ', ' + str(self.Cm[i]) + ', ' + str(self.gamma_m[i])
            material_string.append(back_stress_string)
        if self.Q > 0:
            material_string.append('\t*Cyclic Hardening, parameters')
            material_string.append('\t\t' + str(sy0) + ', ' + str(self.Q) + ',' + str(self.b))
        return material_string

    def umat_depvar(self):
        return 7 + (self.back_stresses+1)*6

    def umat_parameters(self):
        k2 = (self.yield_multi*((1-self.fM)*self.sy0A + self.fM*self.sy0M) - (1-self.fM)*self.sy0A)/self.fM/self.sy0M
        parameters = [self.E, self.v, self.sy0M*k2, self.sy0A,  self.Q, self.b, self.gamma_m.shape[0]]
        kinematic_hardening_params = []
        for C, g in zip(self.Cm*self.yield_multi, self.gamma_m):
            kinematic_hardening_params += [C, g]
        return parameters + kinematic_hardening_params + [self.sde, self.R1, self.R2, self.dV, self.Ms, self.Mss,
                                                          self.k, self.a1, self.a2, self.a3, self.beta, self.alpha,
                                                          self.n, self.g0, self.g1, self.g2, self.g_mean, self.g_std]

    def sy0(self, fm):
        return fm*self.sy0M + (1 - fm)*self.sy0A


test_material = ElasticPlasticTransformMaterial(E=200e3, v=0.3, sy0M=1000000., sy0A=485, Q=0*180., b=100.,
                                                Cm=np.array([135e3, 700e3, 50e3]),
                                                gamma_m=np.array([950., 500., 50.]), a=[0.056, 0.028, 0.],
                                                Ms=169, name='testMaterial', Mss=[0.8, 22., 485], fM=0.8,
                                                beta=0, alpha=4., n=4., sde=0.04, g0=60, g1=176000, g2=5.2,
                                                g_mean=57, g_std=350)

neu_sehitoglu = ElasticPlasticTransformMaterial(E=203.3e3, v=0.3, sy0M=813., sy0A=420., Q=0*2100., b=100.,
                                                Cm=np.array([335485, 245783, 6853]),
                                                gamma_m=np.array([1016.5, 185, 0.]),
                                                a=np.array([0.056/3, 0.028, 0.]),
                                                Ms=169, name='NeuSehitoglu', Mss=[0.65, 22., 485], fM=0.65,
                                                beta=1., alpha=4., n=4., sde=0.04, g0=60, g1=176000, g2=5.2,
                                                g_mean=57, g_std=350)
print(neu_sehitoglu.Mss)

hazar_et_al = ElasticPlasticTransformMaterial(E=200.5e3, v=0.27, sy0M=1016, sy0A=420., Q=180., b=100.,
                                              Cm=np.array([135e3, 700e3, 50e3]),
                                              gamma_m=np.array([950., 500., 50.]),
                                              a=0.01*np.array([0, 0, 0.])/3,
                                              Ms=220, name='SKF', Mss=-48.194154926470915, fM=0.78,
                                              beta=282.0966412080984, alpha=141.507365886, n=4., sde=0.04,
                                              g0=0.,  # 3.077651873747456,
                                              # g1=68.83381914607745, g2=0, g_mean=0, g_std=1.,
                                              g1=0, g2=0, g_mean=0, g_std=1.,
                                              fsb0=0.12230665)
hazar_et_al.dV = 0.023364550877
hazar_et_al.R1 = 0.0154161628962
hazar_et_al.R2 = 0.0073451035740

SS2506 = ElasticPlasticTransformMaterial(E=200.e3, v=0.273, sy0M=1099.7293833, sy0A=420., Q=0., b=0.,
                                         Cm=np.array([3.45681604e+04, 1.65655117e+05, 5.11479494e+05]),
                                         gamma_m=np.array([4.17222046e+01, 3.05659804e+02, 7.58674760e+02]),
                                         a=0.04174*np.array([1., 0., 0.]),
                                         Ms=166.69903883197608, name='SS2506', Mss=-43.44, fM=0.8,
                                         beta=306.6016112399211, alpha=115.6731692918583, n=4., sde=0.0,
                                         g0=3.077651873747456,
                                         g1=68.83381914607745, g2=5, g_mean=0, g_std=29.5540022577844,
                                         fsb0=0.22758101605717712, yield_multi=1.)

SS2506.k = 0.0132
SS2506.R1 = 0.025
SS2506.R2 = 0.01
SS2506.dV = 0.037


"""
a1=0.0020086602368106056, 
g_std=29.5540022577844, 
beta=306.6016112399211, 
R1=0.018685671366484635, 
R2=0.005105251428246899, 
fsb0=0.22758101605717712, 
dV=0.022435285414492856, 
alpha=115.6731692918583, 
g0=115.6731692918583, 
g1=68.83381914607745,
"""
# Hazar et. al
# Temperature   beta
# 22            ??
# 75            4.92067433e+02

if __name__ == '__main__':
    print(neu_sehitoglu.umat_parameters())
    print(-np.log(1-0.78)/hazar_et_al.k - hazar_et_al.Ms - hazar_et_al.a1*800 + 22)
