from __future__ import print_function, division
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

        self.M_sigma = 22
        self.M_d = 383.585

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
                                                          self.n, self.g0, self.g1, self.g2, self.g_mean, self.g_std,
                                                          self.M_sigma, self.M_d]

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
                                                a=np.array([0.03501554463596976, 0.013017505521745059, 0.]),
                                                Ms=129.1, name='NeuSehitoglu', Mss=-112.59175344345212, fM=0.65,
                                                beta=0., alpha=4., n=4., sde=0.04, g0=60, g1=176000, g2=5.2,
                                                g_mean=57, g_std=350)
neu_sehitoglu.R1 = 0.02013530291125219
neu_sehitoglu.R2 = 0.02013530291125219
neu_sehitoglu.dV = 0.037438411828539964
neu_sehitoglu.k = 0.008300252500000004


hazar_et_al = ElasticPlasticTransformMaterial(E=200.5e3, v=0.27, sy0M=1016, sy0A=420., Q=180., b=100.,
                                              Cm=np.array([135e3, 700e3, 50e3]),
                                              gamma_m=np.array([950., 500., 50.]),
                                              a=0.0084967*np.array([1, 0, 0.]),
                                              Ms=220., name='SKF', Mss=-114.08, fM=0.78,
                                              beta=841.893, alpha=129.5, n=4., sde=0.04,
                                              g0=-1.918/2,  # 3.077651873747456,
                                              # g1=68.83381914607745, g2=0, g_mean=0, g_std=1.,
                                              g1=5.18, g2=1.918/2., g_mean=0, g_std=1.,
                                              fsb0=0.12948)
                                              # fsb0=0.23 )
hazar_et_al.dV = 0.023364550877
hazar_et_al.R1 = 0.0198377
hazar_et_al.R2 = 0.00533789

SS2506 = ElasticPlasticTransformMaterial(E=195.e3, v=0.27, sy0M=1137.77878007, sy0A=420., Q=0., b=0.,
                                         Cm=np.array([43496.06492213, 438169.64328922, 542518.46435554]),
                                         gamma_m=np.array([52.29850275, 483.57664975, 2541.45441662]),
                                         a=np.array([0.009337669592, 2*2e-5,
                                                     81/16*2e-7]),
                                         Ms=197.75201617791663, name='SS2506', Mss=-148.38998859637542, fM=0.8,
                                         beta=0*841.893, alpha=129.5, n=4., sde=0.04,
                                         g0=-1.918/2,  # 3.077651873747456,
                                         # g1=68.83381914607745, g2=0, g_mean=0, g_std=1.,
                                         g1=5.18, g2=1.918/2., g_mean=0, g_std=1.,
                                         fsb0=0.12948)

SS2506.k = 0.014793663333333332
SS2506.R1 = 8.76558745e-03
SS2506.R2 = 2.11646233e-02
SS2506.dV = 0.030435492416000003
# SS2506.dV = 0.037


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
    ms = -np.log(0.2)/SS2506.k - SS2506.Ms --53.40614314338714 + 22
    print(ms, SS2506.dV/3)
    print(SS2506.a1, SS2506.a2, SS2506.a3)
