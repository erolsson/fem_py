import numpy as np


# noinspection PyPep8Naming
class ElasticPlasticTransformMaterial:
    def __init__(self, E, v, sy0M, sy0A, Q, b, Cm, gamma_m, a, Ms, name, uniaxial_data, fM, beta, alpha, n, sde):
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
        self.R2 = 0.02

        self.dV = 0.037

        self.k = 0.017

        self.sde = sde

        self.back_stresses = Cm.shape[0]
        self.name = name

        self.Mss = (-1./self.k*np.log(1-uniaxial_data[0]) - self.Ms - uniaxial_data[2]*(self.a1 + self.a2 + 2*self.a3/27) +
                    uniaxial_data[1])

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
        if self.back_stresses > 0:
            return 3 + (self.gamma_m.shape[0]+1)*6
        return 3

    def umat_parameters(self):
        parameters = [self.E, self.v, self.sy0M, self.sy0A,  self.Q, self.b, self.gamma_m.shape[0]]
        kinematic_hardening_params = []
        for C, g in zip(self.Cm, self.gamma_m):
            kinematic_hardening_params += [C, g]
        return parameters + kinematic_hardening_params + [self.R1, self.R2, self.dV, self.Ms, self.Mss, self.k,
                                                          self.a1, self.a2, self.a3, self.beta, self.alpha, self.n,
                                                          self.sde]


test_material = ElasticPlasticTransformMaterial(E=200e3, v=0.3, sy0M=1000000., sy0A=485, Q=0*180., b=100.,
                                                Cm=np.array([135e3, 700e3, 50e3]),
                                                gamma_m=np.array([950., 500., 50.]), a=[0.056, 0.028, 0.],
                                                Ms=169, name='testMaterial', uniaxial_data=[0.8, 22., 485], fM=0.8,
                                                beta=0, alpha=4., n=4., sde=0.04)

neu_sehitoglu = ElasticPlasticTransformMaterial(E=200e3, v=0.3, sy0M=1000., sy0A=1000., Q=0*2100., b=100.,
                                                Cm=np.array([15432, 281622, 470894]),
                                                gamma_m=np.array([5., 236., 2301.]), a=np.array([0.056/3, 0.028, 0.]),
                                                Ms=169, name='NeuSehitoglu', uniaxial_data=[0.65, 22., 485], fM=0.65,
                                                beta=800, alpha=4., n=4., sde=0.04)

if __name__ == '__main__':
    print test_material.umat_parameters()
