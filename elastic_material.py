import numpy as np


class ElasticMaterial:
    # noinspection PyPep8Naming
    def __init__(self, E, v):
        self.E = E
        self.v = v

        C1 = self.E*(1-self.v)/(1-2*self.v)/(1+self.v)
        C2 = self.E*self.v/(1-2*self.v)/(1+self.v)
        C3 = self.E/2/(1+self.v)

        self.C = np.array([[C1, C2, C2, 0,  0,  0],
                           [C2, C1, C1, 0,  0,  0],
                           [C2, C2, C1, 0,  0,  0],
                           [0,  0,  0,  C3, 0,  0],
                           [0,  0,  0,  0,  C3, 0],
                           [0,  0,  0,  0,  0,  C3]])

    # noinspection PyPep8Naming
    def C_el(self):
        return self.C

    def tangent(self):
        return self.C

    def update(self, strain):
        pass
