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
                           [C2, C1, C2, 0,  0,  0],
                           [C2, C2, C1, 0,  0,  0],
                           [0,  0,  0,  C3, 0,  0],
                           [0,  0,  0,  0,  C3, 0],
                           [0,  0,  0,  0,  0,  C3]])
        print self.C
        self._stress = np.zeros(6)

    # noinspection PyPep8Naming
    def tangent(self):
        return self.C

    def update(self, strain):
        self._stress = np.dot(self.C, strain)

    def stress(self):
        return self._stress
