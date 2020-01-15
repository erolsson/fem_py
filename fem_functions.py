import numpy as np

from elastic_material import ElasticMaterial


# noinspection PyPep8Naming
class GaussPoint:
    def __init__(self, coordinates, B, detJ, material):
        self.coordinates = coordinates
        self.B = B
        self.detJ = detJ
        self.material = material


class Node:
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)


class Element8:
    def __init__(self, nodes, material, material_parameters):
        self._pos = np.array([[-1, -1, -1],
                              [1,  -1, -1],
                              [1,   1, -1],
                              [-1,  1, -1],
                              [-1, -1,  1],
                              [1,  -1,  1],
                              [1,   1,  1],
                              [-1,  1,  1]])

        self._nodes = nodes
        self._gauss_points = []
        self._K = np.zeros((24, 24))
        self._fint = np.zeros(24)
        for i in [-1, 1]:
            for j in [-1, 1]:
                for k in [-1, 1]:
                    coordinates = (i/np.sqrt(3), j/np.sqrt(3), k/np.sqrt(3))
                    b_matrix = self.B(*coordinates)
                    jacobian = np.linalg.det(self.J(*coordinates))
                    self._gauss_points.append(GaussPoint(coordinates=coordinates, B=b_matrix, detJ=jacobian,
                                                         material=material(**material_parameters)))
        self.update(np.zeros(24))

    # noinspection PyPep8Naming
    def B(self, xi, eta, zeta):
        B_matrix = np.zeros((6, 24))
        jacobian = self.J(xi, eta, zeta)
        d = self.d_matrix(xi, eta, zeta)
        for i in range(8):
            dx = np.linalg.solve(jacobian, d[:, i])
            for j in range(3):
                B_matrix[j, 3*i+j] = dx[j]
            B_matrix[3, 3*i] = dx[1]
            B_matrix[3, 3*i+1] = dx[0]

            B_matrix[4, 3*i] = dx[2]
            B_matrix[4, 3*i+2] = dx[0]

            B_matrix[5, 3*i + 1] = dx[2]
            B_matrix[5, 3*i + 2] = dx[1]

        return B_matrix

    def d_matrix(self, xi, eta, zeta):
        d = np.zeros((3, 8))
        for i in range(8):
            d[0, i] = (1. + eta*self._pos[i, 1])*(1. + zeta*self._pos[i, 2])*self._pos[i, 0]/8
            d[1, i] = (1. + xi*self._pos[i, 0])*(1. + zeta*self._pos[i, 2])*self._pos[i, 1]/8
            d[2, i] = (1. + xi*self._pos[i, 0])*(1. + eta*self._pos[i, 1])*self._pos[i, 2]/8
        return d

    # noinspection PyPep8Naming
    def J(self, xi, eta, zeta):
        pos_matrix = np.zeros((8, 3))
        for i in range(3):
            pos_matrix[:, i] = [n.coordinates[i] for n in self._nodes]
        return np.dot(self.d_matrix(xi, eta, zeta), pos_matrix)

    def update(self, nodal_displacements):
        self._K *= 0
        self._fint = 0
        B21 = np.zeros(3)
        B24 = np.zeros(3)
        for i, gp in enumerate(self._gauss_points):
            strain = np.dot(gp.B, nodal_displacements)
            gp.material.update(strain)
            B21 += gp.B[[2, 4, 5], 20]**2
            B24 += gp.B[[2, 4, 5], 23]**2
            self._K += np.dot(gp.B.transpose(), np.dot(gp.material.tangent(), gp.B))*gp.detJ
            self._fint += np.dot(gp.B.transpose(), gp.material.stress())*gp.detJ

    # noinspection PyPep8Naming
    def Kt(self):
        return self._K

    def f_int(self):
        return self._fint


if __name__ == '__main__':
    node_list = [Node(coordinates=[0, 0, 0]),
                 Node(coordinates=[1, 0, 0]),
                 Node(coordinates=[1, 1, 0]),
                 Node(coordinates=[0, 1, 0]),
                 Node(coordinates=[0, 0, 1]),
                 Node(coordinates=[1, 0, 1]),
                 Node(coordinates=[1, 1, 1]),
                 Node(coordinates=[0, 1, 1])]

    element = Element8(nodes=node_list, material=ElasticMaterial, material_parameters={'E': 1., 'v': 0.0})
    print element.J(0, 0, 0)

