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
    def __init__(self, nodes, material):
        self.pos = np.array([[-1, -1, -1],
                             [1,  -1, -1],
                             [1,   1, -1],
                             [-1,  1, -1],
                             [-1, -1,  1],
                             [1,  -1,  1],
                             [1,   1,  1],
                             [-1,  1,  1]])

        self.nodes = nodes
        self.gauss_points = []
        for i in [-1, 1]:
            for j in [-1, 1]:
                for k in [-1, 1]:
                    coordinates = (i/np.sqrt(3), j/np.sqrt(3), k/np.sqrt(3))
                    b_matrix = self.B(*coordinates)
                    jacobian = np.linalg.det(self.J(*coordinates))
                    self.gauss_points.append(GaussPoint(coordinates=coordinates, B=b_matrix, detJ=jacobian,
                                                        material=material()))

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
            d[0, i] = (1.+eta*self.pos[i, 1])*(1.+zeta*self.pos[i, 2])*self.pos[i, 0]/8
            d[1, i] = (1.+xi*self.pos[i, 0])*(1.+zeta*self.pos[i, 2])*self.pos[i, 1]/8
            d[2, i] = (1.+xi*self.pos[i, 0])*(1.+eta*self.pos[i, 1])*self.pos[i, 2]/8
        return d

    # noinspection PyPep8Naming
    def J(self, xi, eta, zeta):
        pos_matrix = np.zeros((8, 3))
        for i in range(3):
            pos_matrix[:, i] = [n.coordinates[i] for n in self.nodes]
        return np.dot(self.d_matrix(xi, eta, zeta), pos_matrix)

    def update(self, nodal_displacements):
        for gp in self.gauss_points:
            strain = np.dot(gp.B, nodal_displacements)
            gp.material.update(strain)

    # noinspection PyPep8Naming
    def Kt(self):




if __name__ == '__main__':
    node_list = [Node(coordinates=[0, 0, 0]),
                 Node(coordinates=[1, 0, 0]),
                 Node(coordinates=[1, 1, 0]),
                 Node(coordinates=[0, 1, 0]),
                 Node(coordinates=[0, 0, 1]),
                 Node(coordinates=[1, 0, 1]),
                 Node(coordinates=[1, 1, 1]),
                 Node(coordinates=[0, 1, 1])]

    element = Element8(nodes=node_list, ElasticMaterial)
    print element.B(0, 0, 0)
