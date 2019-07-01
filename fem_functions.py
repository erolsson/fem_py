import numpy as np


class Node:
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)


class Element8:
    def __init__(self, nodes):
        self.pos = np.array([[-1, -1, -1],
                             [1,  -1, -1],
                             [1,   1, -1],
                             [-1,  1, -1],
                             [-1, -1,  1],
                             [1,  -1,  1],
                             [1,   1,  1],
                             [-1,  1,  1]])

        self.nodes = nodes

    def B(self, xi, eta, zeta):
        B_matrix = np.zeros((6, 24))

    def d_matrix(self, xi, eta, zeta):
        d = np.zeros((3, 8))
        for i in range(8):
            d[0, i] = (1.+eta*self.pos[i, 1])*(1.+zeta*self.pos[i, 2])*self.pos[i, 0]/8
            d[1, i] = (1.+xi*self.pos[i, 0])*(1.+zeta*self.pos[i, 2])*self.pos[i, 1]/8
            d[2, i] = (1.+xi*self.pos[i, 0])*(1.+eta*self.pos[i, 1])*self.pos[i, 2]/8
        return d

    def J(self, xi, eta, zeta):
        pos_matrix = np.zeros((8, 3))
        for i in range(3):
            pos_matrix[:, i] = [n.coordinates[i] for n in self.nodes]
        return np.dot(self.d_matrix(xi, eta, zeta), pos_matrix)


if __name__ == '__main__':
    node_list = [Node(coordinates=[0, 0, 0]),
                 Node(coordinates=[1, 0, 0]),
                 Node(coordinates=[1, 1, 0]),
                 Node(coordinates=[0, 1, 0]),
                 Node(coordinates=[0, 0, 1]),
                 Node(coordinates=[1, 0, 1]),
                 Node(coordinates=[1, 1, 1]),
                 Node(coordinates=[0, 1, 1])]

    element = Element8(nodes=node_list)
    print element.J(0, 0, 0)
