from fem_functions import Node, Element8
from elastic_material import ElasticMaterial
from nonlinear_solver import NewtonRahpson

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

    solver = NewtonRahpson(node_list, [element])
    print solver.solve()
