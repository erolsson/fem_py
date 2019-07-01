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

    element = Element8(nodes=node_list, material=ElasticMaterial, material_parameters={'E': 200.e3, 'v': 0.3})

    solver = NewtonRahpson(node_list, [element])
    solver.set_displacement_bc(nodes=[0, 3, 4, 7], components=[0])
    solver.set_displacement_bc(nodes=[0, 1, 4, 5], components=[1])
    solver.set_displacement_bc(nodes=[0, 1, 2, 3], components=[2])

    solver.f[14] = 1.
    solver.f[17] = 1.
    solver.f[20] = 1.
    solver.f[23] = 1.
    solver.solve()
    print solver.u
