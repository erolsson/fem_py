import numpy as np

from fem_functions import Node, Element8
# from elastic_material import ElasticMaterial
from ss2506material import SS2506
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
    SS2506(E=200e3, v=0.3, a=[0.056, 0.028, 0.], dV=0.037, R_params=[0.02, 0.02], sy0_au=400, fM=0.8)
    element = Element8(nodes=node_list, material=SS2506,
                       material_parameters={'E': 200e3, 'v': 0.3,
                                            'a': [0.056, 0.028, 0.],
                                            'dV': 0.037, 'R_params': [0.02, 0.02], 'sy0_au': 400,
                                            'fM': 0.8})

    solver = NewtonRahpson(node_list, [element])
    solver.set_displacement_bc(nodes=[0, 3, 4, 7], components=[0])
    solver.set_displacement_bc(nodes=[0, 1, 4, 5], components=[1])
    solver.set_displacement_bc(nodes=[0, 1, 2, 3], components=[2])
    s_max = 1000
    nodal_forces = np.linspace(0, s_max/4, 1000)
    for f in nodal_forces:
        solver.f[14] = f
        solver.f[17] = f
        solver.f[20] = f
        solver.f[23] = f
        solver.solve()
        print solver.u[20]
