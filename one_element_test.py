import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from fem_functions import Node, Element8
from materials.SS2506 import test_material
from transformation_material import TransformationMaterial
from nonlinear_solver import NewtonRahpson


def one_element_simulation(material_parameters, increments, exx=None, eyy=None, ezz=None, pxx=None, pyy=None, pzz=None,
                           **_):
    node_list = [Node(coordinates=[0, 0, 0]),
                 Node(coordinates=[1, 0, 0]),
                 Node(coordinates=[1, 1, 0]),
                 Node(coordinates=[0, 1, 0]),
                 Node(coordinates=[0, 0, 1]),
                 Node(coordinates=[1, 0, 1]),
                 Node(coordinates=[1, 1, 1]),
                 Node(coordinates=[0, 1, 1])]

    element = Element8(nodes=node_list, material=TransformationMaterial,
                       material_parameters={'parameters': material_parameters, 'fM': 0.8})

    solver = NewtonRahpson(node_list, [element])
    solver.set_displacement_bc(nodes=[0, 3, 4, 7], components=[0])
    solver.set_displacement_bc(nodes=[0, 1, 4, 5], components=[1])
    solver.set_displacement_bc(nodes=[0, 1, 2, 3], components=[2])

    strains = np.zeros((increments, 3))
    stresses = np.zeros((increments, 3))

    old_exx = 0
    old_eyy = 0
    old_ezz = 0

    for i in range(increments):
        if exx is not None:
            solver.set_displacement_bc(nodes=[1, 2, 5, 6], components=[0], value=exx[i, 1] - old_exx)
            old_exx = exx[i, 1]
        if eyy is not None:
            solver.set_displacement_bc(nodes=[2, 3, 6, 7], components=[1], value=eyy[i, 1] - old_eyy)
            old_eyy = eyy[i, 1]
        if ezz is not None:
            solver.set_displacement_bc(nodes=[4, 5, 6, 7], components=[2], value=ezz[i, 1] - old_ezz)
            old_ezz = ezz[i, 1]

        if pxx is not None:
            solver.set_force_bc(nodes=[1, 2, 5, 6], components=[0], value=pxx[i, 1]/4)
        if pyy is not None:
            solver.set_force_bc(nodes=[2, 3, 6, 7], components=[1], value=pyy[i, 1]/4)
        if pzz is not None:
            solver.set_force_bc(nodes=[4, 5, 6, 7], components=[2], value=pzz[i, 1]/4)

        solver.solve()
        for comp in range(3):
            strains[i, comp] = solver.u[6*3 + comp]
            stresses[i, comp] = solver.reaction_forces()[6*3 + comp]*4
    return strains, stresses


if __name__ == '__main__':
    matplotlib.style.use('classic')
    plt.rc('text', usetex=True)
    plt.rc('font', serif='Computer Modern Roman')
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                      'monospace': ['Computer Modern Typewriter']})

    time = np.linspace(0., 1., 1000)
    strain_z = np.zeros((1000, 2))
    strain_z[:, 0] = time
    strain_z[:, 1] = 0.02*np.sin(2*np.pi*time)
    e, s = one_element_simulation(1000, test_material, ezz=strain_z)

    plt.plot(e[:, 2], s[:, 2])
    plt.figure()
    plt.plot(e[:, 2], e[:, 0])
    plt.show()
