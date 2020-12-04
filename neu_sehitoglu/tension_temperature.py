import os

import numpy as np

import matplotlib.pyplot as plt

from fem_py.materials.transformation_materials import neu_sehitoglu
from fem_py.abaqus_material_test.material_test import one_element_abaqus
from fem_py.abaqus_material_test.one_element_input_file import BC


def run_sim_from_experiment(name, temperature, simulation_dir, umat_file, boundary_conditions):
    e, s, _, _ = one_element_abaqus(simulation_dir, material=neu_sehitoglu,
                                    boundary_conditions=boundary_conditions,
                                    simulation_name=name,
                                    temperature=np.array([[0, temperature], [1, temperature]]),
                                    user_subroutine=umat_file,
                                    max_increment=0.01)
    return e, s


def main():
    umat_file = os.path.expanduser('~/python_projects/fem_py/transformation_subroutines/transformation_subroutine.o')
    simulation_dir = os.path.expanduser('~/python_projects/fem_py/abaqus_material_test/neu_sehitoglu')

    neu_sehitoglu.sde = 0
    temp = 22
    for direction, color in zip(['compression', 'tension'], ['r', 'b']):
        data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig2_' + direction),
                             delimiter=',')
        plt.figure(1)
        plt.plot(data[:, 0], data[:, 1], color + '*')
        sim_name = 'stress_22_' + direction
        if direction == 'tension':
            boundary_conditions = [BC(amplitude=[(0., 0.), (1., data[-1, 1])], direction='z', mode='stress')]
        else:
            boundary_conditions = [BC(amplitude=[(0., 0.), (1., -data[-1, 1])], direction='z', mode='stress')]
        strain, stress = run_sim_from_experiment(sim_name, temp, simulation_dir, umat_file, boundary_conditions)

        plt.plot(np.abs(strain[:, 2]), np.abs(stress[:, 2]), color, lw=2)

        plt.figure(2)
        dv = np.sum(strain[:, 0:3], 1) - stress[:, 2]/neu_sehitoglu.E*(1-2*neu_sehitoglu.v)
        data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig3_' + direction),
                             delimiter=',')
        plt.plot(data[:, 0], data[:, 1], color + '*')
        plt.plot(np.abs(strain[:, 2]), dv, color, lw=2)

    plt.figure(4)
    data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig8'), delimiter=',')
    plt.plot(data[:, 0], data[:, 1], '-b*')

    data = np.genfromtxt(os.path.expanduser('~/phase_transformations/neu_sehitoglu/fig7'), delimiter=',')
    plt.figure(3)
    plt.plot(2*data[:, 0], data[:, 1], '-b*')
    print(data[0, 1]/data[0, 0])
    plt.plot([0, 0.02], [0, 4E3/2/(1+0.3)])
    plt.figure(1)
    plt.plot(np.sqrt(3)*data[:, 0], data[:, 1]*np.sqrt(3))
    plt.plot([0, 0.01], [0, 2E3])

    sim_name = 'stress_22_torsion'
    boundary_conditions = [BC(amplitude=[(0., 0.), (1., data[-1, 1])], direction='y', mode='stress'),
                           BC(amplitude=[(0., 0.), (1., -data[-1, 1])], direction='z', mode='stress')]
    strain, stress = run_sim_from_experiment(sim_name, temp, simulation_dir, umat_file, boundary_conditions)
    gamma = strain[:, 1] - strain[:, 2]
    tau = (stress[:, 1] - stress[:, 2])/2
    plt.figure(3)
    plt.plot(gamma, tau)

    plt.figure(4)
    plt.plot(tau, np.sum(strain[:, 0:3], axis=1))

    plt.show()


if __name__ == '__main__':
    main()
