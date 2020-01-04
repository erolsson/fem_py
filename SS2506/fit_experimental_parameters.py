from __future__ import print_function
import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt

from neu_sehitoglu.fit_uniaxial_behaviour import kinematic_hardening, kinematic_params_residual

# Plotting curves with no transformation
core_data = np.genfromtxt('experimental_data/tension_data_core_BA.dat')[:-10, :]
case_data = -np.genfromtxt('experimental_data/compression_data_case_EN.dat', delimiter=',')
num_back_stresses = 3
kinematic_parameters = np.zeros((2, num_back_stresses*2))
hardness = np.array([450., 770.])
starting_params = [15432, 5., 281622, 236, 470894, 2301]
yield_stresses = np.zeros(2)
E = 205E3
yield_stress_offset = 1e-4
for i, data_set in enumerate([core_data, case_data]):
    strain = data_set[:, 0]
    stress = data_set[:, 1]
    plt.figure(0)
    plt.plot(strain, stress, lw=3)
    epl = strain - stress/E
    sy0 = np.interp(yield_stress_offset, epl, stress)
    print(sy0)

    sy = stress[epl > yield_stress_offset]
    yield_stresses[i] = sy0
    epl = epl[epl > yield_stress_offset]
    plt.figure(1)
    plt.plot(epl, sy)
    params = fmin(kinematic_params_residual, starting_params, args=(epl, sy - sy0), maxfun=1e6, maxiter=1e6)
    kinematic_parameters[i, :] = abs(params)
    alpha = kinematic_hardening(epl, params)
    plt.plot(epl, sy0 + alpha)
    print(sy0, params)

epl = np.linspace(0, 10e-2, 1000)

for hv in np.linspace(450, 800, 50):
    kin_parameters = np.zeros(2*num_back_stresses)
    for i in range(2*num_back_stresses):
        kin_parameters[i] = np.interp(hv, hardness, kinematic_parameters[:, i])
    sy = np.interp(hv, hardness, yield_stresses)
    stress = 0*epl + sy + kinematic_hardening(epl, kin_parameters)
    strain = epl + stress/E
    plt.figure(0)
    plt.plot(strain, stress)

syA = np.array([0, 420.])
syM = np.array([yield_stresses[0], (yield_stresses[1] - 0.2*syA[1])/0.8])
syA[0] = syA[1]*syM[0]/syM[1]

print(" --- Initial yield stresses ---")
print("Austenite:")
syA_b = np.diff(syA)/(770-450)
syA_a = syA[0] - syA_b*450
print('\t', syA_a, syA_b)
print("Martensite:")
print(syM)
syM_b = np.diff(syM)/(770-450)
syM_a = syM[0] - syM_b*450
print('\t', syM_a, syM_b)
print(" --- Hardening ---")
C = kinematic_parameters[:, 0::2]
g = kinematic_parameters[:, 1::2]
C_b = (C[1, :] - C[0, :])/(770-450)
g_b = (g[1, :] - g[0, :])/(770-450)
C_a = C[0, :] - C_b*450
g_a = g[0, :] - g_b*450
print("Cm")
print("\t", C_a, C_b)
print("gm")
print("\t", g_a, g_b)

plt.figure(2)
strain_data = np.genfromtxt('experimental_data/transversal_strain_tension_EN.dat', delimiter=',')
plt.plot(strain_data[:, 0], strain_data[:, 1])

plt.figure(3)
case_data = -np.genfromtxt('experimental_data/compression_data_case_EN.dat', delimiter=',')
plt.plot(case_data[:, 0], case_data[:, 1])
case_data = np.genfromtxt('experimental_data/tension_data_case_EN.dat', delimiter=',')
plt.plot(case_data[:, 0], case_data[:, 1])

plt.show()
