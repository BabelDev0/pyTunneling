import numpy as np
import matplotlib.pyplot as plt

hb = 1  # h bar in atomic units
eVtoHa = 0.03674  # eV to hartree
AtoB = 1.88973  # Angstrom to Bohr
VeV = 0.100  # Barrier height in eV
V = VeV * eVtoHa  # Barrier height in Hartree
m0 = 1  # electron mass in a.u.
m1 = m0 * 0.067  # effective mass in the well
m2 = m0 * 0.067  # effective mass in the barrier
L1A = 50  # thickness of the first barrier in Ang
L2A = 50  # thickness of the well in Ang
L3A = 50  # thickness of the second barrier in Ang
I1 = 0
I2 = L1A * AtoB
I3 = (L1A + L2A) * AtoB
I4 = (L1A + L2A + L3A) * AtoB

E = 0
num_points = 10000
x = np.zeros(num_points)
T = np.zeros(num_points)

for n in range(num_points):
    E += V / num_points
    k1 = np.sqrt(2 * m1 * E) / hb
    k2 = np.sqrt(2 * m2 * (V - E)) / hb

    M1 = np.array([[1, 1], [1j * k1, -1j * k1]])
    M2 = np.array([[1, 1], [k2, -k2]])
    M3 = np.array([[np.exp(k2 * I2), np.exp(-k2 * I2)], [k2 * np.exp(k2 * I2), -k2 * np.exp(-k2 * I2)]])
    M4 = np.array([[np.exp(1j * k1 * I2), np.exp(-1j * k1 * I2)],
                   [1j * k1 * np.exp(1j * k1 * I2), -1j * k1 * np.exp(-1j * k1 * I2)]])
    M5 = np.array([[np.exp(1j * k1 * I3), np.exp(-1j * k1 * I3)],
                   [1j * k1 * np.exp(1j * k1 * I3), -1j * k1 * np.exp(-1j * k1 * I3)]])
    M6 = np.array([[np.exp(k2 * I3), np.exp(-k2 * I3)], [k2 * np.exp(k2 * I3), -k2 * np.exp(-k2 * I3)]])
    M7 = np.array([[np.exp(k2 * I4), np.exp(-k2 * I4)], [k2 * np.exp(k2 * I4), -k2 * np.exp(-k2 * I4)]])
    M8 = np.array([[np.exp(1j * k1 * I4), np.exp(-1j * k1 * I4)],
                   [1j * k1 * np.exp(1j * k1 * I4), -1j * k1 * np.exp(-1j * k1 * I4)]])

    M = np.linalg.inv(M1) @ M2 @ np.linalg.inv(M3) @ M4 @ np.linalg.inv(M5) @ M6 @ np.linalg.inv(M7) @ M8
    x[n] = E / eVtoHa
    T[n] = 1 / (M[0, 0] * np.conj(M[0, 0]))

plt.plot(x, T)
plt.xlabel('Energy (Hartree)')
plt.ylabel('Transmission Probability')
plt.title('Transmission Probability through a Double Barrier System')
plt.show()
