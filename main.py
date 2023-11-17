import numpy as np
import matplotlib.pyplot as plt

# Constants
HBAR = 6.582119514e-16  # eV*s
MASS_E = 5.6856300620910916e-30  # eV*s^2/nm^2

L = 2  # Barrier width
V0 = 1  # Barrier height
NB = 2  # Number of barriers
LF = 2  # Free space width
MOD = L + LF  # Periodicity of the barrier
LM = (L * NB) + ((NB - 1) * LF)  # Last barrier position
E = 0.1  # Electron energy eV

def barrier_function(x):
    if x >= LM:
        return 0
    if (x % MOD) <= (L - 1):
        return V0
    else:
        return 0

def calculate_transfer_matrix():
    k = np.emath.sqrt(2 * MASS_E * (E - V0) / HBAR ** 2)
    return np.array([[np.cos(k * L), -1j * np.sin(k * L) / k],
                     [-1j * k * np.sin(k * L), np.cos(k * L)]])

def overall_transfer_matrix():
    overall_matrix = np.identity(2, dtype=complex)
    transfer_matrix = calculate_transfer_matrix()
    for _ in range(NB):
        overall_matrix = np.dot(transfer_matrix, overall_matrix)
    return overall_matrix

def reflection_probability(x):
    overall_matrix = overall_transfer_matrix()
    reflection_prob = np.abs(overall_matrix[1, 0] / overall_matrix[0, 0])**2
    return reflection_prob * barrier_function(x)

def main():
    # create x values range from 0 to 10 step to 1
    x_values = np.arange(0, 10, 0.1)
        
    plt.plot(x_values, [barrier_function(x) for x in x_values])
    plt.show()
    

if __name__ == '__main__':
    main()
