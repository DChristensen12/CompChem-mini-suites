# -*- coding: utf-8 -*-

## Part 1-2: Quantum Solutions and Bose-Einstein Statistics of Harmonic Oscillator

# -----  Parameters ------ #
N = 100              # number of grid points
x_max = 5.0          # domain for x*
x = np.linspace(-x_max, x_max, N)
dx = x[1] - x[0]

# ===== Kinetic energy operator (finite difference method) ===== #
T = -np.diag(np.ones(N-1), -1) + 2*np.diag(np.ones(N), 0) - np.diag(np.ones(N-1), 1)
T *= 1 / dx**2

# ===== Potential energy operator ===== #
V = np.diag(x**2)

# ===== Hamiltonian in reduced units (scaled by ħω/2) ===== #
H = T + V

# ----- Diagonalize the Hamiltonian ----- #
e_vals, e_vecs = eigh(H)
energies = e_vals[:3]           # Three lowest eigenvalues
wavefuncs = e_vecs[:, :3]       # Corresponding eigenfunctions

# ===== Normalize wavefunctions and compute probability densities ===== #
prob_densities = []
for i in range(3):
    psi = wavefuncs[:, i]
    norm = np.trapz(np.abs(psi)**2, x)
    prob_densities.append(np.abs(psi)**2 / norm)

# ===== Plots probability densities ===== #
plt.figure(figsize=(8, 5))
for n, prob in enumerate(prob_densities):
    plt.plot(x, prob, label=f"$n={n}$, $E^*={energies[n]:.3f}$")

plt.title("Lowest Three Quantum Harmonic Oscillator States\n(Dimensionless Units)")
plt.xlabel("Dimensionless position $x^*$")
plt.ylabel("Probability density")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#  ----- Bose-Einstein average energy function  ----- #
def bose_avg_energy(energies, T_star):
    energies = np.array(energies)
    weights = np.exp(-energies / T_star)
    Z = np.sum(weights)
    avg_E = np.sum(energies * weights) / Z
    return avg_E

# ===== Computing average energy over a range of dimensionless temperatures ===== #
T_stars = np.linspace(0.1, 5.0, 100)
avg_Es = [bose_avg_energy(energies, T) for T in T_stars]

# ===== Plotting average energy vs. temperature ===== #
plt.figure(figsize=(7, 4))
plt.plot(T_stars, avg_Es, label="Quantum Avg. Energy", linewidth=2)
plt.xlabel("Dimensionless Temperature $T^* = k_B T / (\hbar \omega / 2)$")
plt.ylabel("Average Energy $\langle E^* \\rangle$")
plt.title("Bose-Einstein Average Energy of Harmonic Oscillator")
plt.grid(True)
plt.tight_layout()
plt.show()

# ===== Outputting eigenvalues for reference ===== #
for n, E in enumerate(energies):
    print(f"n = {n}, E* = {E:.6f}")
