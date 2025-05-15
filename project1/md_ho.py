# -*- coding: utf-8 -*-

## Part 1-3: Molecular Dynamics Simulation and Thermodynamic Properties of the Classical Oscillator
# ----- Atomic units (a.u.) ----- #
mass = 1.0      # mass in atomic units
omega = 1.0     # angular frequency in atomic units
k_B = 1.0       # Boltzmann constant in atomic units

# ----- Velocity Verlet with Langevin thermostat ----- #
def velocity_verlet_langevin(T, dt, steps, gamma=1.0):
    """
    Integrate a 1D harmonic oscillator with a Langevin thermostat.
    Returns arrays of total energies and positions.
    """
    x = 0.0
    v = 0.0
    x_vals = []
    E_vals = []

    # Noise amplitude for stochastic force: sqrt(2*gamma*k_B*T*dt)/mass
    sigma = np.sqrt(2 * gamma * k_B * T * dt) / mass

    for _ in range(steps):
        # 1) Deterministic half-step velocity update
        a = -omega**2 * x
        v += 0.5 * a * dt

        # 2) Position update
        x += v * dt

        # 3) Complete deterministic velocity update
        a_new = -omega**2 * x
        v += 0.5 * a_new * dt

        # 4) Langevin thermostat: damping + random kick
        v = v * (1 - gamma * dt) + sigma * np.random.normal()

        # 5) Record position and total energy
        x_vals.append(x)
        E_kin = 0.5 * mass * v**2
        E_pot = 0.5 * mass * omega**2 * x**2
        E_vals.append(E_kin + E_pot)

    return np.array(E_vals), np.array(x_vals)

# ----- Estimate average energy and heat capacity with error bars ----- #
def estimate_Cv(E_vals, T):
    """
    Compute:
      • mean total energy
      • heat capacity via fluctuations
      • standard error of the mean energy
    """
    E = np.array(E_vals)
    mean_E = np.mean(E)
    var_E = np.var(E)
    Cv = var_E / (k_B * T**2)
    std_err = np.std(E) / np.sqrt(len(E))
    return Cv, mean_E, std_err

# ----- Simulation parameters ----- #
temps = [0.5, 1.0, 2.0, 4.0]  # temperatures to sample
dt = 0.01                    # time step
steps = 100_000              # number of MD steps per temperature

# ===== Run simulations across temperatures ===== #
results = []
for T in temps:
    E_vals, x_vals = velocity_verlet_langevin(T, dt, steps)
    Cv, avg_E, err_E = estimate_Cv(E_vals, T)
    results.append((T, avg_E, err_E, Cv))

# ===== Display results as DataFrame ===== #
df = pd.DataFrame(results,
                  columns=["T", "Average E", "Error E", "Heat Capacity Cv"])
print(df)

# ===== Plot histogram of x for T = 1.0 and compare with canonical distribution ===== #
T_plot = 1.0
E_vals_plot, x_vals_plot = velocity_verlet_langevin(T_plot, dt, steps)

# Compute normalized histogram (PDF)
hist, bins = np.histogram(x_vals_plot, bins=100, density=True)
bin_centers = 0.5 * (bins[1:] + bins[:-1])

# ----- Canonical (Boltzmann) distribution for harmonic oscillator ----- #
canonical = np.exp(-0.5 * mass * omega**2 * bin_centers**2 / (k_B * T_plot))
canonical /= np.trapz(canonical, bin_centers)

# ----- Error bars from Poisson statistics on histogram bin counts ----- #
bin_counts, _ = np.histogram(x_vals_plot, bins=100)
pdf_errors = np.sqrt(bin_counts) / (len(x_vals_plot) * (bins[1] - bins[0]))

# ===== Plotting ===== #
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers,
             hist,
             yerr=pdf_errors,
             fmt='o',
             alpha=0.5,
             label='Simulated PDF')
plt.plot(bin_centers,
         canonical,
         label='Canonical Distribution',
         linewidth=2)
plt.title("Position Probability Density at $T = 1.0$")
plt.xlabel("Position $x$")
plt.ylabel("Probability Density")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ===== Fit Gaussian Distribution ===== #
# Compute sample mean and standard deviation of the trajectory
mu    = np.mean(x_vals_plot)
sigma = np.std(x_vals_plot)

# Build Gaussian PDF at the bin centers
gauss_pdf = (1.0/(sigma * np.sqrt(2*np.pi))) * \
             np.exp(-0.5 * ((bin_centers - mu)/sigma)**2)

# ===== Plot: Histogram + Gaussian Fit + Canonical Curve ===== #
plt.figure(figsize=(8,5))

# 1) Underlying histogram of the data (for context)
plt.hist(x_vals_plot,
         bins=100,
         density=True,
         alpha=0.2,
         color='tab:blue',
         label='Histogram of $x$')

# 2) Smooth Gaussian fit
plt.plot(bin_centers,
         gauss_pdf,
         linestyle='-',
         color='tab:green',
         linewidth=2,
         label=f'Gaussian Fit\n($\\mu={mu:.2f},\\,\\sigma={sigma:.2f}$)')

# 3) Canonical (Boltzmann) distribution
plt.plot(bin_centers,
         canonical,
         linestyle='--',
         color='tab:orange',
         linewidth=2,
         label='Canonical Distribution')

plt.title("Gaussian Fit vs. Canonical PDF at $T=1.0$")
plt.xlabel("Position $x$")
plt.ylabel("Probability Density")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
