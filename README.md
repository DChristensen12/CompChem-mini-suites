# Computational Chemistry Mini‑Projects Suite

A curated collection of four mini‑projects that showcase core techniques in modern computational chemistry and molecular simulation.  The suite is designed for instructional use and can be tackled sequentially or à la carte.

> **Goal**  Show how quantum mechanics, statistical mechanics, and molecular dynamics work together to predict molecular properties across temperature and phase.

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Project 1 – Quantum & Classical Harmonic Oscillators](#project-1)
3. [Project 2 – Quantum Chemistry of the Hydrogen Molecule](#project-2)
4. [Project 3 – Path‑Integral Monte Carlo for H₂](#project-3)
5. [Project 4 – Phase Diagram of Lennard‑Jones Fluids](#project-4)
6. [Reproducibility & Testing](#reproducibility--testing)
7. [Citing & Acknowledgements](#citing--acknowledgements)

---

## Getting Started

```bash
# Clone the repository
$ git clone https://github.com/<your‑user>/compchem-mini-suite.git
$ cd compchem-mini-suite

# Create & activate a fresh environment (conda or venv)
$ conda env create -f environment.yml
$ conda activate compchem-mini

# Run notebooks/examples
$ jupyter lab
```

**Dependencies**

* Python ≥ 3.11
* `numpy`, `scipy`, `sympy`, `matplotlib`, `pandas`
* `numba` (JIT acceleration)
* `ase` (optional—atomic visualization)
* `pytest` (unit tests)

The *environment.yml* file pins exact versions to guarantee identical results on any machine.

---

## Project 1 – Quantum & Classical Harmonic Oscillators  <a id="project-1"></a>

> *Comparative Study of Quantum and Classical Harmonic Oscillators*

Directory: **`project1/`**

| Part | Notebook / Script            | What You Do                                                                                                                                                      | Key Outputs                                                        |              |                                    |
| ---- | ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ------------ | ---------------------------------- |
| 1‑1  | `1‑1_dimensionless.ipynb`    | Derive a dimension‑less Hamiltonian by scaling $x \to \xi = x/x_0$.                                                                                              | Closed‑form expression for $ \hat H = \tfrac12 (p^2 + \xi^2)$.     |              |                                    |
| 1‑2  | `qm_harmonic.py`             | Build a finite basis of HO eigenfunctions, diagonalise with `scipy.linalg.eigh`, compute first three energies, and evaluate Bose–Einstein averages.              | Table of eigenvalues; plots of (                                   | \psi\_n(\xi) | ^2); $\langle E(T)\rangle$ vs *T*. |
| 1‑3  | `md_ho.py`                   | Classical MD via velocity Verlet + Langevin thermostat.  Uses atomic units ($m = \omega = 1$).  Block averaging gives energies, heat capacities, and error bars. | `results/ho_md_*` CSV files; PDF plots of energy & position PDFs.  |              |                                    |
| 1‑4  | `compare_qm_classical.ipynb` | Cross‑compare $\langle E\rangle$ and $C_V$ from Parts 1‑2 & 1‑3; comment on equipartition breakdown.                                                             | Combined plot highlighting quantum/classical disparity at low *T*. |              |                                    |

**Implementation Highlights**

* **Symbolic Derivation** (Part 1‑1) uses SymPy to keep algebra transparent.
* **Basis Truncation** (Part 1‑2) exposes convergence; you can tune `N_basis` in the config block.
* **Thermostat Choice** (Part 1‑3) defaults to Langevin (γ = 0.5), but any Nose–Hoover chain can be plugged in.
* **Error Bars** handled via block‑averaging; block length chosen from the autocorrelation time estimated on the fly.

---

## Project 2 – Quantum Chemistry of the Hydrogen Molecule  <a id="project-2"></a>

> *Quantum Mechanical Investigation of the H₂ Covalent Bond*

Directory: **`project2/`**

| Part | Notebook / Script     | Core Tasks                                                                                                                                                                      | Deliverables                                      |
| ---- | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------- |
| 2‑1  | `hf_solver.py`        | Build a *restricted* Hartree–Fock code from scratch.  Integral engine uses Obara–Saika recurrence; STO‑3G & STO‑4G parameters fetched at runtime from Basis‑Set‑Exchange API.   | Converged SCF energies at multiple bond lengths.  |
| 2‑2  | `hf_optim.ipynb`      | Profile the naive integrals (`cProfile`), then accelerate with `numba` & BLAS.                                                                                                  | Bar chart of speedups; timing table (old vs new). |
| 2‑3  | `pes_scan.py`         | Sweep 1.0–3.0 Å (≤ 0.05 Å spacing), record total energy.  `plot_pes.ipynb` plots the curve and extracts $R_e$ & dissociation energy with `scipy.optimize.curve_fit`.            | PNG of PES; printed $R_e$ and $D_e$.              |
| 2‑4  | `potential_fit.ipynb` | Fit Morse & Lennard‑Jones to the ab‑initio energies.  Uses Levenberg–Marquardt and reports $\chi^2$.  Discuss discrepancies and propose multi‑body or perturbative corrections. | Overlaid PES plot + fit parameters table.         |

**Implementation Highlights**

* **Integral Caching** To avoid $O(N^4)$ recomputation, integrals are stored in HDF5.
* **DIIS Acceleration** ensures robust SCF convergence (<15 cycles typical).
* **Vectorised Scan** Bond distances are distributed to CPUs via `multiprocessing.Pool`.

---

## Project 3 – Path‑Integral Monte Carlo for H₂  <a id="project-3"></a>

> *Finite‑Temperature Quantum Statistics with Ring‑Polymer Mapping*

Directory: **`project3/`**

| Part | File                             | Focus                                                                                                                                                                  | Output                                               |
| ---- | -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| 3‑1  | `virial_estimator.md`            | Paper‑style derivation of the Virial energy estimator for a Morse potential.                                                                                           | Markdown document with equation derivations.         |
| 3‑2  | `pimc_core.py`                   | PIMC engine: `RingPolymer` class, periodic box (L = 10 a.u.), Morse potential, single‑bead Metropolis moves.                                                           | Trajectory `.npz`; log of acceptance ratios.         |
| 3‑3  | `pimc_core.py` (module extended) | Collective whole‑polymer translations; acceptance calculated from action difference.                                                                                   | Updated acceptance stats; faster energy convergence. |
| 3‑4  | `analysis.ipynb`                 | Runs at 1000 K & 2000 K; chooses step size to target 30 % acceptance.  Uses blocking to compute $\langle E \rangle\pm\sigma$.  Computes radial PDF P(r) via histogram. | Energy vs T table; P(r) plots for both temps.        |

**Implementation Highlights**

* **Autotuned Step Size** First 5 × 10³ moves adapt Δ to hit desired acceptance, then freeze.
* **Periodic Images** Minimum‑image convention wrapped in `vectorized_distance()` for speed.
* **Blocking & Jackknife** Used for error bars on mean energy and histogram counts.

---

## Project 4 – Phase Diagram of Lennard‑Jones Fluids  <a id="project-4"></a>

> *Mapping Phase Behaviour and Transport in a Simple Fluid Model*

Directory: **`project4/`**

| Part | Script                                 | Objective                                                                                                           | Artifacts                                   |
| ---- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| 4‑1  | `lj_system.py`                         | Base `LJSystem` class with cell‑list neighbour search (cutoff = 3σ).                                                | Module ready for import by later parts.     |
| 4‑2  | `integrators.py`                       | Velocity‑Verlet in reduced units (Δt⋆ = 1e‑3).                                                                      | Verified energy drift < 10⁻⁴ ε per 10 τ.    |
| 4‑3  | `thermostat.py`, `equilibration.ipynb` | Nosé–Hoover chain thermostat; run at ρσ³ = 0.6, T⋆ = 1.5 for 10 τ.                                                  | Plot of E\_kin, E\_pot, E\_tot stabilising. |
| 4‑4  | `phase_scan.py`                        | Grid of (T⋆, ρσ³) points (table in spec).  Runs 20 τ starting from FCC & liquid.  Computes RDF `g(r)` for phase ID. | `phase_map.csv`; RDF plots.                 |
| 4‑5  | `analysis_phase_boundary.md`           | Discuss hysteresis, propose interface‑pinning & thermodynamic integration strategy.                                 | Markdown report.                            |
| 4‑6  | `diffusion.py`                         | Computes D from MSD (Einstein) and VACF (Green–Kubo) at ρσ³ = 0.1, 0.6.                                             | Table & plot of D vs ρ.                     |

**Implementation Highlights**

* **Cell Lists** `numba`‑accelerated double loop gives $O(N)$ force evaluation.
* **Structure Factor** `scipy.fft` option to characterise long‑range order (liquid vs solid).
* **Diffusion Analysis** Linear regression on MSD after ballistic regime; VACF integrated with trapezoidal rule.

---

## Reproducibility & Testing  <a id="reproducibility--testing"></a>

```bash
# Run end‑to‑end test suite
$ pytest -q
```

Each project folder contains *unit* tests (e.g., force consistency, energy conservation) and *regression* tests comparing against stored reference data.

CI is configured via GitHub Actions (see `.github/workflows/ci.yml`).  Every push runs the notebooks in headless mode and uploads artifacts.

---

## Citing & Acknowledgements  <a id="citing--acknowledgements"></a>

Please cite this repository as:

```
Christensen, Daniel; *Computational Chemistry Mini‑Projects Suite*, 2025
```

Special thanks to *Basis Set Exchange* for open access to basis‑set parameters and to the open‑source community for `numpy`, `scipy`, and `numba`.
