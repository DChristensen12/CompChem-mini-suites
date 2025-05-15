# Computational Chemistry Mini‑Projects Suite

A teaching‑oriented collection of **four stand‑alone mini‑projects**—plus a companion *all‑in‑one* Jupyter notebook—for exploring core ideas in quantum chemistry, statistical mechanics, and molecular dynamics.

```
repo‑root/
│
├─ notebook/                  # one‑stop interactive walkthrough
│   └─ CompChem-mini-suites_DC.ipynb
│
├─ project1/                  # Harmonic oscillator code & docs
├─ project2/                  # Hydrogen molecule HF study
├─ project3/                  # Path‑integral Monte Carlo
├─ project4/                  # Lennard‑Jones phase diagram
│
└─ environment.yml            # reproducible conda env
```

Choose your adventure:

1. **Run everything interactively** in the notebook folder, or
2. **Dive into each project directory** for modular scripts/notebooks, unit tests, and detailed READMEs.

---

## Quick Start

```bash
# Clone the repo and install dependencies
$ git clone https://github.com/<your-user>/compchem-mini-suite.git
$ cd compchem-mini-suite
$ conda env create -f environment.yml
$ conda activate compchem-mini

# Option A – launch the master notebook
$ jupyter lab notebook/CompChem-mini-suites_DC.ipynb

# Option B – run a specific project
$ cd project2
$ python hf_solver.py --help  # or open the local notebook
```

### Software Stack

This suite relies on the **standard scientific‑Python toolkit** (NumPy, SciPy, Matplotlib, etc.) plus a handful of extras for performance and visualisation. *All* exact requirements—including pin‑exact versions—are captured in **`environment.yml`**, so you don’t need to hunt each one down manually.  If you add new packages, just update that file and re‑export.

---

## Table of Contents

*Jump within this README or open the corresponding folder/section.*

| Project                                                | In‑Doc Anchor      | Folder Link              |
| ------------------------------------------------------ | ------------------ | ------------------------ |
| Project 1 – Quantum & Classical Harmonic Oscillators   | [jump](#project-1) | [`project1/`](project1/) |
| Project 2 – Quantum Chemistry of the Hydrogen Molecule | [jump](#project-2) | [`project2/`](project2/) |
| Project 3 – Path‑Integral Monte Carlo for H₂           | [jump](#project-3) | [`project3/`](project3/) |
| Project 4 – Phase Diagram of Lennard‑Jones Fluids      | [jump](#project-4) | [`project4/`](project4/) |

---

## Project 1 – Quantum & Classical Harmonic Oscillators  <a id="project-1"></a>

**Folder:** `project1/` · **Notebook section:** *Project 1* (within the master notebook)

| Part | File                         | Tasks                                                                                                                                                                                                                                                                                                                  | Key Outputs                                                                         |
| ---- | ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| 1‑1  | `1-1_dimensionless.ipynb`    | Derive a dimensionless Hamiltonian $\hat H = \tfrac12(p^2 + \xi^2)$.                                                                                                                                                                                                                                                   | Markdown derivation, symbolic check.                                                |
| 1‑2  | `qm_harmonic.py`             | Solve the time‑independent Schrödinger equation for the dimensionless Hamiltonian using a finite basis; **report the three lowest eigenvalues**, plot their probability densities $\lvert\psi_n(\xi)\rvert^2$, and, from those states, compute the Bose–Einstein average energy $⟨E(T)⟩$ over a range of temperatures. | Table/CSV of eigenvalues; PNG/PDF wave‑function plots; energy‑vs‑temperature curve. |
| 1‑3  | `md_ho.py`                   | Velocity‑Verlet MD (a.u., γ = 0.5); Langevin thermostat; block averages & error bars.                                                                                                                                                                                                                                  | Energies vs T; heat‑capacity plot; position PDF.                                    |
| 1‑4  | `compare_qm_classical.ipynb` | Overlay quantum vs classical $⟨E$ and $C_V$.                                                                                                                                                                                                                                                                           | Combined figure with commentary.                                                    |

*Highlights:* SymPy for algebra, `numba` for MD inner loop, on‑the‑fly autocorrelation to choose block length.

---

## Project 2 – Quantum Chemistry of the Hydrogen Molecule  <a id="project-2"></a>

**Folder:** `project2/` · **Notebook section:** *Project 2*

| Part | File                  | Description                                                                               |
| ---- | --------------------- | ----------------------------------------------------------------------------------------- |
| 2‑1  | `hf_solver.py`        | Restricted Hartree–Fock with STO‑3G & STO‑4G; integral engine uses Obara–Saika; DIIS SCF. |
| 2‑2  | `hf_optim.ipynb`      | Profile naïve integrals, accelerate with `numba` & BLAS; timing comparison.               |
| 2‑3  | `pes_scan.py`         | Sweep H–H distance, plot PES, extract $R_e$ & $D_e$.                                      |
| 2‑4  | `potential_fit.ipynb` | Fit Morse & Lennard‑Jones potentials; discuss discrepancies & improvement ideas.          |

---

## Project 3 – Path‑Integral Monte Carlo for H₂  <a id="project-3"></a>

**Folder:** `project3/` · **Notebook section:** *Project 3*

| Part | File                  | Focus                                                             |
| ---- | --------------------- | ----------------------------------------------------------------- |
| 3‑1  | `virial_estimator.md` | Derivation of Virial energy estimator for Morse potential.        |
| 3‑2  | `pimc_core.py`        | Ring‑polymer PIMC; single‑bead moves; periodic box (L = 10 a.u.). |
| 3‑3  | *same*                | Add whole‑polymer translation moves; improved convergence.        |
| 3‑4  | `analysis.ipynb`      | Runs at 1000 K & 2000 K; energy & radial PDF with error bars.     |

---

## Project 4 – Phase Diagram of Lennard‑Jones Fluids  <a id="project-4"></a>

**Folder:** `project4/` · **Notebook section:** *Project 4*

| Part | File                                   | Objective                                                     |
| ---- | -------------------------------------- | ------------------------------------------------------------- |
| 4‑1  | `lj_system.py`                         | Cell‑list LJ potential (cutoff 3σ); periodic images.          |
| 4‑2  | `integrators.py`                       | Velocity‑Verlet (Δt⋆ = 1 × 10⁻³) in reduced units.            |
| 4‑3  | `thermostat.py`, `equilibration.ipynb` | Nosé–Hoover chain; equilibrate at ρσ³ = 0.6, T⋆ = 1.5.        |
| 4‑4  | `phase_scan.py`                        | Grid of (T⋆, ρσ³) points; identify phase via RDF.             |
| 4‑5  | `analysis_phase_boundary.md`           | Discuss hysteresis & enhanced sampling for solid–liquid line. |
| 4‑6  | `diffusion.py`                         | Diffusion coefficients from MSD & VACF at two densities.      |

---

## Reproducibility & Testing

Each directory ships with unit tests (`pytest`) that verify forces, energy conservation, and SCF convergence. GitHub Actions CI runs both the stand‑alone scripts and the master notebook in headless mode to catch regressions early.

---

## License

**No explicit license has been chosen yet.** Until a LICENSE file is added, all rights are reserved by the author. If you would like to reuse or adapt this material, please open an issue or contact the repository owner to discuss licensing terms.

---

## Citation

> Christensen, D.; **Computational Chemistry Mini‑Projects Suite**, 2025.

If you employ this material in teaching or research, please cite the repository. Thanks to the open‑source community whose tools make this work possible.
