# Vibrational Frequencies of N₂ and O₂ Molecules

This repository contains code and data for simulating the vibrational motion of diatomic molecules (N₂ and O₂) using Symplectic Euler and Velocity Verlet time integration schemes. The project includes the production of separation and energy data, computation of wavenumbers, and comparison between purely vibrational motion and rotational/spin motion.

## Table of Contents

1. [Simulation Code](#simulation-code)
2. [Report Analysis](#report-analysis)
3. [Files Included](#files-included)
4. [Example Usage](#example-usage)


## Simulation Code

The main simulation code is **`simulate_particle.py`**, which performs a time integration of two particles interacting via the Morse potential. It produces:

* **Time-dependent separation** of the two atoms (plotted in Python)
* **Total energy** of the system (plotted in Python)
* **Wavenumber** (numerically calculated from the separation peaks and printed in the console)

The simulation supports two integration schemes:

1. **Symplectic Euler**
2. **Velocity Verlet**

The particles themselves are represented by the **`Particle3D`** class (in `particle3D.py`), which handles:

* Position and velocity updates
* Kinetic energy calculation
* Momentum calculation

### Input Data

Input files specify particle properties and initial conditions. Included are:

* `oxygen.dat` / `oxygen_spin.dat`
* `nitrogen.dat` / `nitrogen_spin.dat`

The spin files introduce initial velocities in multiple dimensions (x and y), simulating rotational effects, whereas the non-spin files simulate purely vibrational motion (velocity only in x).

### Output Data

Running a simulation produces `.dat` files containing:

* Time
* Separation between particles
* Total energy

These `.dat` files are produced by the code when running simulations; they are not included in the repository. Examples:

* `verlet.dat` (Velocity Verlet output)
* `euler.dat` (Symplectic Euler output)

The Python code also produces plots of separation vs time and total energy vs time (and x and y positions of the particles vs time) for immediate visualisation.


## Report Analysis

Further analysis for the report was performed in Excel (not included in this repository). This includes:

* **Time-step estimation**: identifying the optimal δt for ≤0.5% relative frequency error
* **Energy inaccuracy assessment**: computing ΔE = max(E) - min(E)
* **Further analysis of the simulation results**: e.g., comparisons of spin vs no spin, multi-dimensional vs 1D motion, and preparation of figures for the report

### Key Observations from the Report

* Velocity Verlet consistently showed lower energy inaccuracy than Symplectic Euler.

* Wavenumbers for purely vibrational motion are higher than for spin/rotational motion.

* Comparison to experimental values:

  * Nitrogen: ~2359 cm⁻¹
  * Oxygen: ~1580 cm⁻¹

* Simulations with spin produce slightly lower wavenumbers due to multi-dimensional motion increasing oscillation periods.

### Using Output for Analysis

The `.dat` files produced by the simulation can be plotted directly in Python or imported into Excel for further analysis:

* **Separation vs time** (Python plot)
* **Total energy vs time** (Python plot)
*  **Peak detection in Python**: used to calculate vibrational periods and wavenumbers from the separation data
* Comparisons between purely vibrational and spin simulations (`*_spin.dat`)

## Files Included

| File                   | Description                                             |
| ---------------------- | ------------------------------------------------------- |
| `simulate_particle.py` | Main simulation script                                  |
| `particle3D.py`        | Class definition for Particle3D                         |
| `oxygen.dat`           | Initial conditions for O₂ (vibrational only)            |
| `oxygen_spin.dat`      | Initial conditions for O₂ (including rotational motion) |
| `nitrogen.dat`         | Initial conditions for N₂ (vibrational only)            |
| `nitrogen_spin.dat`    | Initial conditions for N₂ (including rotational motion) |

**Note:** `.dat` output files (e.g., `verlet.dat`, `euler.dat`) are produced when running simulations and are not included in the repository.


## Example Usage

Run a simulation for oxygen without spin using Velocity Verlet:

```bash
%run simulate_particle.py oxygen.dat verlet verlet.dat
```

Or using Symplectic Euler:

```bash
%run simulate_particle.py oxygen.dat euler euler.dat
```

For spin simulations, replace the input file:

```bash
%run simulate_particle.py oxygen_spin.dat verlet verlet_spin.dat
```

The generated `.dat` files can then be used for further analysis or plotting.
