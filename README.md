# Physics Animations & Numerical Optics

This repository collects small Python projects and a final course project demonstrating physics concepts through simulation and visualization.

---

## Projects

### 1. Wave Superposition
- Visualizes constructive and destructive interference of sinusoidal waves.
- Demonstrates how frequency, phase, and amplitude determine resulting waveforms.
- Includes: `wave_superposition.py`, `wave_superposition.mp4`

### 2. Particle Dynamics
- Simulates particle trajectories under Newtonian dynamics.
- Demonstrates motion from defined initial positions and velocities.
- Includes: `particle_dynamics.py`, `particle_dynamics.mp4`

### 3. Numerical Optics (Undergraduate Final Project, Physics 325)
- Full write-up: [`numerical_optics.pdf`](numerical_optics.pdf)
- Solves the **1D paraxial wave equation** using multiple numerical methods:
  - Rayleigh–Sommerfeld integral formulation
  - Leapfrog finite-difference scheme
  - Forward Euler method
  - Crank–Nicolson method (unconditionally stable)
- Discusses boundary conditions, stability issues, and computational tradeoffs.
- Python implementations included in the appendix of the report.

---

## Usage
Clone the repo and run any script:

```bash
git clone https://github.com/hschn58/Physics.git
cd Physics/Wave_Superposition
python3 wave_superposition.py
```

Scripts for the optical numerical methods are recorded at the end of [`numerical_optics.pdf`](numerical_optics.pdf)

---

## Example Outputs

| Wave Superposition | Particle Dynamics | Rayleigh–Sommerfeld Integral |
|--------------------|-------------------|-------------------------------|
| <img src="Example_Media/wave_superposition.gif" alt="Wave Superposition" width="250"> | <img src="Example_Media/2D_particle_gas.gif" alt="Particle Dynamics" width="250"> | <img src="Example_Media/rayleigh_sommerfeld_integral.png" alt="Rayleigh–Sommerfeld Integral" width="250"> |
