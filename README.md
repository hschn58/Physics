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

### 3. Numerical Optics (Physics 325 Final Project)
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
python wave_superposition.py
Animations of wave superposition and particle kinetics
```

---

## Example Outputs

### Wave Superposition
![Wave Superposition](2D_particle_gas.gif)

### Particle Dynamics
![Particle Dynamics](particle_dynamics.gif)
