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

### 4. Anharmonic Oscillator (Undergraduate Final Project, Physics 311)
- Notebook: [`Anharmonic_Oscillator.ipynb`](Anharmonic_Oscillator.ipynb)
- Studies the **nonlinear driven pendulum**, including:
  - Lagrangian → Hamiltonian → exact nonlinear equation of motion
  - Amplitude-dependent period (elliptic integral derivation)
  - Parametric driving at ω = ω₀ and why divergence does not occur
  - Stabilization of the inverted pendulum at high-frequency drive
  - Chaotic motion and sensitivity to initial conditions
- Combines analytic derivations with numerical simulation (RK4) and visualizations.

---

## Requirements

The scripts require Python 3.x and the following external libraries:

- [NumPy](https://numpy.org/)  
- [Matplotlib](https://matplotlib.org/)  
- [SciPy](https://scipy.org/)  

You can install them with:

```bash
pip install numpy matplotlib scipy
```

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
