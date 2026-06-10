# Euler vs RK4 Projectile Motion Visualization

## Overview

For this project I wanted to show methods for simulating projectile motion under constant acceleration, usually gravitational. The simulation compares three solutions:

1. Euler's Method
2. Fourth-Order Runge-Kutta (RK4)
3. Exact Analytical Solution

The resulting trajectories were exported as VTK files and visualized in VisIt. The visualization highlights how numerical error accumulates in Euler's Method while RK4 remains nearly identical to the exact solution. To run this simulation you must have Python and csv, os, and math modules installed.

---


## Running the Simulation

Run the program:

```bash
python3 euler_step.py
```

Example inputs:
The terminal will show examples for every input the user needs to provide so if they are confused there is guidance.

```text
Initial x position: 0
Initial y position: 0

Initial x velocity: 5
Initial y velocity: 15

x acceleration: 0
y acceleration: -9.8

Timestep size: 0.1
Number of steps: 100
```

The program generates:

```text
vtk_output/
    particle0000.vtk
    particle0001.vtk
    ...
    particles.visit
```

---

## Visualization in VisIt

1. Open `particles.visit`.
2. Add a Pseudocolor plot using `solution_type`.
3. Add a Vector plot using `velocity_vector`.
4. Draw the plots.
5. Press play, the animation should start.

---

## Results

The visualization demonstrates that Euler's Method accumulates numerical error over time, causing its trajectory to diverge from the exact solution. RK4 produces significantly higher accuracy and closely matches the analytical solution throughout the simulation. I know this is true because the RK4 solution overlaps the exact solution, proving they are the same and the Euler's method is slightly incorrect. 

---
