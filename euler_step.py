# 2D Euler vs RK4 vs Exact Solution and VTK output for VisIt

import csv
import os
import math


def rk4_step(x, y, vx, vy, ax, ay, dt):
    def deriv(state):
        sx, sy, svx, svy = state
        return [svx, svy, ax, ay]

    s = [x, y, vx, vy]

    k1 = deriv(s)
    k2 = deriv([s[i] + 0.5 * dt * k1[i] for i in range(4)])
    k3 = deriv([s[i] + 0.5 * dt * k2[i] for i in range(4)])
    k4 = deriv([s[i] + dt * k3[i] for i in range(4)])

    return [
        s[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        for i in range(4)
    ]


def write_vtk_file(filename, euler_traj, rk4_traj, exact_traj, rk4_vx, rk4_vy, t):
    r = 0.09

    total_trail_points = len(euler_traj) + len(rk4_traj) + len(exact_traj)
    total_points = total_trail_points + 12

    euler_particle_start = total_trail_points
    rk4_particle_start = euler_particle_start + 4
    exact_particle_start = rk4_particle_start + 4

    xe, ye, euler_error = euler_traj[-1]
    xr, yr, rk4_error = rk4_traj[-1]
    xt, yt = exact_traj[-1]

    with open(filename, "w") as file:
        file.write("# vtk DataFile Version 3.0\n")
        file.write("Euler vs RK4 vs Exact projectile motion\n")
        file.write("ASCII\n")
        file.write("DATASET POLYDATA\n")

        file.write(f"POINTS {total_points} float\n")

        for x, y, error in euler_traj:
            file.write(f"{x} {y} 0\n")

        for x, y, error in rk4_traj:
            file.write(f"{x} {y} 0\n")

        for x, y in exact_traj:
            file.write(f"{x} {y} 0\n")

        # Euler, RK4, and Exact particle squares
        for cx, cy in [(xe, ye), (xr, yr), (xt, yt)]:
            file.write(f"{cx-r} {cy-r} 0\n")
            file.write(f"{cx+r} {cy-r} 0\n")
            file.write(f"{cx+r} {cy+r} 0\n")
            file.write(f"{cx-r} {cy+r} 0\n")

        # Three trajectory lines
        if len(euler_traj) > 1:
            file.write(
                f"LINES 3 "
                f"{len(euler_traj) + 1 + len(rk4_traj) + 1 + len(exact_traj) + 1}\n"
            )

            file.write(str(len(euler_traj)))
            for i in range(len(euler_traj)):
                file.write(f" {i}")
            file.write("\n")

            rk4_offset = len(euler_traj)
            file.write(str(len(rk4_traj)))
            for i in range(len(rk4_traj)):
                file.write(f" {rk4_offset + i}")
            file.write("\n")

            exact_offset = len(euler_traj) + len(rk4_traj)
            file.write(str(len(exact_traj)))
            for i in range(len(exact_traj)):
                file.write(f" {exact_offset + i}")
            file.write("\n")

        # Three particle squares
        file.write("POLYGONS 3 15\n")
        file.write(
            f"4 {euler_particle_start} {euler_particle_start+1} "
            f"{euler_particle_start+2} {euler_particle_start+3}\n"
        )
        file.write(
            f"4 {rk4_particle_start} {rk4_particle_start+1} "
            f"{rk4_particle_start+2} {rk4_particle_start+3}\n"
        )
        file.write(
            f"4 {exact_particle_start} {exact_particle_start+1} "
            f"{exact_particle_start+2} {exact_particle_start+3}\n"
        )

        file.write(f"POINT_DATA {total_points}\n")

        # 1 = Euler, 2 = RK4, 3 = Exact
        file.write("SCALARS solution_type float 1\n")
        file.write("LOOKUP_TABLE default\n")
        for _ in euler_traj:
            file.write("1\n")
        for _ in rk4_traj:
            file.write("2\n")
        for _ in exact_traj:
            file.write("3\n")
        for _ in range(4):
            file.write("1\n")
        for _ in range(4):
            file.write("2\n")
        for _ in range(4):
            file.write("3\n")

        file.write("SCALARS error float 1\n")
        file.write("LOOKUP_TABLE default\n")
        for x, y, error in euler_traj:
            file.write(f"{error}\n")
        for x, y, error in rk4_traj:
            file.write(f"{error}\n")
        for _ in exact_traj:
            file.write("0\n")
        for _ in range(4):
            file.write(f"{euler_error}\n")
        for _ in range(4):
            file.write(f"{rk4_error}\n")
        for _ in range(4):
            file.write("0\n")

        file.write("SCALARS time float 1\n")
        file.write("LOOKUP_TABLE default\n")
        for _ in range(total_points):
            file.write(f"{t}\n")

        # RK4 velocity vector.
        # Put the arrow only on the RK4 particle square points.
        # All other points get a zero vector so they do not make arrows.
        file.write("VECTORS velocity_vector float\n")

        for _ in euler_traj:
            file.write("0 0 0\n")

        for _ in rk4_traj:
            file.write("0 0 0\n")

        for _ in exact_traj:
            file.write("0 0 0\n")

        # Euler square: no arrow
        for _ in range(4):
            file.write("0 0 0\n")

        # RK4 square: visible arrow
        for _ in range(4):
            file.write(f"{rk4_vx} {rk4_vy} 0\n")

        # Exact square: no arrow
        for _ in range(4):
            file.write("0 0 0\n")


def run_simulation():
    print("=== Euler vs RK4 vs Exact Solution ===")
    print("Example values:")
    print("  Initial x position: 0")
    print("  Initial y position: 0")
    print("  Initial x velocity: 5")
    print("  Initial y velocity: 15")
    print("  x acceleration: 0")
    print("  y acceleration: -9.8")
    print("  Timestep size: 0.05")
    print("  Number of steps: 120\n")

    x0 = float(input("Initial x position (example: 0): "))
    y0 = float(input("Initial y position (example: 0): "))
    vx0 = float(input("Initial x velocity (example: 5): "))
    vy0 = float(input("Initial y velocity (example: 15): "))
    ax = float(input("x acceleration (example: 0): "))
    ay = float(input("y acceleration (example: -9.8): "))
    dt = float(input("Timestep size (example: 0.05): "))
    num_steps = int(input("Number of steps (example: 120): "))

    output_dir = "vtk_output"
    os.makedirs(output_dir, exist_ok=True)

    results = []
    vtk_files = []

    euler_traj = []
    rk4_traj = []
    exact_traj = []

    xe, ye = x0, y0
    vxe, vye = vx0, vy0

    xr, yr = x0, y0
    vxr, vyr = vx0, vy0

    t = 0.0

    print("\nResults:")
    print("time\teuler_error\trk4_error")

    for step in range(num_steps + 1):
        xt = x0 + vx0 * t + 0.5 * ax * t * t
        yt = y0 + vy0 * t + 0.5 * ay * t * t

        euler_error = math.sqrt((xe - xt) ** 2 + (ye - yt) ** 2)
        rk4_error = math.sqrt((xr - xt) ** 2 + (yr - yt) ** 2)

        euler_traj.append((xe, ye, euler_error))
        rk4_traj.append((xr, yr, rk4_error))
        exact_traj.append((xt, yt))

        print(f"{t:.2f}\t{euler_error:.6f}\t{rk4_error:.6f}")

        results.append([
            step, t,
            xe, ye,
            xr, yr,
            xt, yt,
            euler_error, rk4_error
        ])

        vtk_name = f"particle{step:04d}.vtk"
        vtk_path = os.path.join(output_dir, vtk_name)
        write_vtk_file(vtk_path, euler_traj, rk4_traj, exact_traj, vxr, vyr, t)
        vtk_files.append(vtk_name)

        # Euler update
        vxe = vxe + ax * dt
        vye = vye + ay * dt
        xe = xe + vxe * dt
        ye = ye + vye * dt

        # RK4 update
        xr, yr, vxr, vyr = rk4_step(xr, yr, vxr, vyr, ax, ay, dt)

        t = t + dt

    with open("euler_results.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "step", "time",
            "euler_x", "euler_y",
            "rk4_x", "rk4_y",
            "exact_x", "exact_y",
            "euler_error", "rk4_error"
        ])
        writer.writerows(results)

    visit_path = os.path.join(output_dir, "particles.visit")
    with open(visit_path, "w") as file:
        for vtk_file in vtk_files:
            file.write(vtk_file + "\n")

    print("\nSaved simulation data to euler_results.csv")
    print(f"Saved VTK timestep files in {output_dir}/")
    print(f"Saved VisIt file: {visit_path}")


run_simulation()