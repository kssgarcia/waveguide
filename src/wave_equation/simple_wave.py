"""
Time-Dependent Wave Equation Solver with Custom FEM Implementation

This script implements a comprehensive finite element method solver for the 2D
time-dependent wave equation from scratch. It demonstrates how to:

- Generate structured triangular meshes for rectangular domains
- Assemble consistent mass and stiffness matrices using linear basis functions
- Apply homogeneous Dirichlet boundary conditions (fixed edges)
- Solve the time-dependent wave equation: ∂²u/∂t² - c²∇²u = f(x,y,t)
- Implement stable time integration using the Newmark method
- Handle initial conditions for displacement and velocity
- Visualize wave propagation with:
  * 2D contour plots at multiple time snapshots
  * 3D surface plots showing wave evolution
  * Animated visualizations of wave propagation

Key Features:
- Newmark time integration (β=0.25, γ=0.5) for unconditional stability
- Consistent mass matrix formulation for accurate dynamics
- Customizable source functions and initial conditions
- Flexible mesh refinement for convergence studies
- Comprehensive visualization tools for wave analysis

@author: kssgarcia
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def generate_refined_mesh(refinement_level, x_range=(0, 1), y_range=(0, 1)):
    """
    Generate a refined triangular mesh for a rectangular domain

    Parameters:
    refinement_level: int - number of subdivisions per edge (minimum 1)
    x_range: tuple - (x_min, x_max) for the domain in x direction
    y_range: tuple - (y_min, y_max) for the domain in y direction

    Returns:
    nodes: array of node coordinates
    elements: list of triangular elements (node indices)
    """
    n = refinement_level + 1  # number of nodes per edge

    # Generate nodes
    x = np.linspace(x_range[0], x_range[1], n)
    y = np.linspace(y_range[0], y_range[1], n)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

    # Generate triangular elements
    elements = []
    for i in range(refinement_level):
        for j in range(refinement_level):
            # Lower triangle
            v0 = i * n + j
            v1 = (i + 1) * n + j
            v2 = i * n + (j + 1)
            elements.append([v0, v1, v2])

            # Upper triangle
            v0 = (i + 1) * n + j
            v1 = (i + 1) * n + (j + 1)
            v2 = i * n + (j + 1)
            elements.append([v0, v1, v2])

    return nodes, elements

def get_boundary_nodes(nodes, x_range=(0, 1), y_range=(0, 1), tol=1e-10):
    """
    Identify nodes on all four boundaries of the rectangular domain

    Parameters:
    nodes: array of node coordinates
    x_range: tuple - (x_min, x_max) for the domain in x direction
    y_range: tuple - (y_min, y_max) for the domain in y direction
    tol: float - tolerance for boundary detection

    Returns:
    list of boundary node indices (all four boundaries)
    """
    x_min, x_max = x_range
    y_min, y_max = y_range

    boundary_nodes = []
    for i, (x, y) in enumerate(nodes):
        # Check if node is on any of the four boundaries
        if (abs(x - x_min) < tol or  # Left boundary
            abs(x - x_max) < tol or  # Right boundary
            abs(y - y_min) < tol or  # Bottom boundary
            abs(y - y_max) < tol):   # Top boundary
            boundary_nodes.append(i)

    return boundary_nodes

def local_matrices(p1, p2, p3):
    # Build local matrices using linear basis functions
    mat = np.array([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
    area = 0.5 * np.abs(np.linalg.det(mat))

    # Stiffness matrix (gradient terms)
    grads = np.linalg.inv(mat)[1:, :]
    K_local = area * grads.T @ grads

    # Mass matrix (consistent mass matrix)
    M_local = area / 12 * (np.ones((3, 3)) + np.eye(3))

    return K_local, M_local


def solve_wave_fem(refinement_level=8, T=2.0, dt=0.01, c=1.0):
    """
    Solve the wave equation ∂²u/∂t² - c²∇²u = f(x,y,t) using FEM

    Parameters:
    refinement_level: int - mesh refinement level
    T: float - final time
    dt: float - time step size
    c: float - wave speed
    """

    # Source function (time-dependent)
    def f(x, y, t):
        # Gaussian source that oscillates
        r_sq = (x - 0.0)**2 + (y - 0.0)**2
        return 10 * np.exp(-50 * r_sq) * np.sin(5 * np.pi * t) * (t < 0.5)

    # Initial conditions
    def u0(x, y):
        r_sq = x**2 + y**2
        return np.exp(-20 * r_sq)  # Gaussian centered at (0,0) with width parameter 20

    def v0(x, y):
        # Initial velocity
        return np.zeros_like(x)

    # STEP 1: Generate refined mesh
    x_range, y_range = (-0.5, 0.5), (-0.5, 0.5)
    nodes, elements = generate_refined_mesh(refinement_level, x_range=x_range, y_range=y_range)
    N = len(nodes)

    print(f"Refinement level: {refinement_level}")
    print(f"Number of nodes: {N}")
    print(f"Number of elements: {len(elements)}")
    print(f"Time step: {dt}, Final time: {T}")

    # STEP 2: Initialize global mass matrix M and stiffness matrix K
    M = np.zeros((N, N))  # Mass matrix
    K = np.zeros((N, N))  # Stiffness matrix

    # STEP 4: Assembly loop
    for element in elements:
        i, j, k = element
        p1, p2, p3 = nodes[i], nodes[j], nodes[k]
        K_local, M_local = local_matrices(p1, p2, p3)

        # Add to global system
        for a, global_a in enumerate([i, j, k]):
            for b, global_b in enumerate([i, j, k]):
                K[global_a, global_b] += K_local[a, b]
                M[global_a, global_b] += M_local[a, b]

    # STEP 5: Apply Dirichlet BCs (u = 0 on boundary nodes)
    boundary_nodes = get_boundary_nodes(nodes, x_range=x_range, y_range=y_range)
    for node in boundary_nodes:
        K[node, :] = 0
        K[:, node] = 0
        K[node, node] = 1
        M[node, :] = 0
        M[:, node] = 0
        M[node, node] = 1

    # STEP 6: Set initial conditions
    x_coords = nodes[:, 0]
    y_coords = nodes[:, 1]

    u_curr = u0(x_coords, y_coords)  # u^n
    v_curr = v0(x_coords, y_coords)  # velocity at t^n

    # Apply boundary conditions to initial conditions
    for node in boundary_nodes:
        u_curr[node] = 0
        v_curr[node] = 0

    # STEP 7: Time integration using Newmark method
    # Parameters for Newmark method (average acceleration)
    beta = 0.25
    gamma = 0.5

    # Effective matrices for Newmark method
    dt2 = dt * dt
    M_eff = M + beta * dt2 * c**2 * K

    # Storage for solution history
    time_steps = int(T / dt) + 1
    times = np.linspace(0, T, time_steps)
    u_history = np.zeros((time_steps, N))
    u_history[0, :] = u_curr

    # Initial acceleration (at t=0)
    def compute_force_vector(t):
        b = np.zeros(N)
        for element in elements:
            i, j, k = element
            p1, p2, p3 = nodes[i], nodes[j], nodes[k]

            # Centroid of element
            xc = (p1[0] + p2[0] + p3[0]) / 3
            yc = (p1[1] + p2[1] + p3[1]) / 3
            area = 0.5 * np.abs(np.linalg.det(
                np.array([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
            ))

            fval = f(xc, yc, t)
            b_local = fval * area / 3 * np.ones(3)

            for a, global_a in enumerate([i, j, k]):
                b[global_a] += b_local[a]

        # Apply boundary conditions
        for node in boundary_nodes:
            b[node] = 0

        return b

    F0 = compute_force_vector(0)
    a_curr = np.linalg.solve(M, F0 - c**2 * K @ u_curr)

    # Apply boundary conditions to acceleration
    for node in boundary_nodes:
        a_curr[node] = 0

    print("Starting time integration...")

    # Time stepping loop
    for n in range(1, time_steps):
        t = times[n]

        # Compute force vector at current time
        F = compute_force_vector(t)

        # Predict displacement and velocity
        u_pred = u_curr + dt * v_curr + 0.5 * dt2 * (1 - 2*beta) * a_curr
        v_pred = v_curr + dt * (1 - gamma) * a_curr

        # Solve for acceleration at next time step
        rhs = F - c**2 * K @ u_pred
        a_next = np.linalg.solve(M_eff, rhs)

        # Correct displacement and velocity
        u_next = u_pred + beta * dt2 * a_next
        v_next = v_pred + gamma * dt * a_next

        # Apply boundary conditions
        for node in boundary_nodes:
            u_next[node] = 0
            v_next[node] = 0
            a_next[node] = 0

        # Store solution and update
        u_history[n, :] = u_next
        u_curr = u_next.copy()
        v_curr = v_next.copy()
        a_curr = a_next.copy()

        if time_steps >= 10 and n % (time_steps // 10) == 0:
            print(f"Progress: {100*n/time_steps:.1f}% (t = {t:.2f})")

    return nodes, elements, times, u_history

def plot_snapshots(nodes, elements, times, u_history, snapshot_times=[0.0, 0.5, 1.0, 1.5]):
    """
    Plot solution at specific time snapshots
    """
    from matplotlib.tri import Triangulation

    triangles = np.array(elements)
    x = nodes[:, 0]
    y = nodes[:, 1]
    triang = Triangulation(x, y, triangles)

    # Find closest time indices
    snapshot_indices = []
    for t_snap in snapshot_times:
        idx = np.argmin(np.abs(times - t_snap))
        snapshot_indices.append(idx)

    # Calculate number of rows and columns for subplots
    n_snapshots = len(snapshot_times)
    n_cols = 3  # 3 columns
    n_rows = (n_snapshots + n_cols - 1) // n_cols  # Ceiling division

    # Create subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()

    vmin = np.min(u_history)
    vmax = np.max(u_history)

    for i, idx in enumerate(snapshot_indices):
        ax = axes[i]
        cs = ax.tricontourf(triang, u_history[idx], levels=20, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.triplot(triang, 'k-', alpha=0.2, linewidth=0.3)
        ax.set_title(f"t = {times[idx]:.3f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect('equal')
        plt.colorbar(cs, ax=ax)

    # Hide unused subplots
    for i in range(len(snapshot_indices), len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.show()

def plot_3d_snapshots(nodes, elements, times, u_history, snapshot_times=[0.0, 0.5, 1.0, 1.5]):
    """
    Plot 3D surface plots at specific time snapshots
    """
    triangles = np.array(elements)
    x = nodes[:, 0]
    y = nodes[:, 1]

    # Find closest time indices
    snapshot_indices = []
    for t_snap in snapshot_times:
        idx = np.argmin(np.abs(times - t_snap))
        snapshot_indices.append(idx)

    # Calculate number of rows and columns for subplots
    n_snapshots = len(snapshot_times)
    n_cols = 3  # 3 columns
    n_rows = (n_snapshots + n_cols - 1) // n_cols  # Ceiling division

    # Create 3D subplots
    fig = plt.figure(figsize=(18, 6*n_rows))

    vmin = np.min(u_history)
    vmax = np.max(u_history)

    for i, idx in enumerate(snapshot_indices):
        ax = fig.add_subplot(n_rows, n_cols, i+1, projection='3d')

        # Create 3D surface plot
        surf = ax.plot_trisurf(x, y, u_history[idx], triangles=triangles,
                              cmap='RdBu_r', vmin=vmin, vmax=vmax,
                              alpha=0.9, linewidth=0.1, edgecolors='k')

        ax.set_title(f"Wave Solution at t = {times[idx]:.3f}", fontsize=12)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u(x, y, t)")

        # Set consistent z-limits for all subplots
        ax.set_zlim(vmin, vmax)

        # Add colorbar
        fig.colorbar(surf, ax=ax, shrink=0.8, aspect=10)

    plt.tight_layout()
    plt.show()

def create_3d_animation(nodes, elements, times, u_history, save_file=None, elevation=30, azimuth=45):
    """
    Create a 3D animation of the wave propagation
    """

    triangles = np.array(elements)
    x = nodes[:, 0]
    y = nodes[:, 1]

    # Set up the figure and 3D axis
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Determine z-limits
    vmin = np.min(u_history)
    vmax = np.max(u_history)

    # Initial plot
    surf = ax.plot_trisurf(x, y, u_history[0], triangles=triangles,
                          cmap='RdBu_r', vmin=vmin, vmax=vmax,
                          alpha=0.9, linewidth=0.1, edgecolors='k')

    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_zlim(vmin, vmax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x, y, t)")
    ax.view_init(elev=elevation, azim=azimuth)

    # Add colorbar
    cbar = plt.colorbar(surf, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label("u(x, y, t)")

    def animate_3d(frame):
        ax.clear()

        surf = ax.plot_trisurf(x, y, u_history[frame], triangles=triangles,
                              cmap='RdBu_r', vmin=vmin, vmax=vmax,
                              alpha=0.9, linewidth=0.1, edgecolors='k')

        ax.set_xlim(np.min(x), np.max(x))
        ax.set_ylim(np.min(y), np.max(y))
        ax.set_zlim(vmin, vmax)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u(x, y, t)")
        ax.set_title(f"3D Wave Propagation - t = {times[frame]:.3f}")
        ax.view_init(elev=elevation, azim=azimuth)

        return surf,

    # Create animation with fewer frames for better performance
    skip = max(1, len(times) // 50)  # Show about 50 frames max
    frames = range(0, len(times), skip)

    ani = FuncAnimation(fig, animate_3d, frames=frames, interval=100, blit=False, repeat=True)

    if save_file:
        ani.save(save_file, writer='pillow', fps=10)
        print(f"3D Animation saved to {save_file}")

    plt.show()
    return ani

# Example usage
if __name__ == "__main__":
    # Solve the wave equation
    print("Solving 2D Wave Equation with FEM...")
    nodes, elements, times, u_history = solve_wave_fem(
        refinement_level=24,  # Mesh refinement
        T=1.6,               # Final time
        dt=0.005,            # Time step
        c=1.0                # Wave speed
    )

    # Print solution statistics
    print("Solution Statistics:")
    print(f"Max displacement: {np.max(np.abs(u_history)):.6f}")
    print(f"Final time: {times[-1]:.3f}")

    # Create 2D snapshots from 0 to 1.6 with 0.2 step
    snapshot_times = np.arange(0, 1.8, 0.2)  # [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
    print(f"Creating 2D snapshots at times: {snapshot_times}")
    plot_snapshots(nodes, elements, times, u_history, snapshot_times)

    # Create 3D snapshots from 0 to 1.6 with 0.2 step
    print(f"Creating 3D snapshots at times: {snapshot_times}")
    plot_3d_snapshots(nodes, elements, times, u_history, snapshot_times)

    # Create 3D animation (uncomment to see animation)
    print("Creating 3D animation...")
    create_3d_animation(nodes, elements, times, u_history, save_file='wave_3d_animation.gif')
    # create_3d_animation(nodes, elements, times, u_history, save_file="wave_3d_animation.gif")

# @author: kssgarcia
