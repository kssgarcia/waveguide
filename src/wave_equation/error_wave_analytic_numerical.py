"""
2D Wave Equation Solver and Analytical Comparison

This script solves the 2D wave equation on a rectangular domain using:
- Finite Element Method (FEM) via a custom solver
- Analytical solution via Fourier series expansion

It demonstrates how to:
- Solve the wave equation numerically using FEM
- Compute the analytical solution at arbitrary node points
- Compare numerical and analytical solutions via L2 error norms
- Visualize solution snapshots in 2D and 3D over time

The wave equation solved is: ∂²u/∂t² = c² ∇²u on a unit square
with zero Dirichlet boundary conditions and a Gaussian initial condition.

@author: kssgarcia
"""
# %%
from simple_wave import solve_wave_fem, plot_snapshots, plot_3d_snapshots
import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
sin = np.sin
cos = np.cos

def _compute_a_mn_by_trapz(f, Nm, Nn, a=1.0, b=1.0, Nx_int=200, Ny_int=200):
    """
    Numerically compute a_mn = 4/(a*b) * \int_0^a \int_0^b f(x,y) sin(m*pi*x/a) sin(n*pi*y/b) dx dy
    using a trapezoidal rule on an Nx_int x Ny_int grid.
    Returns array shape (Nm, Nn).
    """
    x_int = np.linspace(0, a, Nx_int)
    y_int = np.linspace(0, b, Ny_int)
    X_int, Y_int = np.meshgrid(x_int, y_int, indexing='xy')  # shape (Ny_int, Nx_int)

    F = f(X_int, Y_int)  # vectorized

    a_mn = np.zeros((Nm, Nn), dtype=float)
    for m in range(1, Nm + 1):
        sx = np.sin(m * pi * X_int / a)  # shape (Ny_int, Nx_int)
        for n in range(1, Nn + 1):
            sy = np.sin(n * pi * Y_int / b)
            integrand = F * sx * sy  # elementwise
            # integrate first over x (axis=1), then over y (axis=0)
            I = np.trapz(np.trapz(integrand, x_int, axis=1), y_int, axis=0)
            a_mn[m - 1, n - 1] = 4.0 / (a * b) * I

    return a_mn

def analytical_u_at_nodes(nodes, times,
                          f=None,
                          Nm=10, Nn=10,
                          a=1.0, b=1.0, c=1.0,
                          Nx_int=200, Ny_int=200):
    """
    nodes: array (Nnodes, 2) with coordinates. Assumes coordinates are centered in [-a/2,a/2] x [-b/2,b/2].
           If your nodes already lie in [0,a]x[0,b], set nodes_centered=False and skip the shift.
    times: 1D array-like of times
    f: function f(x,y) vectorized; default is a Gaussian same as your original code
    Nm,Nn: truncation indices for Fourier series
    returns U: array shape (len(times), len(nodes))
    """
    nodes = np.asarray(nodes)
    times = np.asarray(times)
    Nnodes = nodes.shape[0]
    T = times.shape[0]

    if f is None:
        f = lambda X, Y: np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.5) ** 2))

    # --- map nodes from centered coordinates [-a/2, a/2] to domain [0,a] ---
    x_nodes = nodes[:, 0] + a / 2.0
    y_nodes = nodes[:, 1] + b / 2.0

    # --- compute coefficients a_mn once ---
    a_mn = _compute_a_mn_by_trapz(f, Nm, Nn, a=a, b=b, Nx_int=Nx_int, Ny_int=Ny_int)  # shape (Nm, Nn)

    # --- precompute mode shapes at node points ---
    # Sx[m-1, i] = sin(m*pi*x_i / a)
    Sx = np.array([np.sin(m * pi * x_nodes / a) for m in range(1, Nm + 1)])  # (Nm, Nnodes)
    Sy = np.array([np.sin(n * pi * y_nodes / b) for n in range(1, Nn + 1)])  # (Nn, Nnodes)

    # Sprod[m,n,node] = sin(m*pi*x_node/a)*sin(n*pi*y_node/b)
    Sprod = Sx[:, None, :] * Sy[None, :, :]  # shape (Nm, Nn, Nnodes)

    # angular frequencies
    m_idx = np.arange(1, Nm + 1)[:, None]
    n_idx = np.arange(1, Nn + 1)[None, :]
    omega = pi * c * np.sqrt((m_idx ** 2) / (a ** 2) + (n_idx ** 2) / (b ** 2))  # shape (Nm, Nn)

    # coef_time[m,n,t] = a_mn[m,n] * cos(omega[m,n] * t)
    times = np.asarray(times)
    coef_time = a_mn[:, :, None] * np.cos(omega[:, :, None] * times[None, :])  # (Nm, Nn, T)

    # contract over m,n to produce U[t, node]
    # tensordot: sum_{m,n} coef_time[m,n,t] * Sprod[m,n,node] -> shape (T, Nnodes)
    U = np.tensordot(coef_time, Sprod, axes=([0, 1], [0, 1]))  # (T, Nnodes)

    return U


nodes, elements, times, u_numerical = solve_wave_fem(
    refinement_level=24,  # Mesh refinement
    T=1.6,               # Final time
    dt=0.2,            # Time step
    c=1.0                # Wave speed
)
snapshot_times = np.arange(0, 1.8, 0.2)  
u_analytic = analytical_u_at_nodes(nodes, times, Nm=10, Nn=10, a=1.0, b=1.0, c=1.0)

# %% Plot solutions

# Numerical
plot_snapshots(nodes, elements, times, u_numerical, snapshot_times)
plot_3d_snapshots(nodes, elements, times, u_numerical, snapshot_times)

# Analytic
plot_snapshots(nodes, elements, times, u_analytic, snapshot_times)
plot_3d_snapshots(nodes, elements, times, u_analytic, snapshot_times)

# %% Compute error
error = u_numerical - u_analytic  # shape (T, Nnodes)
error_norm = np.linalg.norm(error, axis=1) / np.sqrt(error.shape[1])  # L2 per timestep

print("L2 error per snapshot:", error_norm)

# --- Plot error evolution over time ---
plt.figure(figsize=(6,4))
plt.plot(times, error_norm, marker='o')
plt.xlabel('Time')
plt.ylabel('L2 error')
plt.title('Error between FEM and analytical solution over time')
plt.grid(True)
plt.show()

plot_snapshots(nodes, elements, times, error, snapshot_times)
