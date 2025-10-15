"""
Wave Eigenvalue Problem Solver with Custom FEM Implementation

This script solves the 2D Wave eigenvalue problem using a custom finite
element implementation. It demonstrates how to:
- Generate structured triangular meshes for rectangular waveguide domains
- Assemble sparse stiffness and mass matrices using linear basis functions
- Apply Dirichlet boundary conditions on top and bottom walls only
- Solve the generalized eigenvalue problem: Au = λMu
- Extract and visualize the computed eigenmodes
- Compare numerical results with analytical solutions

This implementation is designed for waveguide mode analysis where waves
propagate freely in the x-direction but are constrained by walls at y = ±1.

@author: kssgarcia
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from matplotlib.tri import Triangulation

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
    n = refinement_level + 1
    x = np.linspace(x_range[0], x_range[1], n)
    y = np.linspace(y_range[0], y_range[1], n)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

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
    mat = np.array([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
    area = 0.5 * np.abs(np.linalg.det(mat))

    grads = np.linalg.inv(mat)[1:, :]
    A_local = area * grads.T @ grads # Local stiffness matrix

    M_local = area / 12 * (np.ones((3, 3)) + np.eye(3)) # Local Mass matrix
    return A_local, M_local

def solve_wave_eigenproblem(refinement_level=4, num_eigs=6):

    x_range, y_range = (-0.5, 0.5), (-0.5, 0.5)
    nodes, elements = generate_refined_mesh(refinement_level, x_range=x_range, y_range=y_range)
    N = len(nodes)

    A = lil_matrix((N, N))
    M = lil_matrix((N, N))

    for element in elements:
        i, j, k = element
        p1, p2, p3 = nodes[i], nodes[j], nodes[k]
        A_local, M_local = local_matrices(p1, p2, p3)

        for a, global_a in enumerate([i, j, k]):
            for b, global_b in enumerate([i, j, k]):
                A[global_a, global_b] += A_local[a, b]
                M[global_a, global_b] += M_local[a, b]

    # Fixed: Now get boundary nodes for all four sides
    boundary_nodes = set(get_boundary_nodes(nodes, x_range=x_range, y_range=y_range))
    all_nodes = set(range(N))
    interior_nodes = sorted(list(all_nodes - boundary_nodes))

    print(f"Total nodes: {N}")
    print(f"Boundary nodes: {len(boundary_nodes)}")
    print(f"Interior nodes: {len(interior_nodes)}")

    # Convert to CSR format
    A = A.tocsr()
    M = M.tocsr()

    # Restrict to interior system
    A_int = A[interior_nodes, :][:, interior_nodes]
    M_int = M[interior_nodes, :][:, interior_nodes]

    # Solve eigenvalue problem
    eigvals, eigvecs_reduced = eigsh(A_int, k=num_eigs, M=M_int, sigma=0, which='LM')

    # Recover full eigenvectors with zeros at boundary
    eigvecs = np.zeros((N, num_eigs))
    for i in range(num_eigs):
        eigvecs[interior_nodes, i] = eigvecs_reduced[:, i]

    return eigvals, eigvecs, nodes, elements

def plot_eigenmode(nodes, elements, u, mode_number, refinement_level):
    triangles = np.array(elements)
    x = nodes[:, 0]
    y = nodes[:, 1]

    triang = Triangulation(x, y, triangles)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Solution contours
    triang = Triangulation(x, y, triangles)
    cs = ax1.tricontourf(triang, u, levels=20, cmap='viridis')
    ax1.triplot(triang, 'k-', alpha=0.3, linewidth=0.5)
    fig.colorbar(cs, ax=ax1, label="u(x, y)")
    ax1.set_title(f"Eigenmode {mode_number+1} (Refinement: {refinement_level})")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect('equal')

    # Plot 2: 3D surface
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_trisurf(x, y, u, triangles=triangles, cmap='viridis', alpha=0.8)
    ax2.set_title(f"Eigenmode {mode_number+1} - 3D Surface")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("u(x, y)")

    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    refinement_level = 24
    num_eigs = 6

    eigvals, eigvecs, nodes, elements = solve_wave_eigenproblem(refinement_level, num_eigs)

    print(f"\n{'='*50}")
    print(f"Refinement Level: {refinement_level}")
    print(f"Computed Eigenvalues (k²):")
    for i, val in enumerate(eigvals):
        print(f"Mode {i+1:2d}:  k² ≈ {val:.6f}")

    print("\nExpected (analytical) eigenvalues for unit square:")
    print("π² ≈ 9.8696, 2π² ≈ 19.7392, 5π² ≈ 49.3480, ...")

    # Plot selected eigenmode
    # mode_to_plot = 0  # 0 = first mode
    # u = eigvecs[:, mode_to_plot]
    # plot_eigenmode(nodes, elements, u, mode_to_plot, refinement_level)

    # Plot all eigenmode
    for i in range(num_eigs):
        u = eigvecs[:, i]
        plot_eigenmode(nodes, elements, u, i, refinement_level)
# @author: kssgarcia
