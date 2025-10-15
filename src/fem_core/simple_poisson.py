"""
Simple Poisson Equation Solver with Custom FEM Implementation

This script implements a basic finite element method solver for the 2D Poisson
equation from scratch. It demonstrates how to:
- Generate structured triangular meshes for the unit square domain
- Assemble global stiffness matrices using linear basis functions
- Apply homogeneous Dirichlet boundary conditions on all boundaries
- Solve the Poisson equation: -∇²u = f(x,y) with u = 0 on ∂Ω
- Visualize solutions with contour plots and 3D surfaces
- Perform convergence analysis with mesh refinement

This is an educational implementation showing the fundamental steps of FEM
for elliptic PDEs without using external FEM libraries.

@author: kssgarcia
"""
# %%
import numpy as np
import matplotlib.pyplot as plt

def generate_refined_mesh(refinement_level):
    """
    Generate a refined triangular mesh for the unit square [0,1] x [0,1]

    Parameters:
    refinement_level: int - number of subdivisions per edge (minimum 1)

    Returns:
    nodes: array of node coordinates
    elements: list of triangular elements (node indices)
    """
    n = refinement_level + 1  # number of nodes per edge

    # Generate nodes
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
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

def get_boundary_nodes(nodes, tol=1e-10):
    """
    Identify nodes on the boundary of the unit square
    """
    boundary_nodes = []
    for i, (x, y) in enumerate(nodes):
        if (abs(x) < tol or abs(x - 1) < tol or
            abs(y) < tol or abs(y - 1) < tol):
            boundary_nodes.append(i)
    return boundary_nodes

def solve_poisson_fem(refinement_level=1):
    """
    Solve -Δu = 1 on unit square with homogeneous Dirichlet BCs

    Parameters:
    refinement_level: int - mesh refinement level (1 = 2 triangles, 2 = 8 triangles, etc.)
    """

    # STEP 1: Generate refined mesh
    nodes, elements = generate_refined_mesh(refinement_level)
    N = len(nodes)

    print(f"Refinement level: {refinement_level}")
    print(f"Number of nodes: {N}")
    print(f"Number of elements: {len(elements)}")

    # STEP 2: Define f(x, y) source function
    def f(x, y):
        return 1.0  # constant source

    # STEP 3: Initialize global stiffness matrix A and RHS b
    A = np.zeros((N, N))
    b = np.zeros(N)

    # STEP 4: Local stiffness matrix and RHS computation
    def local_stiffness_and_rhs(p1, p2, p3):
        # Build local stiffness matrix using linear basis functions
        # Compute area of triangle
        mat = np.array([
            [1, p1[0], p1[1]],
            [1, p2[0], p2[1]],
            [1, p3[0], p3[1]]
        ])
        area = 0.5 * np.abs(np.linalg.det(mat))

        # Gradients of basis functions
        def grad_phi(p_a, p_b, p_c):
            # Formula from barycentric coordinates
            mat = np.array([
                [1, p_a[0], p_a[1]],
                [1, p_b[0], p_b[1]],
                [1, p_c[0], p_c[1]]
            ])
            inv = np.linalg.inv(mat)
            return inv[1:, :]  # gradients of phi_1, phi_2, phi_3

        grads = grad_phi(p1, p2, p3)  # shape (2, 3)

        # Local stiffness matrix A_local[i,j] = grad(phi_i)·grad(phi_j) * area
        A_local = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                A_local[i, j] = area * np.dot(grads[:, i], grads[:, j])

        # Local RHS vector (using centroid quadrature)
        xc = (p1[0] + p2[0] + p3[0]) / 3
        yc = (p1[1] + p2[1] + p3[1]) / 3
        fval = f(xc, yc)
        b_local = fval * area / 3 * np.ones(3)

        return A_local, b_local

    # STEP 5: Assembly loop
    for element in elements:
        i, j, k = element
        p1, p2, p3 = nodes[i], nodes[j], nodes[k]
        A_local, b_local = local_stiffness_and_rhs(p1, p2, p3)

        # Add to global system
        for a, global_a in enumerate([i, j, k]):
            b[global_a] += b_local[a]
            for b_, global_b in enumerate([i, j, k]):
                A[global_a, global_b] += A_local[a, b_]

    # STEP 6: Apply Dirichlet BCs (u = 0 on boundary nodes)
    boundary_nodes = get_boundary_nodes(nodes)
    for node in boundary_nodes:
        A[node, :] = 0
        A[:, node] = 0
        A[node, node] = 1
        b[node] = 0

    # STEP 7: Solve linear system
    u = np.linalg.solve(A, b)

    return nodes, elements, u

def plot_solution(nodes, elements, u, refinement_level):
    """
    Plot the FEM solution
    """
    from matplotlib.tri import Triangulation

    triangles = np.array(elements)
    x = nodes[:, 0]
    y = nodes[:, 1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Solution contours
    triang = Triangulation(x, y, triangles)
    cs = ax1.tricontourf(triang, u, levels=20, cmap='viridis')
    ax1.triplot(triang, 'k-', alpha=0.3, linewidth=0.5)
    fig.colorbar(cs, ax=ax1, label="u(x, y)")
    ax1.set_title(f"FEM Solution -Δu = 1 (Refinement: {refinement_level})")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect('equal')

    # Plot 2: 3D surface
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_trisurf(x, y, u, triangles=triangles, cmap='viridis', alpha=0.8)
    ax2.set_title(f"3D Surface (Refinement: {refinement_level})")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("u(x, y)")

    plt.tight_layout()
    plt.show()

    # Print some statistics
    print(f"Max solution value: {np.max(u):.6f}")
    print(f"Min solution value: {np.min(u):.6f}")

# Example usage with different refinement levels
if __name__ == "__main__":
    # Test with different refinement levels
    ref_level = 16

    nodes, elements, u = solve_poisson_fem(refinement_level=ref_level)
    plot_solution(nodes, elements, u, ref_level)

    # Compare convergence
    print(f"\n{'='*50}")
    print("Convergence Analysis:")
    print("Refinement | Nodes | Elements | Max u")
    print("-" * 40)

    nodes, elements, u = solve_poisson_fem(refinement_level=ref_level)
    print(f"{ref_level:10d} | {len(nodes):5d} | {len(elements):8d} | {np.max(u):.6f}")

# @author: kssgarcia
