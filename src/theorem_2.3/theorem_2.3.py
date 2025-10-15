# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from matplotlib.tri import Triangulation
from scipy.sparse.linalg import eigsh
import src.grid_generation.gen_mesh_center_y as gen
import dipole_theorem as dip

def find_boundary_edges(elements):
    """
    Find boundary edges by identifying edges that belong to only one triangle
    Returns a set of boundary edges (each edge as a sorted tuple)
    """
    edge_count = {}

    for element in elements:
        # Get the three edges of this triangle
        edges = [
            tuple(sorted([element[0], element[1]])),
            tuple(sorted([element[1], element[2]])),
            tuple(sorted([element[2], element[0]]))
        ]

        for edge in edges:
            edge_count[edge] = edge_count.get(edge, 0) + 1

    # Boundary edges appear in exactly one triangle
    boundary_edges = {edge for edge, count in edge_count.items() if count == 1}

    return boundary_edges

def get_boundary_nodes(nodes, elements, tol=1e-10):
    """
    Identify nodes on different boundaries for applying different boundary conditions

    Returns:
    top_nodes: list of node indices on top boundary (y = b = 1.0)
    bottom_nodes: list of node indices on bottom boundary (y = -b = -1.0)
    left_nodes: list of node indices on left boundary (x = -L = -5.0)
    right_nodes: list of node indices on right boundary (x = L = 5.0)
    obstacle_nodes: list of node indices on obstacle boundary (circular hole)
    """
    nodes_array = np.array(nodes) if not isinstance(nodes, np.ndarray) else nodes
    x, y = nodes_array[:, 0], nodes_array[:, 1]

    # Find all boundary edges using mesh connectivity
    boundary_edges = find_boundary_edges(elements)

    # Get all boundary nodes
    boundary_nodes_set = set()
    for edge in boundary_edges:
        boundary_nodes_set.add(edge[0])
        boundary_nodes_set.add(edge[1])

    boundary_nodes = np.array(list(boundary_nodes_set))

    # Domain boundaries based on mesh geometry
    # From mesh_create.py: L = 5.0, b = 1.0, r = 0.1, a = 0.0
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)

    # Detect outer boundary nodes among boundary nodes
    boundary_x = x[boundary_nodes]
    boundary_y = y[boundary_nodes]

    on_left = (np.abs(boundary_x - x_min) < tol) & (boundary_y > y_min + tol) & (boundary_y < y_max - tol)
    on_right = (np.abs(boundary_x - x_max) < tol) & (boundary_y > y_min + tol) & (boundary_y < y_max - tol)
    on_bottom = np.abs(boundary_y - y_min) < tol   # y = -1.0
    on_top = np.abs(boundary_y - y_max) < tol      # y = 1.0

    # Get node indices for outer boundaries
    top_nodes = boundary_nodes[on_top]
    bottom_nodes = boundary_nodes[on_bottom]
    left_nodes = boundary_nodes[on_left]
    right_nodes = boundary_nodes[on_right]

    # Obstacle nodes are boundary nodes that are NOT on outer boundary
    on_outer_boundary = on_left | on_right | on_bottom | on_top
    on_obstacle = ~on_outer_boundary
    obstacle_nodes = boundary_nodes[on_obstacle]

    return top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes

def apply_boundary_conditions(A, M, nodes, boundary_nodes, bc_type='mixed'):
    """
    Apply boundary conditions to the system matrices

    Parameters:
    A, M: sparse matrices (stiffness and mass)
    nodes: node coordinates
    boundary_nodes: tuple of (top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes)
    bc_type: 'dirichlet_all', 'neumann_all', or 'mixed'

    Returns:
    A_bc, M_bc: modified matrices
    interior_nodes: nodes where the solution is computed
    """
    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes
    N = A.shape[0]

    if bc_type == 'dirichlet_all':
        # Apply Dirichlet BC (u=0) on all boundaries
        all_boundary_nodes = np.concatenate([top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes])
        dirichlet_nodes = np.unique(all_boundary_nodes)

        all_nodes = np.arange(N)
        interior_mask = np.ones(N, dtype=bool)
        interior_mask[dirichlet_nodes] = False
        interior_nodes = all_nodes[interior_mask]

        # Restrict to interior system
        A_bc = A[np.ix_(interior_nodes, interior_nodes)]
        M_bc = M[np.ix_(interior_nodes, interior_nodes)]

    elif bc_type == 'neumann_all':
        # Apply Neumann BC (natural BC, no modification needed)
        interior_nodes = np.arange(N)
        A_bc = A
        M_bc = M

    elif bc_type == 'mixed':
        # Example: Dirichlet on outer boundaries, Neumann on obstacle
        dirichlet_nodes = np.concatenate([top_nodes, bottom_nodes, left_nodes, right_nodes])
        dirichlet_nodes = np.unique(dirichlet_nodes)

        all_nodes = np.arange(N)
        interior_mask = np.ones(N, dtype=bool)
        interior_mask[dirichlet_nodes] = False
        interior_nodes = all_nodes[interior_mask]

        # Restrict to interior + obstacle boundary system
        A_bc = A[np.ix_(interior_nodes, interior_nodes)]
        M_bc = M[np.ix_(interior_nodes, interior_nodes)]

    elif bc_type == 'neuman_dirichlet':
        # Example: Dirichlet on outer boundaries, Neumann on obstacle
        dirichlet_nodes = np.concatenate([left_nodes, right_nodes])
        dirichlet_nodes = np.unique(dirichlet_nodes)

        all_nodes = np.arange(N)
        interior_mask = np.ones(N, dtype=bool)
        interior_mask[dirichlet_nodes] = False
        interior_nodes = all_nodes[interior_mask]

        # Restrict to interior + obstacle boundary system
        A_bc = A[np.ix_(interior_nodes, interior_nodes)]
        M_bc = M[np.ix_(interior_nodes, interior_nodes)]

    else:
        raise ValueError(f"Unknown boundary condition type: {bc_type}")

    return A_bc, M_bc, interior_nodes

def vectorized_local_matrices(vertices_batch):
    """
    Compute local matrices for multiple elements at once
    vertices_batch: (num_elements, 3, 2) array of vertex coordinates
    """
    num_elements = vertices_batch.shape[0]

    # Add column of ones for transformation matrices
    ones = np.ones((num_elements, 3, 1))
    mat_batch = np.concatenate([ones, vertices_batch], axis=2)  # Shape: (num_elements, 3, 3)

    # Compute areas for all elements
    det_batch = np.linalg.det(mat_batch)
    areas = 0.5 * np.abs(det_batch)  # Shape: (num_elements,)

    # Compute gradients for all elements
    inv_mat_batch = np.linalg.inv(mat_batch)
    grads_batch = inv_mat_batch[:, 1:, :]  # Shape: (num_elements, 2, 3)

    # Compute local stiffness matrices
    A_local_batch = np.zeros((num_elements, 3, 3))
    for i in range(num_elements):
        A_local_batch[i] = areas[i] * grads_batch[i].T @ grads_batch[i]

    # Compute local mass matrices (same for all elements of same area)
    M_local_batch = np.zeros((num_elements, 3, 3))
    base_mass = np.ones((3, 3)) + np.eye(3)
    for i in range(num_elements):
        M_local_batch[i] = areas[i] / 12 * base_mass

    return A_local_batch, M_local_batch

def solve_system(A, M, nodes, elements, boundary_nodes,
                     left_nodes, right_nodes,
                     b=1.0,
                     k=1.0, bc_type='mixed', verbose=True):
    """
    Solve (K - B(k)) u = k^2 M u  using Newton on (k,u).
    Inputs A, M are full assembled matrices (csr), before DtN applied.
    left_nodes, right_nodes: arrays of node indices for left/right DtN application.
    Returns converged (k, u) on the reduced system (after applying bc_type).
    """

    # ensure complex dtype
    A = A.tocsr().astype(complex)
    M = M.tocsr().astype(complex)

    A_mod = A.copy().tolil()
    A_mod = A_mod.tocsr()

    # Reduce according to bc_type (this removes Dirichlet nodes)
    A_bc, M_bc, interior_nodes = apply_boundary_conditions(A_mod, M, nodes, boundary_nodes, bc_type)
    A_bc = A_bc.tocsr().astype(complex)
    M_bc = M_bc.tocsr().astype(complex)

    nred = A_bc.shape[0]
    if nred == 0:
        raise RuntimeError("No unknowns after applying BCs.")

    # If first iteration, get initial eigenpair from linearized problem for starting u
    eigvals, eigvecs = eigsh(A_bc, k=10, M=M_bc, sigma=k**2)

    return interior_nodes, A_bc, M_bc, eigvals, eigvecs

def solve_helmholtz_eigenproblem(nodes, elements, bc_type='mixed', b=1.0, k_guess=1.0, max_iter=10, tol=1e-6):
    N = len(nodes)
    nodes_array = np.array(nodes)
    elements_array = np.array(elements)

    # Assemble stiffness/mass
    A = lil_matrix((N, N))
    M = lil_matrix((N, N))
    batch_size = min(1000, len(elements_array))
    for start_idx in range(0, len(elements_array), batch_size):
        end_idx = min(start_idx + batch_size, len(elements_array))
        element_batch = elements_array[start_idx:end_idx]
        vertices_batch = nodes_array[element_batch]
        A_local_batch, M_local_batch = vectorized_local_matrices(vertices_batch)
        for idx, element in enumerate(element_batch):
            i, j, k = element
            indices = np.array([i, j, k])
            for row in range(3):
                for col in range(3):
                    A[indices[row], indices[col]] += A_local_batch[idx, row, col]
                    M[indices[row], indices[col]] += M_local_batch[idx, row, col]
    A = A.tocsr()
    M = M.tocsr()

    # Boundary nodes
    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = get_boundary_nodes(nodes_array, elements_array)
    boundary_nodes = (top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes)

    interior_nodes, A_bc, M_bc, eigvals, eigvecs = solve_system(
        A, M, nodes_array, elements_array, boundary_nodes,
        left_nodes=left_nodes, right_nodes=right_nodes,
        k=k_guess, bc_type=bc_type, verbose=True
    )

    return interior_nodes, boundary_nodes, A_bc, M_bc, eigvals, eigvecs

def plot_boundary_nodes(nodes, elements, boundary_nodes):
    """
    Visualize the different types of boundary nodes
    """
    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes

    nodes_array = np.array(nodes) if not isinstance(nodes, np.ndarray) else nodes
    triangles = np.array(elements) if not isinstance(elements, np.ndarray) else elements
    x, y = nodes_array[:, 0], nodes_array[:, 1]

    triang = Triangulation(x, y, triangles)

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))  # Wider to give more space

    # Plot mesh
    ax.triplot(triang, 'k-', alpha=0.3, linewidth=0.5)

    # Plot different boundary types with different colors
    if len(top_nodes) > 0:
        ax.scatter(x[top_nodes], y[top_nodes], c='red', s=30, label='Top boundary', zorder=5)
    if len(bottom_nodes) > 0:
        ax.scatter(x[bottom_nodes], y[bottom_nodes], c='blue', s=30, label='Bottom boundary', zorder=5)
    if len(left_nodes) > 0:
        ax.scatter(x[left_nodes], y[left_nodes], c='green', s=30, label='Left boundary', zorder=5)
    if len(right_nodes) > 0:
        ax.scatter(x[right_nodes], y[right_nodes], c='orange', s=30, label='Right boundary', zorder=5)
    if len(obstacle_nodes) > 0:
        ax.scatter(x[obstacle_nodes], y[obstacle_nodes], c='purple', s=30, label='Obstacle boundary', zorder=5)

    ax.set_title("Boundary Node Classification")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Place legend outside the plot, on the right
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()

def plot_eigenmode_on_ax(nodes, elements, u_full, ax, title=None):
    triangles = np.array(elements)
    nodes_array = np.array(nodes)
    x, y = nodes_array[:, 0], nodes_array[:, 1]
    triang = Triangulation(x, y, triangles)
    u_real = np.real(u_full)

    cs = ax.tricontourf(triang, u_real, levels=50, cmap='jet')
    ax.triplot(triang, 'k-', alpha=0.2, linewidth=0.1)

    if title is not None:
        ax.set_title(title, fontsize=10)
    ax.set_xlabel("x", fontsize=8)
    ax.set_ylabel("y", fontsize=8)
    ax.set_aspect('equal')
    ax.tick_params(labelsize=7)

    return cs

def Lambda_n(n, b):
    return np.pi**2*n**2/(4*b**2)

# Enhanced main execution
if __name__ == "__main__":
    epsilon = 1e-2
    bc_type = 'neuman_dirichlet'
    L=5.0 # longitud de la guia
    b=1.0 # mitad de la altura de la guia
    xc=0.0 # posicion x del obstaculo
    yc=0.0 # posicion y del obstaculo
    beta=1
    # r=0.2 # radio
    lc=0.01 # Controla la finura del mallado
    k_guess = 3

    mu = dip.dipole(epsilon, beta, xc, yc)
    Lambda1 = Lambda_n(1,b)
    Lambda2 = Lambda_n(2,b)
    sigma = epsilon**2*np.pi**3/b**3 * mu
    k2_analytic = Lambda2-sigma**2
    print(np.sqrt(k2_analytic))

    det_list = []
    kb_list = np.arange(1, 4, 0.05)

    # nodes, elements = gen.mesh_with_obstacle_center(L=L, b=b, xc=xc, yc=yc, r=r, lc=lc)
    nodes, elements = gen.mesh_with_parametric_obstacle_even_x(L=L, b=b, xc=xc, yc=yc, lc=lc, beta=beta, n_points=50, scale=epsilon, plot=True)
    interior_nodes, boundary_nodes, A_bc, M_bc, eigvals, eigvecs_reduced = \
        solve_helmholtz_eigenproblem(nodes, elements, bc_type=bc_type, b=b, k_guess=k_guess)

    # Solve eigenvalue problem
    num_eigs = eigvals.shape[0]
    eigvecs = np.zeros((len(nodes), num_eigs))
    for i in range(num_eigs):
        eigvecs[interior_nodes, i] = eigvecs_reduced[:, i]

    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes

    # Create a figure with 10 subplots (2 rows x 5 columns)
    fig, axes = plt.subplots(5, 2, figsize=(20, 15))
    axes = axes.flatten()  # Flatten to make indexing easier

    for i in range(num_eigs):
        print(f"Mode {i}:")
        print(f"  k = {eigvals[i]**0.5:.6f}")
        print(f"  k² = {eigvals[i]:.6f}")
        print(f"  k² analytic = {k2_analytic:.6f}")

        u_full = eigvecs[:, i].real

        # Plot on the corresponding subplot
        ax = axes[i]
        cs = plot_eigenmode_on_ax(nodes, elements, u_full, ax,
                                title=f"Mode {i}, k={eigvals[i]**0.5:.4f}")

        # Add colorbar to each subplot
        cbar = fig.colorbar(cs, ax=ax)
        cbar.set_label("Re(u)", fontsize=8)
        cbar.ax.tick_params(labelsize=6)

    plt.tight_layout()
    plt.show()

# @author: kssgarcia
