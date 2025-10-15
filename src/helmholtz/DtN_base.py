# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from matplotlib.tri import Triangulation
from scipy.sparse import csr_matrix, bmat
from scipy.sparse.linalg import spsolve, eigs, eigsh
import src.helmholtz.utils.utils as utils
import src.grid_generation.gen_mesh_center_y as gen
from scipy.sparse.linalg import svds

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
        dirichlet_nodes = np.concatenate([top_nodes, bottom_nodes])
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

def vectorized_local_matrices_1d_optimized(vertices_batch):
    """
    Optimized version using vectorized operations throughout
    """
    # Compute element lengths
    lengths = vertices_batch[:, 1] - vertices_batch[:, 0]  # Shape: (num_elements,)

    # Check for valid elements
    if np.any(lengths <= 0):
        raise ValueError("All elements must have positive length")

    # Vectorized stiffness matrices
    # A_local = (1/h) * [[1, -1], [-1, 1]]
    base_stiffness = np.array([[1, -1], [-1, 1]])
    A_local_batch = (1.0 / lengths)[:, None, None] * base_stiffness[None, :, :]

    # Vectorized mass matrices
    # M_local = (h/6) * [[2, 1], [1, 2]]
    base_mass = np.array([[2, 1], [1, 2]])
    M_local_batch = (lengths / 6.0)[:, None, None] * base_mass[None, :, :]

    return A_local_batch, M_local_batch

def assemble_fem_1d_problem_corrected(nodes, elements, side_nodes, k, N_modes, b):
    """
    Correctly implemented 1D FEM problem assembly
    """
    nodes = np.asarray(nodes)
    side_nodes = np.asarray(side_nodes, dtype=int)

    # Local index mapping for boundary nodes
    boundary_index_map = {int(g): idx for idx, g in enumerate(side_nodes)}
    nb = len(side_nodes)

    # Find boundary edges
    all_edges = set()
    for tri in elements:
        if len(tri) == 2:  # 1D elements have 2 nodes
            edge = tuple(sorted([tri[0], tri[1]]))
            all_edges.add(edge)
        else:  # If still using triangular elements, extract edges
            for a, b in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]:
                edge = tuple(sorted([a, b]))
                all_edges.add(edge)

    # Keep only edges on this side
    side_set = set(side_nodes)
    elements_edges = [e for e in all_edges if e[0] in side_set and e[1] in side_set]

    # Initialize global matrices
    n_dof = nb
    A = lil_matrix((n_dof, n_dof))
    M = lil_matrix((n_dof, n_dof))

    # STEP 1: Assemble ALL elements normally (don't skip any!)
    for (ga, gb) in elements_edges:
        # Get node coordinates
        node_coord = nodes[[ga, gb], 1]  # Extract y-coordinates

        # Ensure proper ordering
        if node_coord[1] < node_coord[0]:
            ga, gb = gb, ga
            node_coord = node_coord[::-1]

        # Compute local matrices
        A_local, M_local = vectorized_local_matrices_1d_optimized(np.array([node_coord]))
        A_elem = A_local[0]
        M_elem = M_local[0]

        # Map to local indices
        local_ga = boundary_index_map[ga]
        local_gb = boundary_index_map[gb]
        dofs = [local_ga, local_gb]

        # Assemble into global matrices (ALWAYS do this)
        for i in range(2):
            for j in range(2):
                A[dofs[i], dofs[j]] += A_elem[i, j]
                M[dofs[i], dofs[j]] += M_elem[i, j]

    # STEP 2: Create RHS vector
    if isinstance(b, (int, float)):
        rhs = np.full(n_dof, b, dtype=float)
    else:
        rhs = np.zeros(n_dof, dtype=float)

    # STEP 3: Apply Dirichlet BCs AFTER assembly
    # Find nodes with min/max y coordinates
    y_coords = nodes[side_nodes, 1]
    min_y_idx = np.argmin(y_coords)  # Local index
    max_y_idx = np.argmax(y_coords)  # Local index

    bc_nodes = [min_y_idx, max_y_idx]
    bc_values = [0.0, 0.0]  # Homogeneous Dirichlet BCs (u = 0 at boundaries)

    # Apply BCs by modifying assembled matrices
    for node_idx, bc_val in zip(bc_nodes, bc_values):
        # Zero out row and column
        A[node_idx, :] = 0
        A[:, node_idx] = 0
        M[node_idx, :] = 0
        M[:, node_idx] = 0

        # Set diagonal and RHS
        A[node_idx, node_idx] = 1.0
        M[node_idx, node_idx] = 0.0  # Mass matrix should be 0 for BC nodes
        rhs[node_idx] = bc_val

    # Convert to CSR format
    A = A.tocsr()
    M = M.tocsr()

    # For eigenvalue problem, we need to solve the generalized eigenvalue problem
    # A * u = lambda * M * u
    # But we need to handle the zero mass matrix entries from BCs

    # Find free (non-BC) nodes for eigenvalue computation
    free_nodes = np.setdiff1d(np.arange(n_dof), bc_nodes)

    if len(free_nodes) > 0:
        # Extract submatrices for free nodes only
        A_free = A[np.ix_(free_nodes, free_nodes)]
        M_free = M[np.ix_(free_nodes, free_nodes)]

        # Solve eigenvalue problem on free nodes
        if N_modes > len(free_nodes):
            N_modes = len(free_nodes) - 1

        eigvals_free, eigvecs_free = eigsh(A_free, k=N_modes, M=M_free, which='SM')

        # Reconstruct full eigenvectors
        eigvecs = np.zeros((n_dof, N_modes))
        eigvecs[free_nodes, :] = eigvecs_free
        eigvals = eigvals_free
        for j in range(eigvecs.shape[1]):
            norm = np.sqrt(eigvecs[:, j].T @ (M @ eigvecs[:, j]))
            eigvecs[:, j] /= norm
    else:
        # If no free nodes, return empty arrays
        eigvals = np.array([])
        eigvecs = np.zeros((n_dof, 0))

    return A, M, eigvals, eigvecs

def assemble_DtN_corrected(nodes, elements, side_nodes, k, N_modes=10, b=1.0, eps=1e-14):
    """
    Assemble DtN discrete matrix define in point 4 of the chat
    """
    A, M, eigvals, eigvecs = assemble_fem_1d_problem_corrected(nodes, elements, side_nodes, k, N_modes, b)

    if len(eigvals) == 0:
        # Return zero matrix if no eigenvalues found
        nb = len(side_nodes)
        return csr_matrix((nb, nb), dtype=complex)

    V = eigvecs
    alpha_n = np.lib.scimath.sqrt(k**2 - eigvals)
    diag_i_alpha = 1j * alpha_n
    VtM = (V.T @ M)
    DtN_term = V @ (diag_i_alpha[:, None] * VtM)

    return csr_matrix(DtN_term, dtype=complex)

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
                     N_modes=10, b=1.0,
                     k=1.0, bc_type='mixed', verbose=True, eignum=1):
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

    # assemble B at current k and subtract into A_mod
    B_right = assemble_DtN_corrected(nodes, elements, right_nodes, k, N_modes, b)
    if len(right_nodes) > 0:
        A_mod[np.ix_(right_nodes, right_nodes)] -= B_right

    B_left = assemble_DtN_corrected(nodes, elements, left_nodes, k, N_modes, b)
    if len(left_nodes) > 0:
        A_mod[np.ix_(left_nodes, left_nodes)] += B_left

    A_mod = A_mod.tocsr()

    # Reduce according to bc_type (this removes Dirichlet nodes)
    A_bc, M_bc, interior_nodes = apply_boundary_conditions(A_mod, M, nodes, boundary_nodes, bc_type)
    A_bc = A_bc.tocsr().astype(complex)
    M_bc = M_bc.tocsr().astype(complex)

    nred = A_bc.shape[0]
    if nred == 0:
        raise RuntimeError("No unknowns after applying BCs.")

    # If first iteration, get initial eigenpair from linearized problem for starting u
    # Solve small generalized eigenproblem (closest to k^2)
    try:
        eigvals, eigvecs = eigs(A_bc, k=eignum, M=M_bc, sigma=k**2)
    except Exception:
        eigvals, eigvecs = eigsh(A_bc.real, k=min(6, nred-1), M=M_bc.real, sigma=(k**2 if k!=0 else 0.0))

    u_red = eigvecs[:, 0]
    # normalize: u^H M u = 1
    u_red /= np.sqrt(np.real(u_red.conj().T @ (M_bc @ u_red)))

    # map back to full vector
    u_full = np.zeros(A.shape[0], dtype=complex)
    u_full[interior_nodes] = u_red

    return u_full, interior_nodes, A_bc, M_bc, eigvals, eigvecs

def solve_helmholtz_eigenproblem(nodes, elements, bc_type='mixed', N_modes=5, b=1.0, k_guess=1.0, max_iter=10, tol=1e-6, eignum=1):
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

    u_full, interior_nodes, A_bc, M_bc, eigvals, eigvecs = solve_system(
        A, M, nodes_array, elements_array, boundary_nodes,
        left_nodes=left_nodes, right_nodes=right_nodes,
        k=k_guess, bc_type=bc_type, verbose=True, eignum=eignum
        )

    return u_full, interior_nodes, boundary_nodes, A_bc, M_bc, eigvals, eigvecs

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

def plot_eigenmode(nodes, elements, u_full, title=None):
    triangles = np.array(elements)
    nodes_array = np.array(nodes)
    x, y = nodes_array[:, 0], nodes_array[:, 1]
    triang = Triangulation(x, y, triangles)
    u_real = np.real(u_full)
    fig, ax = plt.subplots(1, 1, figsize=(8,4))
    cs = ax.tricontourf(triang, u_real, levels=50, cmap='jet')
    ax.triplot(triang, 'k-', alpha=0.2, linewidth=0.1)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label("Re(u)")
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect('equal')
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

# Enhanced main execution
if __name__ == "__main__":
    # nodes, elements = utils.load_mesh_meshio("mesh_with_hole")
    # nodes, elements = gen.mesh_with_perturbed_hole(lc=0.5, plot=True)
    bc_type = 'mixed'
    L=5.0 # longitud de la guia
    b=1.0 # mitad de la altura de la guia
    xc=0.0 # posicion x del obstaculo
    yc=0.0 # posicion y del obstaculo
    beta=1
    r=0.125 # radio
    lc=0.02 # Controla la finura del mallado
    k_guess = 6.1313
    N_modes = 20

    det_list = []
    kb_list = np.arange(1, 4, 0.05)

    nodes, elements = gen.mesh_with_two_obstacle(L=L, b=b, xc=0.6/2, yc=yc, r=r, lc=lc, plot=True)

    u_full, interior_nodes, boundary_nodes, A_bc, M_bc, eigvals, eigvecs_reduced = \
        solve_helmholtz_eigenproblem(nodes, elements, bc_type=bc_type, N_modes=N_modes, b=b, k_guess=k_guess, eignum=10)

    num_eigs = eigvals.shape[0]
    eigvecs = np.zeros((len(nodes), num_eigs))
    for i in range(num_eigs):
        eigvecs[interior_nodes, i] = eigvecs_reduced[:, i]

    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes

    # Create a figure with 10 subplots (2 rows x 5 columns)
    fig, axes = plt.subplots(num_eigs, 1, figsize=(20, 15))
    axes = axes.flatten()  # Flatten to make indexing easier

    for i in range(num_eigs):
        print(f"Mode {i}:")
        print(f"  k = {eigvals[i]**0.5:.6f}")
        print(f"  k² = {eigvals[i].real:.6f}")

        u_full = eigvecs[:, i].real

        # Plot on the corresponding subplot
        ax = axes[i]
        cs = plot_eigenmode_on_ax(nodes, elements, u_full, ax,
                                title=f"Mode {i}, k={eigvals[i].real**0.5:.4f}")

        # Add colorbar to each subplot
        cbar = fig.colorbar(cs, ax=ax)
        cbar.set_label("Re(u)", fontsize=8)
        cbar.ax.tick_params(labelsize=6)

    plt.tight_layout()
    plt.show()


# @author: kssgarcia
