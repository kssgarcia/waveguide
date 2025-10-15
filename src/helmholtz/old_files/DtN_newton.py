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

    on_left = np.abs(boundary_x - x_min) < tol     # x = -5.0
    on_right = np.abs(boundary_x - x_max) < tol    # x = 5.0
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

def transverse_modes(y_nodes, N_modes, b):
    phi = np.zeros((len(y_nodes), N_modes))
    for n in range(1, N_modes+1):
        phi[:, n-1] = np.sqrt(1/b) * np.sin(n * np.pi * (y_nodes + b)/(2*b))
    return phi

def assemble_DtN(nodes, boundary_nodes, k, N_modes=5, b=1.0, sign=1):
    """
    sign=1 for right boundary (outgoing wave)
    sign=-1 for left boundary (opposite direction)
    """
    y_nodes = nodes[boundary_nodes,1]
    phi = transverse_modes(y_nodes, N_modes, b)
    beta = np.sqrt(-(np.arange(1, N_modes+1)*np.pi/(2*b))**2 + k**2 + 0*1j)

    num_nodes = len(boundary_nodes)
    B = lil_matrix((num_nodes, num_nodes), dtype=complex)
    for n in range(N_modes):
        B += sign * 1j * beta[n] * (phi[:, n][:, None] @ phi[:, n][None, :])
    return B

def complex_beta(k, lambda_m):
    return np.lib.scimath.sqrt(k**2 - lambda_m)

def assemble_dBdk(nodes, boundary_nodes, k, N_modes=5, b=1.0, sign=1):
    """
    Build derivative dB/dk (size = num_boundary_nodes x num_boundary_nodes) for given k.
    Matches the assemble_DtN convention: B = sum_m i*beta_m * phi_m phi_m^T (collocation).
    Thus dB/dk = sum_m i * (k / beta_m) * phi_m phi_m^T
    """
    y_nodes = nodes[boundary_nodes, 1]
    # same transverse_modes used in assemble_DtN
    phi = transverse_modes(y_nodes, N_modes, b)  # shape (num_nodes_on_boundary, N_modes)

    # transverse eigenvalues for Dirichlet on y in [-b,b] with your definition:
    # You used sin(n * pi * (y + b) / (2b)), so lambda_n = (n*pi/(2b))^2
    lam = (np.arange(1, N_modes+1) * np.pi / (2*b))**2
    beta = complex_beta(k, lam)  # shape (N_modes,)

    # guard tiny beta -> avoid zero division: if beta very small, add small eps
    eps = 1e-12
    beta_safe = beta.copy()
    near_zero = np.abs(beta_safe) < eps
    if np.any(near_zero):
        beta_safe[near_zero] = beta_safe[near_zero] + eps

    # build dB/dk (dense then convert to sparse) — small boundary size so ok
    num_nodes = len(boundary_nodes)
    dB = np.zeros((num_nodes, num_nodes), dtype=complex)
    for m in range(N_modes):
        coeff = 1j * (k / beta_safe[m]) * sign
        vec = phi[:, m][:, None]
        dB += coeff * (vec @ vec.T)

    return csr_matrix(dB)

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

def newton_nep_solve(A, M, nodes, elements, boundary_nodes, 
                     left_nodes, right_nodes,
                     N_modes=10, b=1.0,
                     max_newton_iter=20, tol_res=1e-8, tol_step=1e-8,
                     k_init=1.0, bc_type='mixed', verbose=True):
    """
    Solve (K - B(k)) u = k^2 M u  using Newton on (k,u).
    Inputs A, M are full assembled matrices (csr), before DtN applied.
    left_nodes, right_nodes: arrays of node indices for left/right DtN application.
    Returns converged (k, u) on the reduced system (after applying bc_type).
    """

    # ensure complex dtype
    A = A.tocsr().astype(complex)
    M = M.tocsr().astype(complex)

    # initial linear solve to get starting u and lambda (optional but helpful)
    k = k_init
    if verbose:
        print(f"[Newton] initial k guess = {k_init:.6g}")

    # Newton main loop
    for it in range(max_newton_iter):
        # Build A_mod = A - B_right(k) - B_left(k)
        A_mod = A.copy().tolil()

        # assemble B at current k and subtract into A_mod (same style as your code)
        B_right = assemble_DtN(nodes, right_nodes, k, N_modes, b, sign=1)
        if len(right_nodes) > 0:
            A_mod[np.ix_(right_nodes, right_nodes)] -= B_right

        B_left = assemble_DtN(nodes, left_nodes, k, N_modes, b, sign=-1)
        if len(left_nodes) > 0:
            A_mod[np.ix_(left_nodes, left_nodes)] -= B_left

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
            eigvals, eigvecs = eigs(A_bc, k=1, M=M_bc, sigma=(k**2 if k!=0 else 0.0))
        except Exception:
            eigvals, eigvecs = eigsh(A_bc.real, k=min(6, nred-1), M=M_bc.real, sigma=(k**2 if k!=0 else 0.0))

        # pick eigenpair with smallest residual (or smallest magnitude)
        u_red = eigvecs[:, 0]
        # normalize u: u^H M u = 1
        denom = (u_red.conj().T @ (M_bc @ u_red))
        if np.abs(denom) == 0:
            raise RuntimeError("Initial eigenvector has zero M-norm")
        u_red = u_red / np.sqrt(denom)

        # Newton inner iterations: we treat one Newton step per outer iteration (but we will iterate)
        # Compute residual r = F(k) u = (A_bc - k^2 M_bc) u
        r = (A_bc @ u_red) - (k**2) * (M_bc @ u_red)

        res_norm = np.linalg.norm(r)
        constraint = (u_red.conj().T @ (M_bc @ u_red)) - 1.0
        if verbose:
            print(f"[Newton it {it}] k = {k:.12g}, residual norm = {res_norm:.3e}, constraint = {constraint:.3e}")

        if res_norm < tol_res and abs(constraint) < 1e-10:
            # converged
            if verbose:
                print(f"Converged: k = {k:.12g}, residual = {res_norm:.3e}")
            # build full u vector (size N) with zeros on removed nodes
            u_full = np.zeros(A.shape[0], dtype=complex)
            u_full[interior_nodes] = u_red
            return k, u_full, (A_bc, M_bc, interior_nodes)

        # compute g = dF/dk * u = (-dB/dk - 2k M) u
        dB_right = assemble_dBdk(nodes, right_nodes, k, N_modes, b, sign=1) if len(right_nodes)>0 else csr_matrix((A.shape[0], A.shape[0]), dtype=complex)
        dB_left = assemble_dBdk(nodes, left_nodes, k, N_modes, b, sign=-1) if len(left_nodes)>0 else csr_matrix((A.shape[0], A.shape[0]), dtype=complex)

        # place dB into reduced system (same reduction used above)
        dB_full = csr_matrix((A.shape[0], A.shape[0]), dtype=complex)
        if len(right_nodes) > 0:
            dB_full[right_nodes[:, None], right_nodes] += dB_right
        if len(left_nodes) > 0:
            dB_full[left_nodes[:, None], left_nodes] += dB_left

        # reduce dB to interior
        dB_red = dB_full[interior_nodes][:, interior_nodes]

        # g_red = (-dB_red - 2k M_bc) @ u_red
        g_red = ( -dB_red - 2.0*k * M_bc ) @ u_red

        # Build Newton block system:
        # [A_bc - k^2 M_bc    g_red] [du] = -r
        # [ (u^H M)_row       0   ] [dk]   = -(u^H M u - 1)
        # where (u^H M)_row is row vector shape (1, nred)
        Aop = (A_bc - k**2 * M_bc).tocsr()
        # block matrix using scipy.sparse.bmat
        top_left = Aop
        top_right = csr_matrix(g_red.reshape(-1,1))   # nred x 1
        bottom_left = csr_matrix((u_red.conj().T @ M_bc).reshape(1, -1))  # 1 x nred
        bottom_right = csr_matrix((1,1), dtype=complex)

        # compose block matrix (nred+1 x nred+1)
        block = bmat([[top_left, top_right],
                      [bottom_left, bottom_right]], format='csc')

        # RHS
        rhs_top = -r
        rhs_bottom = -((u_red.conj().T @ (M_bc @ u_red)) - 1.0)   # scalar
        rhs = np.concatenate([rhs_top, np.array([rhs_bottom])], axis=0)

        # solve linear system for [du; dk]
        # Use spsolve (direct)
        try:
            sol = spsolve(block, rhs)
        except Exception as e:
            raise RuntimeError(f"Failed to solve Newton linear system: {e}")

        du = sol[:nred]
        dk = sol[-1]

        # Update
        u_red = u_red + du
        k = k + dk

        # renormalize u_red to satisfy u^H M u = 1
        denom = (u_red.conj().T @ (M_bc @ u_red))
        if np.abs(denom) == 0:
            raise RuntimeError("Zero norm during Newton update")
        u_red = u_red / np.sqrt(denom)

        # check step size
        if np.linalg.norm(du) < tol_step and abs(dk) < tol_step:
            if verbose:
                print(f"Newton step small; stopping at it={it}")
            break

    # If we exit loop without meeting tol, return last iterate with warning
    print("WARNING: Newton did not converge to requested tolerance")
    u_full = np.zeros(A.shape[0], dtype=complex)
    u_full[interior_nodes] = u_red
    return k, u_full, (A_bc, M_bc, interior_nodes)

def solve_helmholtz_eigenproblem(bc_type='mixed', N_modes=5, b=1.0, max_iter=10, tol=1e-6):
    # nodes, elements = utils.load_mesh_meshio("mesh_with_hole")
    nodes, elements = gen.mesh_with_obstacle(lc=0.1)
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

    k_final, u_full, (A_bc, M_bc, interior_nodes) = newton_nep_solve(
        A, M, nodes_array, elements_array, boundary_nodes,
        left_nodes=left_nodes, right_nodes=right_nodes,
        N_modes=N_modes, b=b, max_newton_iter=20, tol_res=1e-8,
        k_init=1.0, bc_type=bc_type, verbose=True
    )

    k2_final = k_final**2
    print(f"Final result: k = {k_final}, k^2 = {k2_final}")

    return k_final, k2_final, u_full, nodes_array, elements_array, boundary_nodes

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

# Enhanced main execution
if __name__ == "__main__":
    bc_type = 'mixed'
    N_modes = 12
    b = 5.0
    k_final, k2_final, u_full, nodes, elements, boundary_nodes = \
        solve_helmholtz_eigenproblem(bc_type=bc_type, N_modes=N_modes, b=b)

    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes

    print(f"Converged k^2 (final) = {k2_final}")
    # Print boundary node counts
    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes
    print("Boundary counts:", len(top_nodes), len(bottom_nodes), len(left_nodes), len(right_nodes), len(obstacle_nodes))
    # Visualize
    plot_boundary_nodes(nodes, elements, boundary_nodes)
    plot_eigenmode(nodes, elements, u_full, title=f"Mode Re(u), k^2 = {np.real(k2_final):.6g}")

# @author: kssgarcia (modified with numpy array optimizations)

# %%
