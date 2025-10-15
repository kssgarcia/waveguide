"""
Helmholtz Eigenvalue Problem Solver with Custom FEM Implementation

This script solves the 2D Helmholtz eigenvalue problem using a custom finite
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
# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from matplotlib.tri import Triangulation
import meshio

def load_mesh_meshio():
    folder = "../../meshes"
    # Read the original mesh
    mesh = meshio.read(f"{folder}/out_r_bigger.msh")

    # Extract node coordinates
    nodes = mesh.points  # Shape: (num_nodes, 3)

    # Extract element connectivity
    elements = mesh.get_cells_type("triangle")  # Shape: (num_elements, 3)

    # Extract physical group information
    print("Physical Groups:")
    if mesh.field_data:
        for name, (tag, dim) in mesh.field_data.items():
            print(f"  {name}: tag={tag}, dimension={dim}")

    # Get physical group tags for elements
    physical_tags = None
    if mesh.cell_data and "gmsh:physical" in mesh.cell_data:
        # Get physical tags for triangular elements
        cell_types = [cell.type for cell in mesh.cells]
        triangle_idx = cell_types.index("triangle") if "triangle" in cell_types else None

        if triangle_idx is not None:
            physical_tags = mesh.cell_data["gmsh:physical"][triangle_idx]
            print(f"\nPhysical tags for elements: {len(physical_tags)} elements")
            print(f"Unique physical tags: {set(physical_tags)}")

    # Alternative way to get cell data
    print(f"\nAvailable cell data keys: {list(mesh.cell_data.keys()) if mesh.cell_data else 'None'}")
    node_result = nodes[:, :2]

    return node_result, elements

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
        dirichlet_nodes = np.concatenate([top_nodes, bottom_nodes, left_nodes, right_nodes])
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

def solve_helmholtz_eigenproblem(num_eigs=6, bc_type='dirichlet_all'):
    """
    Solve Helmholtz eigenvalue problem with flexible boundary conditions

    Parameters:
    num_eigs: number of eigenvalues to compute
    bc_type: 'dirichlet_all', 'neumann_all', or 'mixed'
    """
    nodes, elements = load_mesh_meshio()

    N = len(nodes)
    print(f"Generated {N} nodes and {len(elements)} elements")

    A = lil_matrix((N, N))
    M = lil_matrix((N, N))

    # Convert to numpy arrays for efficient computation
    nodes_array = np.array(nodes) if not isinstance(nodes, np.ndarray) else nodes
    elements_array = np.array(elements) if not isinstance(elements, np.ndarray) else elements

    # Option 1: Vectorized batch assembly (more memory efficient for large meshes)
    batch_size = min(1000, len(elements_array))  # Process elements in batches

    for start_idx in range(0, len(elements_array), batch_size):
        end_idx = min(start_idx + batch_size, len(elements_array))
        element_batch = elements_array[start_idx:end_idx]

        # Get vertices for batch of elements
        vertices_batch = nodes_array[element_batch]  # Shape: (batch_size, 3, 2)

        # Compute local matrices for entire batch
        A_local_batch, M_local_batch = vectorized_local_matrices(vertices_batch)

        # Assembly for current batch
        for idx, element in enumerate(element_batch):
            i, j, k = element
            indices = np.array([i, j, k])

            for row in range(3):
                for col in range(3):
                    A[indices[row], indices[col]] += A_local_batch[idx, row, col]
                    M[indices[row], indices[col]] += M_local_batch[idx, row, col]

    # Get boundary nodes
    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = get_boundary_nodes(nodes_array, elements_array)
    boundary_nodes = (top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes)

    # print("Boundary node counts:")
    # print(f"  Top nodes: {len(top_nodes)}")
    # print(f"  Bottom nodes: {len(bottom_nodes)}")
    # print(f"  Left nodes: {len(left_nodes)}")
    # print(f"  Right nodes: {len(right_nodes)}")
    # print(f"  Obstacle nodes: {len(obstacle_nodes)}")

    # Convert to CSR format
    A = A.tocsr()
    M = M.tocsr()

    # Apply boundary conditions
    A_bc, M_bc, interior_nodes = apply_boundary_conditions(A, M, nodes_array, boundary_nodes, bc_type)

    # print(f"Applying boundary conditions: {bc_type}")
    # print(f"System size after BC: {A_bc.shape[0]} × {A_bc.shape[1]}")

    # Solve eigenvalue problem
    try:
        eigvals, eigvecs_reduced = eigsh(A_bc, k=num_eigs, M=M_bc, sigma=0, which='LM') 
    except Exception as e:
        print(f"Error in eigenvalue computation: {e}")
        # Try with different parameters
        eigvals, eigvecs_reduced = eigsh(A_bc, k=min(num_eigs, A_bc.shape[0]//2),
                                        M=M_bc, sigma=1e-6, which='LM')
    # Recover full eigenvectors
    eigvecs = np.zeros((N, len(eigvals)))
    eigvecs[interior_nodes, :] = eigvecs_reduced

    for i in range(eigvecs.shape[1]):
        eigvecs[:, i] /= np.sqrt(eigvecs[:, i].T @ M @ eigvecs[:, i])

    return eigvals, eigvecs, nodes, elements, boundary_nodes

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


def plot_all_eigenmodes(nodes, elements, eigvecs, eigvals, bc_type='dirichlet_all'):
    """
    Plot all eigenmodes in a single figure using subplots
    """
    # Ensure arrays are numpy arrays
    triangles = np.array(elements) if not isinstance(elements, np.ndarray) else elements
    nodes_array = np.array(nodes) if not isinstance(nodes, np.ndarray) else nodes
    x, y = nodes_array[:, 0], nodes_array[:, 1]

    triang = Triangulation(x, y, triangles)
    
    num_modes = eigvecs.shape[1]
    
    # Calculate grid dimensions for subplots
    if num_modes <= 3:
        rows, cols = 1, num_modes
        figsize = (5*num_modes, 4)
    elif num_modes <= 6:
        rows, cols = 2, 3
        figsize = (15, 8)
    elif num_modes <= 9:
        rows, cols = 3, 3
        figsize = (15, 12)
    else:
        rows, cols = int(np.ceil(num_modes/4)), 4
        figsize = (20, 5*rows)

    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    
    # Handle single subplot case
    if num_modes == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes if hasattr(axes, '__len__') else [axes]
    else:
        axes = axes.flatten()
    
    # Find global min/max for consistent color scaling
    vmin = np.min(eigvecs)
    vmax = np.max(eigvecs)
    
    for i in range(num_modes):
        u = eigvecs[:, i]
        
        # Plot the eigenmode
        cs = axes[i].tricontourf(triang, u, levels=50, cmap='viridis', vmin=vmin, vmax=vmax)
        axes[i].triplot(triang, 'k-', alpha=0.2, linewidth=0.1)
        
        # Add colorbar to each subplot
        cbar = fig.colorbar(cs, ax=axes[i], shrink=0.8)
        cbar.set_label("u(x, y)", fontsize=8)
        
        axes[i].set_title(f"Mode {i+1}: k² = {eigvals[i]:.4f}", fontsize=10)
        axes[i].set_xlabel("x", fontsize=8)
        axes[i].set_ylabel("y", fontsize=8)
        axes[i].set_aspect('equal')
        axes[i].tick_params(labelsize=7)
    
    # Hide unused subplots
    for i in range(num_modes, len(axes)):
        axes[i].set_visible(False)
    
    # Add overall title
    fig.suptitle(f'Helmholtz Eigenmodes - {bc_type.upper()} Boundary Conditions', 
                 fontsize=14, y=0.98)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Enhanced main execution
if __name__ == "__main__":
    print("HELMHOLTZ EIGENVALUE SOLVER - NUMPY ARRAY VERSION")
    print("="*60)

    print("\n" + "="*60)
    print("SOLVING HELMHOLTZ EIGENVALUE PROBLEM")
    print("="*60)

    num_eigs = 6

    # First, show how to access boundary nodes for Neumann BC
    print("\n" + "="*60)
    print("DEMONSTRATING BOUNDARY NODE ACCESS FOR NEUMANN BC")
    print("="*60)

    # Test different boundary conditions ['dirichlet_all', 'neumann_all', 'mixed']
    bc_type = 'dirichlet_all'

    eigvals, eigvecs, nodes, elements, boundary_nodes = \
        solve_helmholtz_eigenproblem(num_eigs, bc_type=bc_type)

    top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes

    print(f"\nComputed Eigenvalues (k²) with {bc_type}:")
    for i, val in enumerate(eigvals):
        print(f"Mode {i+1:2d}:  k² ≈ {val:.6f}")

    # Demonstrate array properties of results
    print("\nResult array properties:")
    print(f"Nodes shape: {np.array(nodes).shape}")
    print(f"Elements shape: {np.array(elements).shape}")
    print(f"Eigenvalues shape: {eigvals.shape}")
    print(f"Eigenvectors shape: {eigvecs.shape}")
    print(f"\nBoundary nodes breakdown:")
    print(f"  Top boundary: {len(top_nodes)} nodes")
    print(f"  Bottom boundary: {len(bottom_nodes)} nodes")
    print(f"  Left boundary: {len(left_nodes)} nodes")
    print(f"  Right boundary: {len(right_nodes)} nodes")
    print(f"  Obstacle boundary: {len(obstacle_nodes)} nodes")

    # Visualize boundary nodes
    plot_boundary_nodes(nodes, elements, boundary_nodes)

    # Plot all eigenmodes in a single figure using subplots
    plot_all_eigenmodes(nodes, elements, eigvecs, eigvals, bc_type)

# @author: kssgarcia (modified with numpy array optimizations)
