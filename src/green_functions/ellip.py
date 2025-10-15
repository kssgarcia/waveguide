import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve
from scipy.io import savemat

import src.grid_generation.gen_p1grid as gen
from fem_core.loc_matrices import locA, locB, locM

def ellip():
    """
    Finite element code for solving a 2nd-order convection-diffusion problem
    using P1 elements.
    """
    # Velocity components in the governing equation
    Vx = 1.0
    Vy = 1.0

    # Generate a triangular grid by uniformly refining a coarse grid
    nref = 3  # level of refinement
    ne, nn, p, conn, gbc = gen.gen_p1grid(nref)

    # Adjust for 0-based indexing (conn comes from gen_p1grid as 1-based)
    # conn = conn - 1

    # Plot the mesh
    plt.figure(1)
    plt.triplot(p[:, 0], p[:, 1], conn)
    plt.title('Finite Element Mesh')
    plt.pause(2)

    # Specify the exact solution and use it for Dirichlet BC
    u_ex = np.zeros(nn)
    for i in range(nn):
        u_ex[i] = EXACT(p[i, 0], p[i, 1])
        if gbc[i, 0] == 1:  # indicator for Dirichlet BC
            gbc[i, 1] = u_ex[i]

    # Initialize global matrix and RHS vector
    # Ag = lil_matrix((nn, nn))  # Use sparse matrix for efficiency
    Ag = np.zeros((nn, nn))
    b = np.zeros(nn)
    
    # Loop over the elements
    for l in range(ne):
        # Get node coordinates for this element
        j = conn[l, 0]
        x1, y1 = p[j, 0], p[j, 1]
        j = conn[l, 1]
        x2, y2 = p[j, 0], p[j, 1]
        j = conn[l, 2]
        x3, y3 = p[j, 0], p[j, 1]

        # Compute local element matrices
        elmA = locA(x1, y1, x2, y2, x3, y3)
        elmB = locB(Vx, Vy, x1, y1, x2, y2, x3, y3)
        elmM = locM(x1, y1, x2, y2, x3, y3)

        # Assemble into global system
        for i in range(3):
            i1 = conn[l, i]
            for j in range(3):
                j1 = conn[l, j]
                Ag[i1, j1] += elmA[i, j] + elmB[i, j]
                b[i1] += elmM[i, j] * SRC(p[j1, 0], p[j1, 1])

    # Convert to CSR format for efficient solving
    # Ag = Ag.tocsr()
    savemat("Ag_py.mat", {"Ag_py": Ag, "b_py": b, "gbc_py": gbc})
    

    # Impose the Dirichlet BC
    for m in range(nn):
        if gbc[m, 0] == 1:
            # b -= Ag[:, m].toarray().flatten() * gbc[m, 1]
            b -= Ag[:, m] * gbc[m, 1]
            Ag[m, :] = 0
            Ag[:, m] = 0
            Ag[m, m] = 1.0
            b[m] = gbc[m, 1]

    plt.figure(1)
    # plt.imshow(Ag.toarray())
    plt.imshow(Ag)
    plt.colorbar()
    # Solve the linear system
    # u_fem = spsolve(Ag, b)
    u_fem = solve(Ag, b)

    # Plot the numerical solution
    plt.figure(2)
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(p[:, 0], p[:, 1], u_fem, triangles=conn, cmap='viridis')
    plt.title('Numerical Solution')

    # Compare with exact solution
    max_err = np.max(np.abs(u_fem - u_ex))
    print(f'max err={max_err}')

    plt.show()

def EXACT(x, y):
    """Exact solution for our model problem."""
    return x**2 * y**2

def SRC(x, y):
    """
    Source function at the RHS of the equation:
    -(u_xx + u_yy) + v_vec * grad u = src
    For u = x^2 * y^2 and v_vec = (1,1)
    """
    return -2*y**2 - 2*x**2 + 2*x*y**2 + 2*y*x**2

# Include the previously converted functions (gen_p1grid, locA, locB, locM) here
# ...

if __name__ == "__main__":
    ellip()
