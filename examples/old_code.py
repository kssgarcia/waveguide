"""
Waveguide Eigenmode Analysis using DOLFINx and SLEPc

This script computes the electromagnetic modes of a rectangular waveguide using
the finite element method. It demonstrates how to:
- Create a rectangular waveguide geometry with DOLFINx
- Set up a generalized eigenvalue problem for the Helmholtz equation
- Use SLEPc for efficient eigenvalue computation
- Apply Dirichlet boundary conditions on waveguide walls
- Visualize the computed eigenmodes using PyVista

The problem solved is: ∇²E + k²E = 0 with E = 0 on the waveguide walls,
where k² are the eigenvalues representing the cutoff frequencies.

@author: kssgarcia
"""
# %%
from mpi4py import MPI
import numpy as np
import ufl
from petsc4py import PETSc
from slepc4py import SLEPc
from dolfinx import mesh, fem, plot
from dolfinx.fem import functionspace
from dolfinx.fem.petsc import assemble_matrix  # Correct import path
import pyvista, pyvistaqt

# 1. Create rectangular waveguide geometry
domain = mesh.create_rectangle(MPI.COMM_WORLD,
    [np.array([0.0, -1.0]), np.array([10.0, 1.0])],
    n=(100, 40), cell_type=mesh.CellType.triangle)
V = functionspace(domain, ("Lagrange", 2))

# 2. Apply Dirichlet BC on top & bottom
def top_bot(x):
    return np.isclose(x[1], -1.0) | np.isclose(x[1], 1.0)

bc = fem.dirichletbc(PETSc.ScalarType(0.0),
    fem.locate_dofs_geometrical(V, top_bot), V)

# 3. Define variational forms: a(u,v) and m(u,v)
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = fem.form(ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx)
m = fem.form(u * v * ufl.dx)

# Use the correct import path for assemble_matrix
A = assemble_matrix(a, bcs=[bc])
A.assemble()
M = assemble_matrix(m, bcs=[bc])
M.assemble()

# 4. Set up SLEPc eigensolver for A x = λ M x
eps = SLEPc.EPS().create(comm=MPI.COMM_WORLD)
eps.setOperators(A, M)
eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)

# 5. Solver settings: find smallest eigenvalues first
eps.setDimensions(nev=6)  # request first six modes
eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
eps.setTolerances(tol=1e-8, max_it=200)

# Configure the spectral transformation to avoid zero pivot issues
st = eps.getST()
st.setType(SLEPc.ST.Type.SHIFT)
st.setShift(0.1)  # Small positive shift to avoid zero eigenvalues

# 6. Solve
eps.solve()
n_conv = eps.getConverged()
print(f"Number of converged eigenpairs: {n_conv}")

# 7. Extract & visualize first few modes
plt_data =  plot.vtk_mesh(V)

for i in range(min(n_conv, 5)):
    eig_vec, _ = A.getVecs()
    eig_val = eps.getEigenpair(i, eig_vec, None)
    print(f"Mode {i}, λ = {eig_val.real:.6f}")

    u_sol = fem.Function(V)
    u_sol.vector.setArray(eig_vec.array.real)

    grid = pyvista.UnstructuredGrid(plt_data[0], plt_data[1], domain.geometry.x)
    grid.point_data["u"] = u_sol.x.array.real

    p = pyvistaqt.BackgroundPlotter()
    p.add_mesh(grid, scalars="u", cmap="viridis", show_edges=True)
    p.view_xy()
    p.add_text(f"Mode {i}, λ={eig_val.real:.3f}", font_size=14)
    p.show()

# @author: kssgarcia

# %%
cells = plt_data[0]
types = plt_data[1]
x = domain.geometry.x
grid = pyvista.UnstructuredGrid(cells, types, x)

# Attach the cells tag data to the pyvita grid

# Create a plotter with two subplots, and add mesh tag plot to the
# first sub-window
subplotter = pyvista.Plotter()
cells, types, x = plot.vtk_mesh(domain)

sub_grid = pyvista.UnstructuredGrid(cells, types, x)
subplotter.add_text("Subset of mesh", font_size=14, color="black", position="upper_edge")
subplotter.add_mesh(sub_grid, show_edges=True, edge_color="black")

if pyvista.OFF_SCREEN:
    subplotter.screenshot(
        "2D_markers.png", transparent_background=transparent, window_size=[2 * figsize, figsize]
    )
else:
    subplotter.show()
# %%
