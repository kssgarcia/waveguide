"""
Helmholtz Equation Solver with External Mesh using DOLFINx

This script solves the 2D Helmholtz equation using DOLFINx with an externally
generated mesh from GMSH. It demonstrates how to:
- Import and use GMSH mesh files (.msh) in DOLFINx
- Set up Dirichlet boundary conditions on waveguide walls
- Solve the Helmholtz eigenvalue problem: -∇²u - k²u = 0
- Export solutions to XDMF format for visualization

The geometry represents a waveguide with walls at y = ±1.0 where u = 0.

@author: kssgarcia
"""
# %%
from dolfinx.io import gmshio
from mpi4py import MPI
import numpy as np

# Leer la malla desde GMSH (.msh)
domain, cell_tags, facet_tags = gmshio.read_from_msh("../../meshes/out_triangles_only.msh", MPI.COMM_WORLD, gdim=3)

# %%
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
import ufl
from petsc4py import PETSc
from dolfinx.io import XDMFFile

# Espacio funcional
V = fem.functionspace(domain, ("CG", 1))

# Parámetro espectral (k^2)
# k0 = 4 * np.pi
k0 = 5

# Condición de Dirichlet (u=0 en paredes de la guía)
def walls(x):
    return np.isclose(x[1], -1.0) | np.isclose(x[1], 1.0)

walls_dofs = fem.locate_dofs_geometrical(V, walls)
bc = fem.dirichletbc(PETSc.ScalarType(0), walls_dofs, V)

# Funciones de prueba y test
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

# Forma débil del problema de Helmholtz
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx - k0**2 * ufl.inner(u, v) * ufl.dx
L = fem.Constant(domain, PETSc.ScalarType(0)) * v * ufl.dx

# Solución
problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "cg", "pc_type": "hypre"})
uh = problem.solve()
assert problem.solver.getConvergedReason() > 0
print(domain)
print(dir(uh))
print(uh.x)

# %%

with XDMFFile(domain.comm, "./solutions/hemlholtz_meshio_dolfinx.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(uh)

# @author: kssgarcia
