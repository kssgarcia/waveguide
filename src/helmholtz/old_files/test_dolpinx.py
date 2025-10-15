# %%
import dolfinx.mesh
import meshio
from mpi4py import MPI
from dolfinx.io import gmshio

# Path to your mesh
mesh_file = "../../meshes/out.msh"

# Load GMSH mesh into dolfinx
msh, cell_tags, facet_tags = gmshio.read_from_msh(mesh_file, MPI.COMM_WORLD, gdim=2)
# %%

import ufl
from dolfinx.fem import Function, FunctionSpace, dirichletbc, locate_dofs_geometrical, form
from dolfinx.fem.petsc import LinearProblem
import numpy as np

deg = 1
V = FunctionSpace(msh, ("Lagrange", deg))

# Trial and test functions
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
k0 = 4 * np.pi

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx - k0**2 * ufl.inner(u, v) * ufl.dx

# Example: Dirichlet BC on left and bottom walls using coordinates
dofs_D = locate_dofs_geometrical(
    V, lambda x: np.isclose(x[0], -5.0) | np.isclose(x[1], -1.0)
)

u_bc = Function(V)
u_bc.interpolate(lambda x: np.zeros_like(x[0]))  # u=0 for Dirichlet
bcs = [dirichletbc(u_bc, dofs_D)]

# Example using facet tags:
# facets_left = facet_tags.find(1)  # 1 = tag of left wall in GMSH
# bcs = [dirichletbc(u_bc, locate_dofs_topological(V, msh.topology.dim-1, facets_left))]

theta = np.pi / 4
A = 1.0

u_exact = Function(V)
u_exact.interpolate(lambda x: A * np.exp(1j * k0 * (np.cos(theta)*x[0] + np.sin(theta)*x[1])))

n = ufl.FacetNormal(msh)
g = -ufl.dot(n, ufl.grad(u_exact))
L = -ufl.inner(g, v) * ufl.ds

uh = Function(V)
problem = LinearProblem(
    a, L, bcs=bcs, u=uh, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
)
_ = problem.solve()

from dolfinx.fem import assemble_scalar

diff = uh - u_exact

H1_diff = assemble_scalar(form(ufl.inner(ufl.grad(diff), ufl.grad(diff)) * ufl.dx))
H1_exact = assemble_scalar(form(ufl.inner(ufl.grad(u_exact), ufl.grad(u_exact)) * ufl.dx))
print("Relative H1 error:", np.sqrt(H1_diff)/np.sqrt(H1_exact))

L2_diff = assemble_scalar(form(ufl.inner(diff, diff) * ufl.dx))
L2_exact = assemble_scalar(form(ufl.inner(u_exact, u_exact) * ufl.dx))
print("Relative L2 error:", np.sqrt(L2_diff)/np.sqrt(L2_exact))
