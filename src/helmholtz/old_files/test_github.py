"""
This code use the package from: https://github.com/Vinc0110/Helmi_FEM/blob/master/examples/waveguide/waveguide.py
It solves the Helmholtz equation use eigenmode boundary treatment (Robin condition)
"""
# %%
import numpy as np
import skfem
from helmi import Helmholtz

x_pts = np.linspace(0, 100, 101)
y_pts = np.linspace(-5, 5, 21)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts)

# Define circular obstacle subdomain
R = 1.0
center = np.array([50.0, 0.0])
mesh = mesh.with_subdomains({
    'air': lambda x: x[0] < 50,
    'plastic': lambda x: x[0] >= 50,
    'obstacle': lambda x: np.sum((x - center[:, None])**2, axis=0) < R**2
})

# Boundaries
mesh = mesh.with_boundaries({
    'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
    'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
    'bound_ymin': lambda x: np.isclose(x[1], y_pts[0]),
    'bound_ymax': lambda x: np.isclose(x[1], y_pts[-1])
})

# FEM setup
element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

k0 = 0.5
eps_air = 1
mu_air = 1
eps_plastic = 2 - 0.1j
mu_plastic = 1
eps_obstacle = 1e6  # Very high to mimic rigid obstacle
mu_obstacle = 1

fem.assemble_subdomains(alpha={
    'air': 1/mu_air,
    'plastic': 1/mu_plastic,
    'obstacle': 1/mu_obstacle
}, beta={
    'air': -k0**2 * eps_air,
    'plastic': -k0**2 * eps_plastic,
    'obstacle': -k0**2 * eps_obstacle
}, f={
    'air': 0,
    'plastic': 0,
    'obstacle': 0
})

# Dirichlet BC for top and bottom, plus obstacle
fem.assemble_boundaries_dirichlet(value={
    'bound_ymin': 0,
    'bound_ymax': 0,
    # Treat obstacle as Dirichlet (rigid)
    'obstacle': 0
})

fem.assemble_boundaries_3rd(gamma={
    'bound_xmin': 1/mu_air * 1j*k0,
    'bound_xmax': 1/mu_plastic * 1j*k0
}, q={
    'bound_xmin': 1/mu_air * 2j*k0,
    'bound_xmax': 0
})

fem.solve()

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as mplt

fig, ax = mplt.subplots(2, 1,figsize=(10, 4))
plot(fem.basis, fem.phi_re, ax=ax[0])
plot(fem.basis, fem.phi_im, ax=ax[1])
ax[0].set_aspect(1)
ax[1].set_aspect(1)
ax[0].set_title('Real Part')
ax[1].set_title('Imaginary Part')
mplt.tight_layout()
mplt.savefig('../../solutions/waveguide_github.png')
mplt.close()