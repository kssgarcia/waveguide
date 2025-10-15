# -*- coding: utf-8 -*-
"""
Computes the analytical solution of the wave equation in a square
for homogeneous Dirichlet boundary conditions on the four sides.

The intial condition is a static Gaussian, u_t(t=0) = 0.

Created on Mon Aug 11 14:51:57 2025
@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import matplotlib.colors as mcolors

pi = np.pi
sin = np.sin
cos = np.cos

a = 1
b = 1
h = 0.01
Nx = int(a/h) + 1
Ny = int(b/h) + 1

x = np.linspace(0, a, Nx)
y = np.linspace(0, b, Ny)

X, Y = np.meshgrid(x,y)

Nm = 20
Nn = 20
c = 1 # wave velocity

def f(x,y):
    """Initial position u(0, x, y) = f(x, y)"""
    return np.exp(-100*((x - 0.5)**2 + (y - 0.5)**2))

# In this case, the initial velocity is zero, so u_t(0, x,y) = g(x, y) = 0
# and the coefficients b_mn = 0.

# coefficient matrix, first row and first column won't be used
a_mn = np.zeros((Nm + 1, Nn + 1))

for m in range(1, Nm + 1):
    for n in range(1, Nn + 1):
        # computing coefficient a_mn using Eq. (11.152) in Olver, page 489
        integrand = lambda y,x: f(x,y)*sin(m*pi*x/a)*sin(n*pi*y/b)
        I, err = dblquad(integrand, 0, a, 0, b)
        a_mn[m, n] = 4/(a*b)*I

np.savez("a_mn_dirichlet", a_mn=a_mn)
        
def u_fun(t):
    """ computes analytical solution at time t
    using the Fourier series at the end of page 488 in Olver, truncated
    to Nm values of m and Nn values of n"""

    u = 0

    for m in range(1, Nm + 1):
        for n in range(1, Nn + 1):
            omega_mn = pi*c*np.sqrt(m**2/a**2 + n**2/b**2)
            u += a_mn[m, n]*cos(omega_mn*t)*sin(m*pi*X/a)*sin(n*pi*Y/b)

    return u

plt.figure(1)
plt.clf()
fig, axs = plt.subplots(3, 3, num=1, subplot_kw={"projection": "3d"})
t_list = 0.2*np.arange(9)

norm = mcolors.Normalize(vmin=-0.1, vmax=0.5)

for el, t in enumerate(t_list):
    print(f't = {t:.1f}')
    u = u_fun(t)
    i = el // 3
    j = el % 3
    ax = axs[i, j]
    # ax.imshow(u)
    surf = ax.plot_surface(X, Y, u, cmap='coolwarm', edgecolor='none', norm=norm)
    ax.set(zlim=(-0.1,1))
    ax.set_title(f'$t = {t:.1f}$')
    ax.set_xticks([])
    ax.set_yticks([])
    # fig.colorbar(surf, ax=axs[i,j])

# Single shared colorbar
fig.colorbar(surf, ax=axs.ravel().tolist(), shrink=0.6, aspect=20, label="Value")
plt.savefig('../../../figures/analytical_dirichlet.pdf')
