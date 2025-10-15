# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 09:46:39 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice

# values in Table 2 of linton1998greens.pdf
d = 1
k = 2
beta = 0 # np.sqrt(2)
M = 200
Lh = 20 # l = 0, 1, ..., L = 2*lh + 1

S = lattice.lattice_sums(d, k, beta, M, Lh)

xi = 0
eta = 0.8
x = 0.1
# y_vals = [0.0, 0.5*d, d]
y_vals = np.linspace(0.0, d, 101)
n = len(y_vals)
G_vals = np.zeros(n, dtype=np.complex128)

for i, y in enumerate(y_vals):
    X = x - xi
    Y = y - eta
    G_vals[i] = lattice.greens_periodic(X, Y, S, k, d)

plt.figure(1)
plt.clf()
plt.plot(y_vals, np.real(G_vals), '-')
# plt.plot(y_vals, np.imag(G_vals), '-')
plt.xlabel("$y$")
plt.ylabel("Re($G$)")
plt.grid(True)

#------------------------------------------------------------
# Green's function that satisfies Neumann boundary conditions
# on a waveguide

S_2d = lattice.lattice_sums(2*d, k, beta, M, Lh)

# def greens_neumann(x, y, xi, eta, S_2d, k, d):
#     """
#     Computes Green's function that satisfies
#     Neumann boundary conditions
#     """
#     X = x - xi
#     term1 = lattice.greens_periodic(X, y - eta, S_2d, k, 2*d)
#     term2 = lattice.greens_periodic(X, y + eta, S_2d, k, 2*d)
#     return term1 + term2

G_neumann = np.zeros(n, dtype=np.complex128)

for i, y in enumerate(y_vals):
    G_neumann[i] = lattice.greens_neumann(x, y, xi, eta, S_2d, k, d)

plt.figure(2)
plt.clf()
plt.plot(y_vals, np.real(G_neumann), '-')
# plt.plot(y_vals, np.imag(G_vals), '-')
plt.xlabel("$y$")
plt.ylabel("Re($G$)")
plt.grid(True)

#------------------------------------------------------------
# Green's function that satisfies Dirichlet boundary conditions
# on a waveguide
G_dirichlet = np.zeros(n, dtype=np.complex128)

for i, y in enumerate(y_vals):
    G_dirichlet[i] = lattice.greens_dirichlet(x, y, xi, eta, S_2d, k, d)

plt.figure(3)
plt.clf()
plt.plot(y_vals, np.real(G_dirichlet), '-')
# plt.plot(y_vals, np.imag(G_dirichlet), '-')
plt.xlabel("$y$")
plt.ylabel("Re($G$)")
plt.grid(True)


# ------------------------------------------------------------
# Green's function that satisfies mixed boundary conditions
# G = 0 on y = 0
# dG/dn = 0 on y = d

beta = np.pi/(2*d)
S_dn = lattice.lattice_sums(2*d, k, beta, M, Lh)
    
G_dn = np.zeros(n, dtype=np.complex128)

for i, y in enumerate(y_vals):
    G_dn[i] = lattice.greens_dir_neu(x, y, xi, eta, S_dn, k, d)

plt.figure(4)
plt.clf()
plt.plot(y_vals, np.real(G_dn), '-')
# plt.plot(y_vals, np.imag(G_vals), '-')
plt.xlabel("$y$")
plt.ylabel("Re($G$)")
plt.grid(True)
