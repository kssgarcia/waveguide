# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 14:31:37 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice
import sys

a = 0.2
b = 1 # -b <= y <= b
h = 0.7 # height of circle center over midline

d = 2*b

pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt

# The contours is a circle of radius a
def rho(t):
    """
    Computes radius as function of angle t for
    a circle of radius a
    """
    return a
    
def rho_p(t):
    return 0

def rho_pp(t):
    return 0
    
epsilon = 1e-6
def Gx(x, y, xi, eta):
    Gm = G(x, y, xi - epsilon, eta)
    Gp = G(x, y, xi + epsilon, eta)
    return (Gp - Gm)/(2*epsilon)

def Gy(x, y, xi, eta):
    Gm = G(x, y, xi, eta - epsilon)
    Gp = G(x, y, xi, eta + epsilon)
    return (Gp - Gm)/(2*epsilon)

def Gx_reg(x, y, xi, eta):
    Gm = G_reg(x, y, xi - epsilon, eta)
    Gp = G_reg(x, y, xi + epsilon, eta)
    return (Gp - Gm)/(2*epsilon)

def Gy_reg(x, y, xi, eta):
    Gm = G_reg(x, y, xi, eta - epsilon)
    Gp = G_reg(x, y, xi, eta + epsilon)
    return (Gp - Gm)/(2*epsilon)

def Gn_w(psi, theta):
    r_psi = rho(psi)
    x = r_psi*cos(psi)
    y = r_psi*sin(psi)
    
    r_theta = rho(theta)
    xi = r_theta*cos(theta)
    eta = r_theta*sin(theta)
    
    r_prime = rho_p(theta)
    xi_p = r_prime*cos(theta) - r_theta*sin(theta)
    eta_p = r_prime*sin(theta) + r_theta*cos(theta)

    w = sqrt(r_theta**2 + r_prime**2)
    
    if abs(psi-theta) > 1e-10:
        return xi_p*Gy(x, y, xi, eta) - eta_p*Gx(x, y, xi, eta)
    else:
        r_2prime = rho_pp(theta)
        Y0_n_w = 1/(4*pi*w**2)*(r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
        G_reg_n_w = xi_p*Gy_reg(x, y, xi, eta) - eta_p*Gx_reg(x, y, xi, eta)
        return Y0_n_w + G_reg_n_w


M = 32
# integration range 0 < theta < pi
# theta = (np.arange(1,M+1)-0.5)*2*np.pi/M
theta = (np.arange(1,M+1)-0.5)*2*np.pi/M

#kb_list = pi*np.array([0.496, 0.497, 0.498, 0.4987, 0.499]) # for a=0.2, b=0.6
kb_list = pi*np.array([0.496, 0.497, 0.498, 0.4987, 0.499])
# kd_list = np.array([1.3, 1.4, 1.45, 1.47, 1.5])
# kd_list = np.array([1.39131]) # almost zero determinant
det_list = np.zeros(len(kb_list), dtype=np.complex128)
for el, kb in enumerate(kb_list):
    print(f'kb={kb}')
    k = kb/b
    S_2d = lattice.lattice_sums(2*d, k, beta=0, M=200, Lh=20)
    G = lambda x, y, xi, eta: lattice.greens_dirichlet(x, y + b + h, xi, eta + b + h, S_2d, k, d)
    G_reg = lambda x, y, xi, eta: lattice.greens_dirichlet_reg(x, y + b + h, xi, eta + b + h, S_2d, k, d)
    Kw = np.zeros((M,M), dtype=np.complex128)
    for i in range(M):
        for j in range(M):
#            print(f"(i,j)=({i},{j})")
            Kw[i,j] = Gn_w(theta[i], theta[j])

    A = np.identity(M) - 4*np.pi/M*Kw
    detA = np.linalg.det(A)
    det_list[el] = detA
    print(detA)
    
# plt.figure(7)
# plt.clf()
# plt.plot(kb_list, det_list,'o-')
# plt.plot([kb_list[0],kb_list[-1]],[0,0],'-r')
# plt.grid(True)

plt.figure(7)
plt.clf()
plt.plot(np.real(det_list), np.imag(det_list), '.-')
# adding labels
for det, kb in zip(det_list, kb_list):
    plt.text(np.real(det), np.imag(det), f'$kb/\pi={kb/pi}$')

plt.plot([0],[0],'or') # red dot at origin
plt.grid(True)
plt.xlabel('Re[det $A$]')
plt.ylabel('Im[det $A$]')

plt.figure(8)
plt.clf()
plt.plot(kb_list/pi, np.real(det_list),'o-')
plt.plot([kb_list[0]/pi,kb_list[-1]/pi],[0,0],'-r')
plt.grid(True)
plt.xlabel('$kd/\pi$')
plt.ylabel('Re[det $A$]')


# eig_v = np.linalg.eig(A)
# eig_val = eig_v.eigenvalues

# plt.figure(8)
# plt.clf()
# plt.scatter(np.real(eig_val), np.imag(eig_val))
# plt.grid(True)

# eig_vec = eig_v.eigenvectors

# plt.figure(9)
# plt.clf()
# plt.plot(np.real(eig_vec[:,0]),'.-')
# plt.grid(True)
