# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 14:31:37 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import yv
from scipy.integrate import quad

aod = 0.5
d = 2
a = aod*d
boa = 1.5
b = boa*a

pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt

# The contours is a circle
def rho(t):
    return a*b/sqrt((a*cos(t))**2 + (b*sin(t))**2)

def rho_p(t):
    return a*b*(a**2 - b**2)*cos(t)*sin(t)/(a**2*cos(t)**2 + b**2*sin(t)**2)**1.5

def rho_pp(t):
    return (a**2*cos(t)**4 + 2*a**2*cos(t)**2*sin(t)**2 - 2*b**2*cos(t)**2*sin(t)**2 - b**2*sin(t)**4)*(a**2 - b**2)*a*b/(a**2*cos(t)**2 + b**2*sin(t)**2)**2.5

Y0 = lambda z: yv(0,z)

def G0(x, y, xi, eta):
    r =  np.sqrt((x - xi)**2 + (y - eta)**2)
    r1 =  np.sqrt((x - xi)**2 + (y + eta - 2*d)**2)
    r2 =  np.sqrt((x + xi)**2 + (y - eta)**2)
    r3 =  np.sqrt((x + xi)**2 + (y + eta - 2*d)**2)
    val = 0.25*(Y0(k*r) + Y0(k*r1) + Y0(k*r2) + Y0(k*r3))
    return val

def G0_reg(x, y, xi, eta):
    r1 =  np.sqrt((x - xi)**2 + (y + eta - 2*d)**2)
    r2 =  np.sqrt((x + xi)**2 + (y - eta)**2)
    r3 =  np.sqrt((x + xi)**2 + (y + eta - 2*d)**2)
    val = 0.25*(Y0(k*r1) + Y0(k*r2) + Y0(k*r3))
    return val


def fint_s(t, x, y, xi, eta):
    if t < 1:
        gamma = -1j*np.sqrt(1-t**2)
    else:
        gamma = np.sqrt(t**2-1)

    numerator = np.exp(-k*gamma*d)*np.cosh(k*gamma*(d-y))*np.cosh(k*gamma*(d-eta))*np.cos(k*x*t)*np.cos(k*xi*t)
    denominator = gamma*np.cosh(k*gamma*d)

    return np.real(numerator/denominator)

fint = np.vectorize(fint_s)

def G1_s(x, y, xi, eta):
    I, err = quad(fint_s, 0, 100*d, args=(x, y, xi, eta), points=[1])
    return I


def G(x, y, xi, eta):
    return G0(x, y, xi, eta) + 2/np.pi*G1_s(x, y, xi, eta)

def G_reg(x, y, xi, eta):
    return G0_reg(x, y, xi, eta) + 2/np.pi*G1_s(x, y, xi, eta)

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

M = 8
theta = (np.arange(1,M+1)-0.5)*np.pi/(2*M)

# kd_list = np.array([1, 1.2, 1.3, 1.4, 1.5])
# det_list = np.zeros(len(kd_list))
# for el, kd in enumerate(kd_list):
#     print(f'kd={kd}')
#     k = kd/d
#     Kw = np.zeros((M,M))
#     for i in range(M):
#         for j in range(M):
# #            print(f"(i,j)=({i},{j})")
#             Kw[i,j] = Gn_w(theta[i], theta[j])

#     A = np.identity(M) - np.pi/M*Kw
#     detA = np.linalg.det(A)
#     det_list[el] = detA
#     print(detA)
    
# plt.figure(7)
# plt.clf()
# plt.plot(kd_list, det_list,'o-')
# plt.plot([kd_list[0],kd_list[-1]],[0,0],'-r')
# plt.grid(True)
    
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

plt.figure(1)
plt.clf()
tt = np.linspace(0,pi/2,500)
dt = tt[1]-tt[0]
#plt.plot(tt,rho(tt))
plt.plot(tt,rho_p(tt))
plt.plot(tt,np.gradient(rho(tt), dt), '--')
plt.plot(tt,rho_pp(tt))
plt.plot(tt,np.gradient(rho_p(tt), dt), '--')
plt.grid(True)

r_theta = rho(tt)
r_prime = rho_p(tt)
r_2prime = rho_pp(tt)

x = r_theta*cos(tt)
y = r_theta*sin(tt)
x_p = np.gradient(x, dt, edge_order=2)
x_pp = np.gradient(x_p, dt, edge_order=2)
y_p = np.gradient(y, dt, edge_order=2)
y_pp = np.gradient(y_p, dt, edge_order=2)

plt.figure(2)
plt.clf()
plt.plot(tt, r_theta*r_2prime - r_theta**2 - 2*r_prime**2, label='right')
plt.plot(tt, x_pp*y_p - y_pp*x_p, '--', label='fin diff')
plt.plot(tt, r_prime*r_2prime - r_theta**2 - 2*r_prime**2, label='wrong')
plt.grid(True)
plt.legend()
