# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 14:31:37 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import yv
from scipy.integrate import quad
import os

a = 1
d = 2
k = 0.6

def rho(theta):
    return a

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


x_min = 0
x_max = 3*d
y_min = 0
y_max = d

h = 0.02
nx = int(x_max/h)+1
ny = int(y_max/h) + 1

x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min,y_max, ny)
Y, X = np.meshgrid(y, x)
# Axes convention
# -----> y
# |
# |
# V x
xi, eta = np.sqrt(2)*np.array([1,1])
G0xy = G0(X, Y, xi, eta)


# integral term

plt.figure(1)
plt.clf()
plt.imshow(G0xy.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
plt.colorbar()
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')

plt.figure(2)
plt.clf()

for j in (np.array([0, 1.3, 1.5, 2])/h).astype('int'):
    plt.plot(x, G0xy[:,j], label = f'y={y[j]}')

plt.legend()
plt.xlabel('x')
plt.grid(True)

plt.figure(3)
plt.clf()
for i in (np.array([0, 1.3, 1.5, 6])/h).astype('int'):
    print(i)
    plt.plot(y, G0xy[i,:], label = f'x={x[i]}')

plt.legend()
plt.xlabel('y')
plt.grid(True)

def fint_s(t, x, y, xi, eta):
    if t < 1:
        gamma = -1j*np.sqrt(1-t**2)
    else:
        gamma = np.sqrt(t**2-1)

    numerator = np.exp(-k*gamma*d)*np.cosh(k*gamma*(d-y))*np.cosh(k*gamma*(d-eta))*np.cos(k*x*t)*np.cos(k*xi*t)
    denominator = gamma*np.cosh(k*gamma*d)

    return np.real(numerator/denominator)

fint = np.vectorize(fint_s)
tt = np.linspace(0,10*d+1e-6,1001)
plt.figure(4)
plt.clf()
ftt = fint(tt, 1, 1, xi, eta)
plt.plot(tt,ftt)
plt.grid(True)

def G1_s(x, y, xi, eta):
    I, err = quad(fint_s, 0, 10*d, args=(x, y, xi, eta), points=[1])
    return I

# G1 = np.vectorize(G1_s)
# G1xy = np.zeros((nx,ny))

# # G1xy = G1(X, Y, xi, eta)
# for i in range(nx): 
#     print(f"i={i} of {nx}")
#     for j in range(ny):
#         G1xy[i,j] = G1_s(x[i], y[j], xi, eta)

# np.savez("G1", G1xy=G1xy)
data=np.load("G1.npz")
G1xy = data['G1xy']
plt.figure(5)
plt.clf()
G01xy = G0xy + 2/np.pi*G1xy
plt.imshow(G01xy.T, origin='lower')
plt.colorbar()

plt.figure(6)
plt.clf()

for j in (np.array([0, 1.3, 1.5, 2])/h).astype('int'):
    plt.plot(x, G01xy[:,j], label = f'y={y[j]}')

plt.legend()
plt.xlabel('x')
plt.grid(True)

plt.figure(7)
plt.clf()
for i in (np.array([0, 1.3, 1.5, 6])/h).astype('int'):
    print(i)
    plt.plot(y, G01xy[i,:], label = f'x={x[i]}')

plt.legend()
plt.xlabel('y')
plt.grid(True)

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

def Gn(psi, theta):
    x = a*np.cos(psi)
    y = a*np.sin(psi)
    xi = a*np.cos(theta)
    eta = a*np.sin(theta)
    xi_p = -a*np.sin(theta)
    eta_p = a*np.cos(theta)

    if abs(psi-theta) > 1e-10:
        return 1/a*(xi_p*Gy(x, y, xi, eta) - eta_p*Gx(x, y, xi, eta))
    else:
        Y0_n = 1/(4*np.pi*a**3)*(-a**2)
        G_reg_n = 1/a*(xi_p*Gy_reg(x, y, xi, eta) - eta_p*Gx_reg(x, y, xi, eta))
        return Y0_n + G_reg_n

M = 8
theta = (np.arange(1,M+1)-0.5)*np.pi/(2*M)

kd_list = np.array([1.2, 1.3,1.39,1.391, 1.3913, 1.39131, 1.39132, 1.3914, 1.392, 1.4, 1.5])
# kd_list = np.array([1.39131]) # almost zero determinant
det_list = np.zeros(len(kd_list))
for el, kd in enumerate(kd_list):
    print(f'kd={kd}')
    k = kd/d
    K = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
#            print(f"(i,j)=({i},{j})")
            K[i,j] = Gn(theta[i], theta[j])

    A = np.identity(M) - np.pi/M*a*K
    detA = np.linalg.det(A)
    det_list[el] = detA
    print(detA)
    
np.savez("kd_det_circle", kd_list=kd_list, det_list=det_list)
    
plt.figure(8)
plt.clf()
plt.plot(kd_list, det_list,'o-')
plt.grid(True)

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","linton1992_det_circle.pdf")
plt.savefig(file_path, bbox_inches='tight')


eig_v = np.linalg.eig(A)
eig_val = eig_v.eigenvalues

plt.figure(9)
plt.clf()
plt.scatter(np.real(eig_val), np.imag(eig_val))
plt.grid(True)

eig_vec = eig_v.eigenvectors

plt.figure(10)
plt.clf()
plt.plot(theta/np.pi, np.real(eig_vec[:,0]),'.-')
plt.xlabel(r"$\theta/\pi$")
plt.grid(True)
