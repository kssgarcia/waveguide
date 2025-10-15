# -*- coding: utf-8 -*-
"""
Solves the two-dimensional wave equation with
homogeneous Neumann boundary conditions
by the method of finite differences and compares
the numerical and analytical solutions.


Created on Fri Aug 15 08:32:28 2025
@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import os

pi = np.pi
sin = np.sin
cos = np.cos
tanh = np.tanh

a = 1
b = 1
h = 0.005
dt = 2.5e-3
Nx = int(a/h) + 1
Ny = int(b/h) + 1

x = np.linspace(0, a, Nx)
y = np.linspace(0, b, Ny)

# X[i, j] = x[i] = i*h
# Y[i, j] = y[j] = j*h
Y, X = np.meshgrid(y, x)

# example function for testing lap_fun
u = sin(2*pi*X/a)*sin(pi*Y/b)

plt.figure(1)
plt.clf()
plt.imshow(u.T, origin='lower', extent=[0, a, 0, b])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')


def lap_fun(u):
    """ computes Laplacian of u by finite differences
    second-order centered 5-point formula
    Assumes homogeneous Dirichlet boundary conditions
    The boundary grid points must not be included in u"""
    center = slice(1,-1)
    lower = slice(0,-2)
    upper = slice(2, None)

    h2 = 1/h**2
    
    lap_u = np.zeros(u.shape)
    lap_u[center, center] = h2*(u[lower, center] +  u[upper, center] + u[center, lower] + u[center, upper] - 4*u[center, center])
    lap_u[0, center] = h2*(u[1, center] + u[0, lower] + u[0, upper] - 4*u[0, center])
    lap_u[-1, center] = h2*(u[-2, center] + u[-1, lower] + u[-1, upper] - 4*u[-1, center])
    lap_u[center, 0] = h2*(u[lower, 0] + u[upper, 0] + u[center, 1] - 4*u[center, 0])
    lap_u[center, -1] = h2*(u[lower, -1] + u[upper, -1] + u[center, -2] - 4*u[center, -1])
    lap_u[0, 0] = h2*(u[1, 0] + u[0, 1] - 4*u[0, 0])
    lap_u[0, -1] = h2*(u[1, -1] + u[0,-2] - 4*u[0, -1])
    lap_u[-1, 0] = h2*(u[-2, 0] + u[-1, 1] - 4*u[-1, 0])
    lap_u[-1, -1] = h2*(u[-2, -1] + u[-1, -2]- 4*u[-1, -1])

    return lap_u

lap_u = np.zeros((Nx, Ny))
lap_u[1:-1, 1:-1] = lap_fun(u[1:-1, 1:-1])
lap_real = -((2*pi/a)**2 + (pi/b)**2)*sin(2*pi*X/a)*sin(pi*Y/b)
dlap = abs(lap_u - lap_real)
print('error=', np.max(np.max(dlap)))

plt.figure(2)
plt.clf()
plt.imshow(lap_u.T, origin='lower', extent=[0, a, 0, b])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

c = 1 # wave velocity
x0 = 2
y0 = 2
s = 2

# load coefficient matrices of analytical solution
data = np.load("a_mn_dirichlet.npz")
a_mn = data['a_mn']
# b_mn = data['b_mn']
Nm, Nn = a_mn.shape

def u_fun(t):
    """ computes analytical solution at time t
    using the Fourier series at the end of page 488 in Olver, truncated
    to Nm values of m and Nn values of n"""

    u = 0

    for m in range(1, Nm):
        for n in range(1, Nn):
            omega_mn = pi*c*np.sqrt(m**2/a**2 + n**2/b**2)
            u += a_mn[m, n]*cos(omega_mn*t)*sin(m*pi*X/a)*sin(n*pi*Y/b)

    return u


def f(x,y):
    """Initial position u(0, x, y) = f(x, y)"""
    return np.exp(-100*((x - 0.5)**2 + (y - 0.5)**2))

def g(x,y):
    """Initial velocity u_t(0, x, y) = g(x, y)"""
    return np.zeros(x.shape)

u0 = f(X,Y)[1:-1, 1:-1] # Excluding borders
u1 = u0 + 0.5*dt**2*c**2*lap_fun(u0) + dt*g(X,Y)[1:-1, 1:-1] # Excluding borders

u_borders = np.zeros((Nx, Ny))
nt = 400
for i in range(nt):
    u = 2*u1 + dt**2*c**2*lap_fun(u1) - u0

    if (i+2)%10 == 0:
        u_borders[1:-1, 1:-1] = u
        plt.figure(1)
        plt.clf()
        t = (i+2)*dt
        plt.plot(x, u_borders[:, int(Ny/2)], label='numerical')
        u_a = u_fun(t)
        plt.plot(x, u_a[:,int(Ny/2)],'--', label='analytical')
        # plt.imshow(u.T, origin='lower')
        # plt.plot(x,u)
        # plt.plot(x,u_fun(t),'--')
        plt.title(rf'$t={t:.3f}, y=0.5$')
        # plt.ylim(-1, 1)
        plt.grid(True)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$u$')
        plt.legend()
        plt.pause(0.1)
        print(f'i={i} of {nt}')
    
    u0 = u1
    u1 = u

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","wave_fin_diff_2D_dirichlet.pdf")
plt.savefig(file_path, bbox_inches='tight')
