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

a = 10
b = 10
h = 0.05
dt = 0.025
Nx = int(a/h) + 1
Ny = int(b/h) + 1

x = np.linspace(0, a, Nx)
y = np.linspace(0, b, Ny)

# X[i, j] = x[i] = i*h
# Y[i, j] = y[j] = j*h
Y, X = np.meshgrid(y, x)

u = cos(2*pi*X/a)*cos(pi*Y/b)

plt.figure(1)
plt.clf()
plt.imshow(u.T, origin='lower', extent=[0, a, 0, b])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')


def lap_fun(u):
    """ computes Laplacian of u by finite differences
    second-order centered 5-point formula"""
    center = slice(1,-1)
    lower = slice(0,-2)
    upper = slice(2, None)

    h2 = 1/h**2
    
    lap_u = np.zeros((Nx, Ny))
    lap_u[center, center] = h2*(u[lower, center] +  u[upper, center] + u[center, lower] + u[center, upper] - 4*u[center, center])
    lap_u[0, center] = h2*(2*u[1, center] + u[0, lower] + u[0, upper] - 4*u[0, center])
    lap_u[-1, center] = h2*(2*u[-2, center] + u[-1, lower] + u[-1, upper] - 4*u[-1, center])
    lap_u[center, 0] = h2*(u[lower, 0] + u[upper, 0] + 2*u[center, 1] - 4*u[center, 0])
    lap_u[center, -1] = h2*(u[lower, -1] + u[upper, -1] + 2*u[center, -2] - 4*u[center, -1])
    lap_u[0, 0] = h2*(2*u[1, 0] + 2*u[0, 1] - 4*u[0, 0])
    lap_u[0, -1] = h2*(2*u[1, -1] + 2*u[0,-2] - 4*u[0, -1])
    lap_u[-1, 0] = h2*(2*u[-2, 0] + 2*u[-1, 1] - 4*u[-1, 0])
    lap_u[-1, -1] = h2*(2*u[-2, -1] + 2*u[-1, -2]- 4*u[-1, -1])

    return lap_u

lap_u = lap_fun(u)
lap_real = -((2*pi/a)**2 + (pi/b)**2)*cos(2*pi*X/a)*cos(pi*Y/b)
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
data = np.load("a_mn_b_mn_theta_5.npz")
a_mn = data['a_mn']
b_mn = data['b_mn']
Nm, Nn = a_mn.shape

def u_fun(t):
    """ computes analytical solution at time t"""
    u = a_mn[0, 0] + b_mn[0, 0]*t

    for m in range(Nm):
        for n in range(Nn):
            if m > 0 or n > 0:
                omega = pi*c*np.sqrt(m**2/a**2 + n**2/b**2)
                u += (a_mn[m, n]*cos(omega*t) + b_mn[m, n]*sin(omega*t))*cos(m*pi*X/a)*cos(n*pi*Y/b)
            
    return u


# wave front normal unit vector
theta = 5*pi/180
nx, ny = cos(theta), sin(theta)

def f(x,y):
    """Initial position u(0, x, y) = f(x, y)"""
    return 0.5*(1-np.tanh(s*((x - x0)*nx + (y - y0)*ny)))

def g(x,y):
    """Initial velocity u_t(0, x, y) = g(x, y)"""
    return 0.5*c*s*(1 - tanh(s*((x - x0)*nx + (y - y0)*ny))**2)

u0 = f(X,Y)
u1 = u0 + 0.5*dt**2*c**2*lap_fun(u0) + dt*g(X,Y)

nt = 280
for i in range(nt):
    u = 2*u1 + dt**2*c**2*lap_fun(u1) - u0

    if (i+2)%10 == 0:
        plt.figure(1)
        plt.clf()
        t = (i+2)*dt
        plt.plot(x, u[:, 0], label='numerical')
        u_a = u_fun(t)
        plt.plot(x, u_a[:,0],'--', label='analytical')
        # plt.imshow(u.T, origin='lower')
        # plt.plot(x,u)
        # plt.plot(x,u_fun(t),'--')
        plt.title(rf'$t={t:.3f}, y=0$')
        plt.xlim(0,10)
        plt.ylim(-0.1,1.1)
        plt.grid(True)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$u$')
        plt.legend()
        plt.pause(0.1)
        print(f'i={i} of {nt}')
    
    u0 = u1
    u1 = u

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","wave_fin_diff_2D_neumann.pdf")
plt.savefig(file_path, bbox_inches='tight')
