# -*- coding: utf-8 -*-
"""
Computes the analytical solution of the wave equation in a rectangle
for homogeneous Neumann boundary conditions on the four sides.

The initial condition is a plane wave with a front that forms a
given angle theta with the y-axis. When theta = 0, the front propagates
to the right almost without distortion until it gets reflected upwardly at
the right side.

Created on Mon Aug 11 14:51:57 2025
@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import matplotlib.colors as mcolors
import os

pi = np.pi
sin = np.sin
cos = np.cos
tanh = np.tanh

a = 10
b = 10
h = 0.1
Nx = int(a/h) + 1
Ny = int(b/h) + 1

x = np.linspace(0, a, Nx)
y = np.linspace(0, b, Ny)

X, Y = np.meshgrid(x,y)

c = 1 # wave velocity

load_flag = True
if load_flag: # load precomputed coefficient matrices
    data = np.load("a_mn_b_mn_theta_5.npz")
    a_mn = data["a_mn"]
    b_mn = data["b_mn"]
    Nm, Nn = a_mn.shape
    
else: # compute coefficient matrices
    Nm = 21
    Nn = 21

    x0 = 2
    y0 = 2
    s = 2

    # wave front normal unit vector
    theta = 5*pi/180
    nx, ny = cos(theta), sin(theta)

    def f(x,y):
        """Initial position u(0, x, y) = f(x, y)"""
        return 0.5*(1-np.tanh(s*((x - x0)*nx + (y - y0)*ny)))

    def g(x,y):
        """Initial velocity u_t(0, x, y) = g(x, y)"""
        return 0.5*c*s*(1 - np.tanh(s*((x - x0)*nx + (y - y0)*ny))**2)

    # plt.figure(1)
    # plt.clf()
    # plt.imshow(f(X,Y), origin='lower', extent=[0, a, 0, b])
    # plt.xlabel(r'$x$')
    # plt.ylabel(r'$y$')
    # plt.colorbar()

    # Computing coefficient matrix
    a_mn = np.zeros((Nm, Nn))
    b_mn = np.zeros((Nm, Nn))


    norm2 = a*b
    # Note the swapping of x and y  below!!!
    I, err = dblquad(lambda y,x: f(x,y), 0, a, 0, b)
    a_mn[0, 0] = I/norm2

    I, err = dblquad(lambda y,x: g(x,y), 0, a, 0, b)
    b_mn[0, 0] = I/norm2

    print('Computing coefficient matrices')
    for m in range(Nm):
        print(f'm={m} of {Nm}')
        print('n=', end='')
        for n in range(Nn):

            if m > 0 or n > 0:
                print(n, end=',')
                if m == 0 or n == 0:
                    norm2 = 0.5*a*b
                else:
                    norm2 = 0.25*a*b

                integrand = lambda y,x: f(x,y)*cos(m*pi*x/a)*cos(n*pi*y/b)
                I, err = dblquad(integrand, 0, a, 0, b)
                a_mn[m, n] = I/norm2

                omega = pi*c*np.sqrt(m**2/a**2 + n**2/b**2)
                integrand = lambda y,x: g(x,y)*cos(m*pi*x/a)*cos(n*pi*y/b)
                I, err = dblquad(integrand, 0, a, 0, b)
                b_mn[m, n] = I/(omega*norm2)

        print()

    # save coefficient matrices
    np.savez("a_mn_b_mn_theta_5", a_mn=a_mn, b_mn=b_mn)

def u_fun(t):
    """ computes analytical solution at time t"""
    u = a_mn[0, 0] + b_mn[0, 0]*t

    for m in range(Nm):
        for n in range(Nn):
            if m > 0 or n > 0:
                omega = pi*c*np.sqrt(m**2/a**2 + n**2/b**2)
                u += (a_mn[m, n]*cos(omega*t) + b_mn[m, n]*sin(omega*t))*cos(m*pi*X/a)*cos(n*pi*Y/b)
            
    return u

t_list = [0, 2, 4, 6, 7.5, 8, 8.5]
plt.figure(2)
plt.clf()
plt.title('From 2D plane wave, y=0')

for t in t_list:
    print(f't={t}')
    u = u_fun(t)
    plt.plot(x, u[0,:], label=rf'$t={t:.2f}$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    
plt.legend()

plt.figure(3)

t_list = np.arange(0, 8.5, 0.5)
for t in t_list:
    u = u_fun(t)
    plt.clf()
    plt.imshow(u, origin='lower', extent=[0, a, 0, b])
    plt.title(rf'$t={t:.2f}$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.colorbar()
    plt.pause(0.2)

plt.figure(4)
plt.clf()
fig, axs = plt.subplots(3, 3, num=4)
norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

t_list = np.arange(9)
for el, t in enumerate(t_list):
    print(f't = {t:.1f}')
    u = u_fun(t)
    i = el // 3
    j = el % 3
    ax = axs[i, j]
    hmap = ax.imshow(u, origin='lower', extent=[0, a, 0, b], norm=norm)
    ax.set_title(f'$t = {t:.1f}$')
    ax.set_xticks([])
    ax.set_yticks([])

fig.colorbar(hmap, ax=axs.ravel().tolist(), shrink=0.6, aspect=20, label="u")

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","wave_analytical_2D_neumann.pdf")
plt.savefig(file_path, bbox_inches='tight')
