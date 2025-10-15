# -*- coding: utf-8 -*-
"""
Computes the analytical solution of the 1D wave equation
with homogeneous Neumann boundary conditions at both ends.

The initial condition is a wave front with a
hyperbolic tangent profile. The front propagates to the
right almost without distortion until it gets refleted upwardly at
right boundary.

Created on Wed Aug 13 14:34:55 2025
@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os

cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

L = 10
h = 0.1
Nx = int(L/h) + 1

x = np.linspace(0, L, Nx)

c = 1
s = 2
x0 = 2
def f(x):
    return 0.5*(1-np.tanh(s*(x - x0)))

def g(x):
    return 0.5*c*s*(1 - tanh(s*(x - x0))**2)


Nn = 20
a = np.zeros(Nn + 1)
b = np.zeros(Nn + 1)

I, err = quad(f, 0, L)
a[0] = 1/L*I

I, err = quad(g, 0, L)
b[0] = 1/L*I

for n in range(1, Nn + 1):
    I, err = quad(lambda x: f(x)*cos(n*pi*x/L), 0, L)
    a[n] = 2/L*I

    I, err = quad(lambda x: g(x)*cos(n*pi*x/L), 0, L)
    b[n] = 2/(n*pi*c)*I

def u_fun(t):
    u = a[0] + b[0]*t
    
    for n in range(1, Nn + 1):
        omega = n*pi*c/L
        u += (a[n]*cos(omega*t) + b[n]*sin(omega*t))*cos(n*pi*x/L)

    return u

plt.figure(2)
plt.clf()
# plt.plot(x,f(x))
plt.grid(True)

t_list = [0, 2, 4, 6, 7, 7.5, 7.8, 8, 8.2, 9, 10]

for t in t_list:
    if t < 10:
        line_style = "-"
    else:
        line_style = "--"
        
    u = u_fun(t)
    plt.plot(x,u, line_style, label=rf'$t={t:.1f}$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')

plt.xlim(0,L)
plt.title('1D solution')
plt.legend(bbox_to_anchor=(1,1))

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","wave_analytical_1D_neumann.pdf")
plt.savefig(file_path, bbox_inches='tight')

