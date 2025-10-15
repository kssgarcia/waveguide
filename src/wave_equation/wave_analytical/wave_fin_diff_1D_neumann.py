# -*- coding: utf-8 -*-
"""
Solves the one-dimensional wave equation with
homogeneous Neumann boundary conditions
by the method of finite differences and compares
the numerical solution with the analytical solution

Created on Thu Aug 14 15:52:32 2025
@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

cos = np.cos
sin = np.sin
pi = np.pi

# import matplotlib
# matplotlib.use('tkagg')
# plt.ion()
# # Display the plot window in non-blocking mode
# plt.show(block=False)


dt = 0.01
dx = 0.01
L = 10
N = int(L/dx) + 1

x = np.linspace(0, L, N)

c = 1
s = 2
x0 = 2
def f(x):
    """u(0,x) = f(x)"""
    return 0.5*(1-np.tanh(s*(x - x0)))

def g(x):
    """u_t(0,x) = g(x)"""
    return 0.5*c*s*(1 - np.tanh(s*(x - x0))**2)

Nn = 20 # number of Fourer series frequencies
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
    """Analytical solution"""
    u = a[0] + b[0]*t
    
    for n in range(1, Nn + 1):
        omega = n*pi*c/L
        u += (a[n]*cos(omega*t) + b[n]*sin(omega*t))*cos(n*pi*x/L)

    return u



sigma = c*dt/dx # stability condition sigma <= 1

B = np.diag(2*(1-sigma**2)*np.ones(N)) + \
    np.diag(sigma**2*np.ones(N-1),1) + \
    np.diag(sigma**2*np.ones(N-1),-1)

B[0,1] = 2*sigma**2
B[-1,-2] = 2*sigma**2

u0 = f(x)
u1 = 0.5*B.dot(u0) + dt*g(x)

nt = int(1e3)
plt.figure(1)
plt.clf()
for i in range(nt):
    u = B.dot(u1) - u0

    if (i+2)%100 == 0:
        plt.clf()
        t = (i+2)*dt
        plt.plot(x,u)
        plt.plot(x,u_fun(t),'--')
        plt.title(rf'$t={t:.3f}$')
        plt.ylim(0,2.1)
        plt.grid(True)
        plt.pause(0.5)
        print(i)
    
    u0 = u1
    u1 = u

