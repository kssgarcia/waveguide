# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 07:21:56 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt

cos = np.cos
sin = np.sin
pi = np.pi

c = 1

def u_analitica(t, x):
    u = 0
    nf = 100
    for k in range(nf):
        n = 2*k + 1
        b_n = 4*(-1)**k/(pi**2*n**2)
        u = u + b_n*cos(n*pi*t)*sin(n*pi*x)
        
    return u
 

N = 100
L = 1
dx = L/N
dt = 0.01

sigma = c*dt/dx # stability condition sigma <= 1

B = np.diag(2*(1-sigma**2)*np.ones(N-1)) + \
    np.diag(sigma**2*np.ones(N-2),1) + \
    np.diag(sigma**2*np.ones(N-2),-1)

def u0_fun(x):
    if x <= 0.5:
        return x
    else:
        return 1 - x

u0_vec = np.vectorize(u0_fun)
x = dx*np.arange(N+1)
u_full = u0_vec(x)

plt.plot(x,u_full)
plt.xlabel("x")
plt.ylabel("u")
plt.ylim(-0.6,0.6)

u0 = u_full[1:-1]
u1 = 0.5*B.dot(u0)

fig = plt.figure()

u_a = u_analitica(0, x)
plt.cla()
plt.plot(x,u_full, label='numerica')
plt.plot(x,u_a, '--', label='analitica')
plt.legend()
plt.title('t=0')
plt.ylim(-0.6, 0.6)
plt.xlabel("x")
plt.ylabel("u")
filename = "num_t_0.0.pdf"
plt.savefig(filename)
plt.show()


nt=125
for i in range(nt):
    u = B.dot(u1) - u0
    if (i+2) % 25 == 0:
        u_full = np.concatenate(([0],u,[0]))
        t = (i+2)*dt
        u_a = u_analitica(t, x)
        plt.cla()
        plt.plot(x,u_full, label='numerica')
        plt.plot(x,u_a, '--', label='analitica')
        plt.legend()
        plt.title('t='+str(t))
        plt.ylim(-0.6, 0.6)
        plt.xlabel("x")
        plt.ylabel("u")
        filename = "num_t_"+ str(t) + ".pdf"
        plt.savefig(filename)
        plt.show()

    u0 = u1
    u1 = u
