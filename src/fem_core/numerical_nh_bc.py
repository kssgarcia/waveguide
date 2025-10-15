# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 11:53:33 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt

cos = np.cos
sin = np.sin
pi = np.pi

c = 1

def u_analitica(t, x):
    u = sin(t)*(1-x)
    nf = 100
    for n in range(1, nf+1):
        npc = n*pi*c
        fac = npc**2 - 1
        w_n = 2/(n*pi)*((-1/npc - 1/(npc*fac))*sin(npc*t) + 1/fac*sin(t))
        u = u + w_n*sin(n*pi*x)
        
    return u
 

N = 100
L = 1
dx = L/N
dt = 0.01

sigma = c*dt/dx # stability condition sigma <= 1

B = np.diag(2*(1-sigma**2)*np.ones(N-1)) + \
    np.diag(sigma**2*np.ones(N-2),1) + \
    np.diag(sigma**2*np.ones(N-2),-1)


x = dx*np.arange(N+1)

u_full = np.zeros(N+1)

plt.plot(x,u_full)
plt.xlabel("x")
plt.ylabel("u")
plt.ylim(-0.6,0.6)

u0 = u_full[1:-1]
u1 = 0.5*B.dot(u0)

fig = plt.figure(1)

u_a = u_analitica(0, x)
plt.cla()
plt.plot(x,u_full, label='numerica')
plt.plot(x,u_a, '--', label='analitica')
plt.legend()
plt.title('t=0')
plt.ylim(-0.6, 0.6)
plt.xlabel("x")
plt.ylabel("u")
filename = "num_nh_bc_t_0.0.pdf"
plt.savefig(filename)
plt.show()

def g(t):
    return sin(t)

nt=600
for i in range(nt):
    u = B.dot(u1) - u0
    ts = (i+1)*dt
    u[0] += sigma**2*g(ts) # nonhomogeneous boundary condition
    if (i+2) % 120 == 0:
        t = (i+2)*dt
        u_full = np.concatenate(([g(t)],u,[0]))
        u_a = u_analitica(t, x)
        plt.cla()
        plt.plot(x,u_full, label='numerica')
        plt.plot(x,u_a, '--', label='analitica')
        plt.legend()
        plt.title('t='+str(t))
        plt.ylim(-1.2, 1.2)
        plt.xlabel("x")
        plt.ylabel("u")
        filename = "num_nh_bc_t_"+ str(t) + ".pdf"
        plt.savefig(filename)
        plt.pause(0.2)

    u0 = u1
    u1 = u
