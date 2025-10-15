# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 15:34:08 2025

@author: agarz
"""


import numpy as np
import matplotlib.pyplot as plt

cos = np.cos
sin = np.sin
pi = np.pi

c = 1

def f(x,t):
    return sin(t)*sin(pi*x)
 
def u_analitica(t, x):
    wt = -1/((pi*c)**3 -pi*c)*sin(pi*c*t) + 1/((pi*c)**2 - 1)*sin(t)
    return wt*sin(pi*x)
 

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

#fig = plt.figure()

u_a = u_analitica(0, x)
plt.cla()
plt.plot(x,u_full, label='numerica')
plt.plot(x,u_a, '--', label='analitica')
plt.legend()
plt.title('t=0')
plt.ylim(-0.6, 0.6)
plt.xlabel("x")
plt.ylabel("u")
filename = "num_t_nh_0.0.pdf"
plt.savefig(filename)
plt.show()

x_short = x[1:-1]
nt=630
for i in range(nt):
    ts = (i+1)*dt
    u = B.dot(u1) + dt**2*f(x_short, ts) - u0 
    if (i+2) % 120 == 0:
        u_full = np.concatenate(([0],u,[0]))
        t = (i+2)*dt
        u_a = u_analitica(t, x)
        plt.cla()
        plt.plot(x,u_full, label='numerica')
        plt.plot(x,u_a, '--', label='analitica')
        plt.legend()
        plt.title('t='+str(t))
        plt.ylim(-0.2, 0.2)
        plt.xlabel("x")
        plt.ylabel("u")
        filename = "num_t_nh_"+ str(t) + ".pdf"
        plt.savefig(filename)
        plt.pause(0.2)

    u0 = u1
    u1 = u
