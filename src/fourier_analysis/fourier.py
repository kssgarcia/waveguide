# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 21:49:04 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
cos = np.cos
sin = np.sin
pi = np.pi

def u_analitica(t, x):
    u = 0
    N = 100
    for k in range(N):
        n = 2*k + 1
        b_n = 4*(-1)**k/(pi**2*n**2)
        u = u + b_n*cos(n*pi*t)*sin(n*pi*x)
    
    return u

x = np.linspace(0,1,501)

tt = np.array([0, 0.25, 0.5, 0.75, 1, 1.25])
labels = ["0", "1/4","1/2","3/4","1","5/4"]

fig = plt.figure()
for t,lab in zip(tt,labels):    

    u = u_analitica(t, x)
    plt.cla()
    label="$t="+lab+"$"
    plt.plot(x,u, label=label)
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("u")
    plt.ylim(-0.6,0.6)
    filename = "t_"+str(t)+".pdf"
    plt.savefig(filename)
    plt.show()
