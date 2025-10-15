# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 10:24:55 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
sin = np.sin

x = np.linspace(0,1,501)


nf = 100
u = 0
for n in range(1, nf+1):
    b_n = 2/(n*pi)
    u = u + b_n*sin(n*pi*x)
    
plt.cla()
plt.plot(x,1-x)
plt.plot(x,u)
plt.savefig("fourier_nh_bc.pdf")
