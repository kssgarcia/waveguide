# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 11:04:21 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel1, yv

x = np.linspace(0.01,20,500)
H0 = lambda z: hankel1(0,z)
Y0 = lambda z: yv(0,z)

plt.figure(1)
plt.plot(x,np.imag(H0(x)))
plt.plot(x,Y0(x),'--')
plt.grid(True)

def gamma_s(t):
    if t <= 1:
        return -1j*np.sqrt(1-t**2)
    else:
        return np.sqrt(t**2-1)
    
gamma = np.vectorize(gamma_s)
    
plt.figure(2)
t = np.linspace(0,2,500)
plt.plot(t, np.real(gamma(t)), label='Re \gamma')
plt.plot(t, np.imag(gamma(t)), label='Im \gamma')
plt.legend()

k = 1

