# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 10:29:49 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

pi = np.pi

mm = np.arange(-2,3)

beta = np.sqrt(2)
k = 2
d = 1
p = 2*pi/d

if not abs(beta) < k:
    print('abs(beta) should be less than k')
    sys.exit(1)


Nk = int((k - beta)/p)
Mk = int((k + beta)/p)

beta_mm = beta + mm*p

angle = np.linspace(0,2*pi)
plt.figure(1)
plt.clf()
plt.plot(k*np.cos(angle), k*np.sin(angle), 'r')
plt.scatter(np.real(beta_mm), np.imag(beta_mm))
for m, beta_m in zip(mm, beta_mm):
    plt.text(np.real(beta_m), np.imag(beta_m), rf"$\beta_{{{m}}}$")
    print(m, beta_m)

plt.axis('equal')
plt.grid(True)

