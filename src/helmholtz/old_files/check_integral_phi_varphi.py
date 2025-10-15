# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 10:02:58 2025

@author: agarz
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

b = 0.5
M = 8
h = 1.0/M

def varphi(y, n):
    lam = n*np.pi/(2*b)
    return 1/np.sqrt(b)*np.sin(lam*(y + b))

y = np.linspace(-b,b, 501)
plt.figure(1)
plt.clf()
plt.plot(y, varphi(y, 2))
plt.grid(True)

y_i = np.arange(-b, b+h/2, h)

def phi(y, i):
    
    if i == 0:
        if y < y_i[1]:
            return 1/h*(y_i[1] - y)
        else:
            return 0
    elif i == M:
        if y > y_i[M-1]:
            return 1/h*(y - y_i[M-1])
        else:
            return 0
    else:
        if y < y_i[i-1]:
            return 0
        elif y < y_i[i]:
            return 1/h*(y - y_i[i-1])
        elif y < y_i[i+1]:
            return 1/h*(y_i[i+1] - y)
        else:
            return 0
        
            


phi_y = np.zeros(y.shape[0])
for i, yp in enumerate(y):
    phi_y[i] = phi(yp, 1)

plt.figure(2)
plt.clf()
plt.plot(y, phi_y,'-')
plt.grid(True)

for i in range(M+1):
    I, err = quad(lambda y: phi(y, i)*varphi(y, 2), -b, b)
    print(I)
    
I, err = quad(lambda y: varphi(y, 2)*varphi(y, 2), -b, b)