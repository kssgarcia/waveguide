# -*- coding: utf-8 -*-
"""
Computes image given to a function phi by DtN operator
according to equation (I.34) of
http://www.cmap.polytechnique.fr/~chesnel/Documents/Waveguides_Invisibility.pdf

Created on Fri Aug 29 17:20:27 2025
@author: agarz
"""
# %%
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

pi = np.pi
sin = np.sin

k = 4

def phi(y): # argument of DtN operator
    return y*(1-y)

fac = np.sqrt(2)
def phi_n(y, n): # basis functions
    return fac*sin(n*pi*y)

y = np.linspace(0,1,101)

N = 100
sum = 0
for n in range(1, N + 1):
    beta_n = np.sqrt(k**2 - n**2*pi**2 + 0*1j)
    c_n, err = quad(lambda y: phi(y)*phi_n(y,n), 0, 1)
    
    sum += 1j*beta_n*c_n*phi_n(y, n)

# image given by DtN operator
Lphi = sum

plt.figure(1)
plt.clf()
plt.plot(y, phi(y), label=r"$\varphi$")
plt.plot(y, np.real(Lphi), label=r"Re[$\Lambda(\varphi)$]")
plt.plot(y, np.imag(Lphi), label=r"Im[$\Lambda(\varphi)$]")
plt.grid(True)
plt.xlabel('$y$')
plt.legend()
