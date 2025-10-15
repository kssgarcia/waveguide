# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 09:00:49 2025

@author: agarz
"""
import numpy as np
pi = np.pi

epsilon = 1e-1
R0 = 1
b = 1.0
# d = 2*b

Lambda1 = pi**2 / (4 * b**2)
# kb_left, kb_right = pi*0.496, pi*0.499

S = pi * R0**2
mu = R0**2

a = 0.6
alpha = pi * a / b
a0 = (2*b/pi)*np.arctan(S/(2*pi*mu))
sigma = (epsilon**2 * pi**2) / (4 * b**3) * (pi * mu * np.sin(alpha/2)**2 - 0.5 * S * np.cos(alpha/2)**2)
k2_analytic = Lambda1-sigma**2

print('k2_analytic=', k2_analytic)

