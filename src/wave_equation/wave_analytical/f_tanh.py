# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 10:12:04 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import os

x = np.linspace(-4,4,501)
s = 2
f = 0.5*(1-np.tanh(s*x))
plt.plot(x,f)
plt.xlabel(r'$x$')
plt.ylabel(r'$u$')
plt.xlim(-4,4)
plt.grid(True)

waves_path = os.environ["WAVES_PATH"]
print(waves_path)
plt.savefig(os.path.join(waves_path,"figures","f_tanh.pdf"))
