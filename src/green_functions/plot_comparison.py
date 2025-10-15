# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 11:39:38 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import os

N = 5
epsilon = np.zeros(N)
s_root = np.zeros(N)
sigma_analytic = np.zeros(N)
sigma_numeric = np.zeros(N)

for i in range(N):
    filename = f'circle_k2_{i+1}.npz'
    print(filename)
    data = np.load(filename)
    epsilon[i] = data['epsilon']
    s_root[i] = data['s_root']
    sigma_analytic[i] = data['sigma_analytic']
    sigma_numeric[i] = data['sigma_numeric']

ix_sorted = np.argsort(epsilon)
epsilon = epsilon[ix_sorted]
s_root = s_root[ix_sorted]
sigma_analytic = sigma_analytic[ix_sorted]
sigma_numeric = sigma_numeric[ix_sorted]

plt.figure(1)
plt.clf()
plt.loglog(epsilon, sigma_numeric, 'o-', lw=2,
           label='$\sigma_{num}$ numeric')
plt.loglog(epsilon, sigma_analytic, '.--', lw=2,
           label='$\sigma$ analytical')
# plt.grid(True)
plt.tick_params(axis='both', labelsize=12)
plt.xlabel(r'$\varepsilon$', fontsize=18)
plt.ylabel(r'$\sigma$', fontsize=18)
plt.legend(fontsize=14)

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"Portugal 2025. AGL","Portugal 2025. AGL",
                         "sigma_analytic_numeric.pdf")
plt.savefig(file_path, bbox_inches='tight')

plt.show()

error_order = np.abs(epsilon**3*np.log(epsilon))
del_sigma = np.abs(sigma_analytic - sigma_numeric)
plt.figure(2)
plt.clf()
str_nerror = r"$|\sigma - \sigma_{num}|$"
str_aerror = r"$|\varepsilon^3\log(\varepsilon)|$"
plt.loglog(epsilon, del_sigma, 'o-', lw=2, label=str_nerror)
plt.loglog(epsilon, error_order, 'o-', lw=2, label=str_aerror)
plt.tick_params(axis='both', labelsize=12)
plt.xlabel(r'$\varepsilon$', fontsize=18)
plt.ylabel('Difference', fontsize=18)
plt.legend(fontsize=14)

file_path = os.path.join(waves_path,"Portugal 2025. AGL","Portugal 2025. AGL",
                         "difference_sigma.pdf")
plt.savefig(file_path, bbox_inches='tight')

plt.show()



