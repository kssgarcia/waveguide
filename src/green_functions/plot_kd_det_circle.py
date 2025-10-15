# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 09:33:00 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import os

data = np.load("kd_det_circle.npz")
kd_list = data["kd_list"]
det_list = data["det_list"]

plt.figure(1)
plt.clf()
plt.plot(kd_list, det_list,'o-')
plt.grid(True)
plt.xlabel(r'$kd$')
plt.ylabel(r'$\det A$')

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"figures","kd_det_circle.pdf")
plt.savefig(file_path, bbox_inches='tight')

plt.figure(2)
plt.clf()
# range=slice(4,8)
plt.plot(kd_list, det_list,'o-')
plt.xlim(1.3909, 1.3915)
plt.ylim(-0.001, 0.002)
plt.grid(True)
plt.xlabel(r'$kd$')
plt.ylabel(r'$\det A$')

file_path = os.path.join(waves_path,"figures","kd_det_circle_detail.pdf")
plt.savefig(file_path, bbox_inches='tight')
