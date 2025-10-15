# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 12:26:17 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# data = np.load("circle_mesh.npz")
data = np.load("square_mesh.npz")
nodes = data['nodes']
elements = data['elements']

triang = tri.Triangulation(nodes[:,0], nodes[:,1], elements)

plt.figure(1)
plt.clf()
plt.triplot(triang, color='k')
