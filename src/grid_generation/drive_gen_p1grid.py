# %%
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 12:04:26 2025

@author: agarz
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import src.grid_generation.gen_p1grid as gen

ne, nn, p, conn, gbc = gen.gen_p1grid(nref=3)

# Create a Triangulation object
triang = mtri.Triangulation(p[:,0], p[:,1], conn)

# Plot the mesh
plt.figure()
plt.triplot(triang, 'ko-') # 'ko-' for black circles at nodes and lines
plt.xlabel('X-coordinate')
plt.ylabel('Y-coordinate')
plt.title('2D Mesh Plot')
plt.grid(True)
plt.show()

