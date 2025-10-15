# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 14:42:12 2025

@author: agarz
"""

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
cos = np.cos
sin = np.sin

import numpy as np

def polygon_fourier_coeffs(vertices, K, period=2*np.pi):
    """
    Compute exact Fourier coefficients of a closed polygonal curve.
    
    Parameters
    ----------
    vertices : array_like of shape (N,2) or complex
        The polygon vertices, ordered and closed. If not closed,
        the function will close it automatically.
        Can be given as Nx2 real array or as 1D complex array.
    K : int
        Maximum Fourier mode to compute. Coefficients for k=-K..K are returned.
    period : float, optional
        Period of the parameterization. Default = 2*pi.
        The polygon is parameterized proportionally to edge length
        (arc-length parameterization).
    
    Returns
    -------
    coeffs : dict
        Dictionary mapping integer k -> complex Fourier coefficient c_k.
        The Fourier series is
            z(t) = sum_{k=-K}^K c_k * exp(i*k*2*pi*t/period),  t in [0,period).
    """
    #verts = np.asarray(vertices, dtype=complex).reshape(-1)
    verts = [x + 1j*y for x, y in vertices]
    if verts[0] != verts[-1]:
        verts = np.append(verts, verts[0])  # close polygon
    
    # edge lengths and cumulative arc-length parameterization
    edges = verts[1:] - verts[:-1]
    lens = np.abs(edges)
    L = np.sum(lens)
    # cumulative parameters scaled to [0,period)
    cumlen = np.concatenate(([0], np.cumsum(lens)))
    tvals = cumlen * (period / L)
    
    coeffs = {}
    for k in range(-K, K+1):
        if k == 0:
            # c0 is just the average over one period
            # integral z(t) dt over [0,period] divided by period
            area = 0.0j
            for n in range(len(edges)):
                z0, z1 = verts[n], verts[n+1]
                dt = tvals[n+1] - tvals[n]
                # integral of linear segment = average value * dt
                area += 0.5 * (z0+z1) * dt
            coeffs[k] = area / period
        else:
            omega = 2*np.pi*k/period
            total = 0.0j
            for n in range(len(edges)):
                z0, z1 = verts[n], verts[n+1]
                dt = tvals[n+1] - tvals[n]
                dz = (z1 - z0)/dt
                t0 = tvals[n]
                # integrals on this edge
                exp0 = np.exp(-1j*omega*t0)
                expd = np.exp(-1j*omega*dt)
                
                I1 = (1 - expd)/(1j*omega)
                I2 = -(dt*expd)/(1j*omega) - (1 - expd)/(omega**2)
                
                contrib = exp0*( z0*I1 + dz*I2 )
                total += contrib
            coeffs[k] = total/period
    return coeffs

def reconstruct_from_coeffs(coeffs, t, period=2*np.pi):
    """
    Evaluate z(t) from Fourier coefficients.
    
    Parameters
    ----------
    coeffs : dict
        Dictionary {k: c_k} with integer modes as keys.
    t : array_like
        Points where to evaluate, in [0, period).
    period : float
        Period of the parameterization (default 2*pi).
    
    Returns
    -------
    z : ndarray of complex
        Reconstructed complex values z(t).
    """
    t = np.asarray(t, dtype=float)
    z = np.zeros_like(t, dtype=complex)
    for k, ck in coeffs.items():
        omega = 2*np.pi*k/period
        z += ck * np.exp(1j*omega*t)
    return z


# Square polygon
# verts = np.array([[1,1], [-1,1], [-1,-1], [1,-1]], dtype=float)
verts = np.array([[0,-1], [1,0], [0,1]], dtype=float)

K=3
coeffs = polygon_fourier_coeffs(verts, K=K)

for k in range(-K,K+1):
    print(f"k={k:2d}, c_k={coeffs[k]}")

# theta = np.linspace(0,2*pi, 501)
# x0 = cos(theta)
# y0 = sin(theta)

# a = 0.8
# r = 1 + a*cos(theta) #+ b*cos(2*theta)
# x = r*cos(theta)
# y = r*sin(theta)

# plt.figure(1)
# plt.clf()
# plt.plot(x0,y0)
# plt.plot(x,y)
# plt.plot([0],[0],'or')
# plt.grid(True)
# plt.axis('equal')

t = np.linspace(0,2*pi,1000)
z = reconstruct_from_coeffs(coeffs, t)
x = np.real(z)
y = np.imag(z)

# normalize area to one
c_arr = [coeffs[key] for key in range(-K, K+1)]
area = np.sum(pi*np.arange(-K, K+1)*np.abs(c_arr)**2)
print(f'area1={area}')

c_arr = c_arr/np.sqrt(area)
area2 = np.sum(pi*np.arange(-K, K+1)*np.abs(c_arr)**2)
print(f'area2={area2}')

n_coeffs = {key: val for key, val in coeffs.items()}

t = np.linspace(0,2*pi,1000)
z = reconstruct_from_coeffs(n_coeffs, t)
x = np.real(z)
y = np.imag(z)

p_verts = np.concatenate((verts,[verts[0]]))

plt.figure(1)
plt.clf()
plt.plot(p_verts[:,0], p_verts[:,1],'.-')
plt.plot(x,y)
plt.plot([0],[0],'or')
plt.grid(True)
plt.axis('equal')

