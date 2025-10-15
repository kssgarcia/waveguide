# -*- coding: utf-8 -*-
"""
Computation of Green's function by lattice sums

Created on Thu Aug 28 12:21:05 2025
@author: agarz
"""
import numpy as np
from scipy.special import hankel1, zeta, jv
from math import factorial
from sympy import bernoulli
import sys

pi =  np.pi
sqrt = np.sqrt
exp = np.exp
cos = np.cos

H0 = lambda z: hankel1(0,z)

# values in Table 2 of linton1998greens.pdf
d = 1
k = 2
beta = np.sqrt(2)
X = 0
Y = 0.01

# values in Table 3 of linton1998greens.pdf
# d = 1
# k = 10
# beta = 5*np.sqrt(2)
# X = 0
# Y = 0.01

# values in Table 4 of linton1998greens.pdf
# d = 1
# k = 2
# beta = 3
# X = 0
# Y = 0.01


theta = np.arctan2(Y, X)
p = 2*pi/d

# ======= sum of images computation Eq. (2.7) in linton1998greens.pdf
def G_images(X, Y):
    r0 = sqrt(X**2 + Y**2)
    G = H0(k*r0)
    M = 5000
    for m in range(1,M+1):
        rm_p = sqrt(X**2 + (Y - m*d)**2)
        rm_m = sqrt(X**2 + (Y + m*d)**2)
        G += H0(k*rm_p)*exp(1j*m*beta*d) + H0(k*rm_m)*exp(-1j*m*beta*d)

    G *= -1j/4
    return G

G = G_images(X, Y)

print('G images=',G)

# ============ spectral representation ==============
def gamma_m(beta_m):
    arg = beta_m**2 - k**2
    if arg >= 0:
        return sqrt(arg)
    else:
        return -1j*sqrt(-arg)

def G_spectral(X, Y):
    gamma_0 = gamma_m(beta)
    aX = abs(X)
    sum = exp(-gamma_0*aX)*exp(1j*beta*Y)/gamma_0

    M = 5000
    for m in range(1, M+1):
        beta_mp = beta + m*p
        gamma_mp = gamma_m(beta_mp)

        beta_mm = beta - m*p
        gamma_mm = gamma_m(beta_mm)

        term_p = exp(-gamma_mp*aX)*exp(1j*beta_mp*Y)/gamma_mp
        # print('term_p=', term_p)
        term_m = exp(-gamma_mm*aX)*exp(1j*beta_mm*Y)/gamma_mm
        # print('term_m=', term_m)
        sum += term_p + term_m

    G = -1/(2*d)*sum
    return G

G = G_spectral(X, Y)
print("G spec=", G)


# sys.exit(0)
# ============ lattice sums computation ============
Lh = 2 # number of coefficients is 2*Lh + 1
# Lh = 3 # Table 3
# Lh = 2 # Table 4
S_in = np.zeros(2*Lh + 1, dtype=np.complex128)

# ------ inefficient computation of coefficients -----------
M = 10000
for l in range(2*Lh + 1):
    print('l=',l)
    sum = 0
    for m in range(1, M + 1):
        sum += hankel1(l, m*k*d)*(exp(1j*m*beta*d) + (-1)**l*exp(-1j*m*beta*d))

    S_in[l] = sum

print(S_in)

# def beta_m(m):
#     return beta + m*p

# def gamma_m(m):
#     return np.sqrt(beta_m(m)**2 - k**2)


#---------- compute S0 -------------------
# below, last term to enable complex number arithmetic
gamma_0 = gamma_m(beta) 
zeta3 = zeta(3)
C = np.euler_gamma

#------- inefficient S0 ----------
S0 = -1 -2*1j/pi*(C + np.log(k/(2*p))) - 2*1j/(gamma_0*d)
sum = 0
M = 5000
for m in range(1, M + 1):
    beta_mp = beta + m*p
    gamma_mp = gamma_m(beta_mp)

    beta_mm = beta - m*p
    gamma_mm = gamma_m(beta_mm)

    sum += 1/gamma_mp + 1/gamma_mm - 2/(p*m)

S0 += -2*1j/d*sum
print('ineff S0=',S0)


S = np.zeros(2*Lh + 1, dtype=np.complex128)

# ------- efficient S0
S0 = -1 -2*1j/pi*(C + np.log(k/(2*p))) - 2*1j/(gamma_0*d) \
    - 2*1j*(k**2 + 2*beta**2)*zeta3/(p**3*d)

sum = 0
M = 80 # Table 2
# M = 300 # Table 3
# M = 200 # Table 4
for m in range(1, M + 1):
    beta_mp = beta + m*p
    gamma_mp = gamma_m(beta_mp)

    beta_mm = beta - m*p
    gamma_mm = gamma_m(beta_mm)

    common = - 1/(p*m) - (k**2 + 2*beta**2)/(2*p**3*m**3)
    term_p = 1/gamma_mp + common
    term_m = 1/gamma_mm + common

    sum += term_p + term_m

S0 += -2*1j/d*sum
print('eff S0=', S0)
S[0] = S0

# preparing to define function that computes theta_m
if not abs(beta) < k:
    print('abs(beta) should be less than k')
    sys.exit(1)

Nk = int((k - beta)/p)
Mk = int((k + beta)/p)
 
def theta_m(m):
    beta_mp = beta + m*p
    up = beta_mp/k
    if abs(up) < 1: # theta_mp must be real
        theta_mp = np.arcsin(up) + 0*1j
    else: # up > 1
        theta_mp = 0.5*pi - np.arccosh(up)*1j

    beta_mm = beta - m*p
    um = beta_mm/k
    if abs(um) < 1: # theta_mm is real
        theta_mm = np.arcsin(um) + 0*1j
    else: # um < -1
        theta_mm = -0.5*pi + np.arccosh(abs(um))*1j

    return beta_mp, theta_mp, beta_mm, theta_mm



theta_0 = np.arcsin(beta/k) + 0*1j # must be real

for l in range(1, Lh + 1):
    # ----------- compute S_{2l} ---------
    S2l = -2*1j*exp(-2*1j*l*theta_0)/(gamma_0*d) \
        -2*1j*(-1)**l/pi*(k/(2*p))**(2*l)*zeta(2*l+1) \
        + 1j/(l*pi)

    sum1 = 0
    for m in range(1, M + 1):
        beta_mp, theta_mp, beta_mm, theta_mm = theta_m(m)
        gamma_mp = gamma_m(beta_mp)
        gamma_mm = gamma_m(beta_mm)
        
        term1 = exp(-2*1j*l*theta_mp)/(gamma_mp*d)
        term2 = exp(2*1j*l*theta_mm)/(gamma_mm*d)
        term3 = -(-1)**l/(m*pi)*(k/(2*m*p))**(2*l)

        sum1 += term1 + term2 + term3

    S2l += -2*1j*sum1

    sum2 = 0
    for m in range(1, l+1):
        num = (-1)**m*2**(2*m)*factorial(l + m - 1)
        deno = factorial(2*m)*factorial(l-m)
        fac = (p/k)**(2*m)
        B2m = bernoulli(2*m, beta/p) # slow sympy function?

        sum2 += num/deno*fac*B2m

    S2l += 1j/pi*sum2
    print('l=',l, 'S_{2l}=', S2l)
    S[2*l] = S2l
    
    # ---------------- Compute S_{2l-1} -------------------
    S2lm1 = 2*1j*exp(-1j*(2*l-1)*theta_0)/(gamma_0*d) \
        + 2*(-1)**l*beta*d*l/(pi**2)*(k/(2*p))**(2*l-1)*zeta(2*l+1)

    sum1 = 0
    for m in range(1, M + 1):
        beta_mp, theta_mp, beta_mm, theta_mm = theta_m(m)
        gamma_mp = gamma_m(beta_mp)
        gamma_mm = gamma_m(beta_mm)

        term1 = exp(-1j*(2*l-1)*theta_mp)/(gamma_mp*d)
        term2 = -exp(1j*(2*l-1)*theta_mm)/(gamma_mm*d)
        term3 = 1j*(-1)**l*beta*d*l/(m**2*pi**2)*(k/(2*m*p))**(2*l-1)

        sum1 += term1 + term2 + term3

    S2lm1 += 2*1j*sum1

    sum2 = 0
    for m in range(l):
        num = (-1)**m*2**(2*m)*factorial(l+m-1)
        deno = factorial(2*m+1)*factorial(l-m-1)
        fac = (p/k)**(2*m+1)
        B2m1 = bernoulli(2*m+1, beta/p)

        sum2 += num/deno*fac*B2m1

    S2lm1 += -2/pi*sum2
    S[2*l-1] = S2lm1
    print('S_{2l-1}=', S2lm1)



r = sqrt(X**2 + Y**2)
sum = H0(k*r)
# S = S_in

for l in range(2*Lh + 1):
    if l == 0:
        epsilon = 1
    else:
        epsilon = 2
    
    sum += epsilon*S[l]*jv(l, k*r)*cos(l*(pi/2-theta))

G = -1j/4*sum

print("G lattice=", G)


    
