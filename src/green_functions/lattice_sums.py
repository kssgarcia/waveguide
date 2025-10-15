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
J0 = lambda z: jv(0,z)

def lattice_sums(d, k, beta, M, Lh):
    """
    Compute lattice sums
    """
    if not abs(beta) < k:
        print('abs(beta) should be less than k')
        sys.exit(1)
    
    p = 2*pi/d

    def gamma_m(beta_m):
        arg = beta_m**2 - k**2
        if arg >= 0:
            return sqrt(arg)
        else:
            return -1j*sqrt(-arg)


    # ============ lattice sums computation ============
    S = np.zeros(2*Lh + 1, dtype=np.complex128)

    #---------- compute S0 -------------------
    gamma_0 = gamma_m(beta) 
    zeta3 = zeta(3)
    C = np.euler_gamma

    S0 = -1 -2*1j/pi*(C + np.log(k/(2*p))) - 2*1j/(gamma_0*d) \
        - 2*1j*(k**2 + 2*beta**2)*zeta3/(p**3*d)

    sum = 0
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
    # print(f'S_0 = {S0}')
    S[0] = S0

    # preparing to define function that computes theta_m
    # Nk = int((k - beta)/p)
    # Mk = int((k + beta)/p)

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
        # print(f'S_{2*l-1} = {S2lm1}')

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
        # print(f'S_{2*l} = {S2l}')
        S[2*l] = S2l

    return S


def greens_periodic(X, Y, S, k, d):
    """
    Computes the periodic Green's function
    """
    r = sqrt(X**2 + Y**2)
    
    if r > 0.99*d:
        print('ERROR: r > 0.99*d')
        sys.exit(1)
        
    sum = H0(k*r)
    theta = np.arctan2(Y, X)
    L = len(S)

    for l in range(L):
        epsilon = 1 if l == 0 else 2
        sum += epsilon*S[l]*jv(l, k*r)*cos(l*(pi/2-theta))

    return -1j/4*sum

def greens_periodic_reg(X, Y, S, k, d):
    """
    Computes the periodic Green's function
    without the singular term H0(k*r)
    """
    r = sqrt(X**2 + Y**2)
    
    if r > 0.99*d:
        print('ERROR: r > 0.99*d')
        sys.exit(1)
        
    sum = J0(k*r) # not adding 1j*Y0(k*r) term
    theta = np.arctan2(Y, X)
    L = len(S)

    for l in range(L):
        epsilon = 1 if l == 0 else 2
        sum += epsilon*S[l]*jv(l, k*r)*cos(l*(pi/2-theta))

    return -1j/4*sum

def greens_neumann(x, y, xi, eta, S_2d, k, d):
    """
    Computes Green's function that satisfies
    Neumann boundary conditions
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_2d, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2*d)
    return term1 + term2

def greens_neumann_reg(x, y, xi, eta, S_2d, k, d):
    """
    Computes Green's function that satisfies
    Neumann boundary conditions without the
    singular term 0.25*Y0(k*r)
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_2d, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2*d)
    return term1 + term2
    
def greens_dirichlet(x, y, xi, eta, S_2d, k, d):
    """
    Computes Green's function that satisfies
    Dirichlet boundary conditions
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_2d, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2*d)
    return term1 - term2

def greens_dirichlet_reg(x, y, xi, eta, S_2d, k, d):
    """
    Computes Green's function that satisfies
    Dirichlet boundary conditions without the
    singular term 0.25*Y0(k*r)
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_2d, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2*d)
    return term1 - term2

def greens_dir_neu(x, y, xi, eta, S_dn, k, d):
    """
    Computes Green's function that satisfies
    a Dirichlet boundary condition for y = 0 and
    a Neumann boundary condition for y = d
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_dn, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_dn, k, 2*d)
    return term1 - term2

def greens_dir_neu_reg(x, y, xi, eta, S_dn, k, d):
    """
    Computes Green's function that satisfies
    a Dirichlet boundary condition for y = 0 and
    a Neumann boundary condition for y = d
    without the singular term 0.25*Y0(k*r)
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_dn, k, 2*d)
    term2 = greens_periodic(X, y + eta, S_dn, k, 2*d)
    return term1 - term2


if __name__ == "__main__":

    # values in Table 2 of linton1998greens.pdf
    d = 1
    k = 2
    beta = np.sqrt(2)
    M = 80
    Lh = 2 # l = 0, 1, ..., L = 2*lh + 1

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

    S = lattice_sums(d, k, beta, M, Lh)
    X = 0
    Y = 0.01
    G = greens_periodic(X, Y, S, k, d)
    print(f"G_lattice = {G}")
    
