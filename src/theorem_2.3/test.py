# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import math
import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice
import dipole_theorem as dip

pi = np.pi

def X(t):
    return epsilon * -np.cos(t) + (beta/2)*np.cos(2*t)

def Y(t):
    return epsilon * np.sin(t) - (beta/2)*np.sin(2*t)

def Xp(t):
    return epsilon * np.sin(t) - beta*np.sin(2*t)

def Yp(t):
    return epsilon * np.cos(t) - beta*np.cos(2*t)

def Xpp(t):
    return epsilon * np.cos(t) - 2*beta*np.cos(2*t)

def Ypp(t):
    return epsilon * -np.sin(t) + 2*beta*np.sin(2*t)

def W(t):
    return np.sqrt(Xp(t)**2 + Yp(t)**2)

def Gn_w(psi, theta, G, G_reg):
    x, y = X(psi), Y(psi)
    xi, eta = X(theta), Y(theta)
    xi_p = Xp(theta)
    eta_p = Yp(theta)
    w = W(theta)

    epsilon = 1e-6
    def Gx(x, y, xi, eta):
        return (G(x, y, xi+epsilon, eta) - G(x, y, xi-epsilon, eta))/(2*epsilon)
    def Gy(x, y, xi, eta):
        return (G(x, y, xi, eta+epsilon) - G(x, y, xi, eta-epsilon))/(2*epsilon)
    def Gx_reg(x, y, xi, eta):
        return (G_reg(x, y, xi+epsilon, eta) - G_reg(x, y, xi-epsilon, eta))/(2*epsilon)
    def Gy_reg(x, y, xi, eta):
        return (G_reg(x, y, xi, eta+epsilon) - G_reg(x, y, xi, eta-epsilon))/(2*epsilon)

    if abs(psi-theta) > 1e-10:
        return xi_p*Gy(x, y, xi, eta) - eta_p*Gx(x, y, xi, eta)
    else:
        r2prime = Xpp(theta)*Yp(theta) - Ypp(theta)*Xp(theta)
        Y0_n_w = r2prime / (4*pi*w**2)
        G_reg_n_w = xi_p*Gy_reg(x, y, xi, eta) - eta_p*Gx_reg(x, y, xi, eta)
        return Y0_n_w + G_reg_n_w

def determinant(kb, a_i):
    k = kb/b
    S_2d = lattice.lattice_sums(2*d, k, beta=0, M=200, Lh=20)
    G = lambda x,y,xi,eta: lattice.greens_dirichlet(x, y+b+a_i, xi, eta+b+a_i, S_2d, k, d)
    G_reg = lambda x,y,xi,eta: lattice.greens_dirichlet_reg(x, y+b+a_i, xi, eta+b+a_i, S_2d, k, d)

    Kw = np.zeros((M,M), dtype=np.complex128)
    for i in range(M):
        for j in range(M):
            Kw[i,j] = Gn_w(theta[i], theta[j], G, G_reg)

    A = np.identity(M) - 4*np.pi/M * Kw
    return np.linalg.det(A)

def bisection(f, a_i, kb_left, kb_right, tol=1e-6, maxiter=50):
    fa, fb = np.real(f(kb_left, a_i)), np.real(f(kb_right, a_i))
    if fa*fb > 0:
        # raise ValueError("No sign change in interval, cannot use bisection.")
        return None

    for _ in range(maxiter):
        mid = 0.5*(kb_left+kb_right)
        fm = np.real(f(mid, a_i))

        if abs(fm) < tol or (kb_right-kb_left)/2 < tol:
            return mid

        if fa*fm < 0:
            kb_right, fb = mid, fm
        else:
            kb_left, fa = mid, fm
    return 0.5*(kb_left+kb_right)

# ===========================================
M = 32 # circle divisions
# theta, to be used by function determinant
theta = (np.arange(1, M+1)-0.5)*2*np.pi/M

b = 1.0
d = 2*b
beta=1e-2
epsilon = 0.001

mu = dip.dipole(epsilon, beta, 0, 0)
Lambda1 = pi**2 / (4 * b**2)
Lambda2 = pi**2 * 2**2 / (4 * b**2)
sigma_analytic = epsilon**2*(np.pi**3/b**3) * mu
s_analytic = - 2*np.log10(sigma_analytic)
k2_analytic = Lambda2 - sigma_analytic**2
kb_analytic = np.sqrt(k2_analytic)*b
# a1 = -beta/12
# a = epsilon*a1
a = 0

print("Lambda1,Lambda2", Lambda1, Lambda2)
print("kb_min,kb_max", np.sqrt(Lambda1-sigma_analytic**2)*b, np.sqrt(Lambda2-sigma_analytic**2)*b)
print(f"sigma_analytic={sigma_analytic}")
print(f"s_analytic={s_analytic}")
print(f"k2_analytic={k2_analytic}")
print(f"kb_analytic={kb_analytic}")

val = math.floor(s_analytic)
s = np.array([val-1, val-0.5, val, val+0.5, val+1])
sigma2 = 10**(-s.astype(float))
kb_test = b*np.sqrt(Lambda2 - sigma2)
det_vals = [np.real(determinant(kb, a)) for kb in kb_test]

convert = lambda s: b*np.sqrt(Lambda2 - 10**(-s))
f = lambda s, a: determinant(convert(s), a)
s_left, s_right = math.floor(s_analytic), math.ceil(s_analytic)+0.1
s_root = bisection(f, a, s_left, s_right)
print(s_root)
s_numeric = 10**(-0.5*s_root)

plt.figure()
plt.plot(s, det_vals, 'o-')
plt.axhline(0, color='r')
plt.xlabel("$-\log_{10}(\sigma^2)$")
plt.ylabel("Re[det A]")
plt.grid(True)
plt.show()

