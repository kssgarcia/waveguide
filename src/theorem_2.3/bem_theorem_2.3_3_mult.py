# %%
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

import math

import dipole_theorem as dip
import lattice_sums as lattice
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi


# For x-axis symmetry: Y(t) should be odd and X(t) should be even
def X(t):
    return epsilon * (np.cos(t) - (beta / 2) * np.cos(2 * t))

def Y(t):
    return epsilon * (np.sin(t) - (beta / 2) * np.sin(2 * t))


def Xp(t):
    return epsilon * (-np.sin(t) + beta * np.sin(2 * t))


def Xpp(t):
    return epsilon * (-np.cos(t) + 2 * beta * np.cos(2 * t))


def Yp(t):
    return epsilon * (np.cos(t) - beta * np.cos(2 * t))


def Ypp(t):
    return epsilon * (-np.sin(t) + 2 * beta * np.sin(2 * t))

def W(t):
    return np.sqrt(Xp(t) ** 2 + Yp(t) ** 2)

def Gn_w(psi, theta, G, G_reg):
    x, y = X(psi), Y(psi)
    xi, eta = X(theta), Y(theta)
    xi_p = Xp(theta)
    eta_p = Yp(theta)
    w = W(theta)

    epsilon = 1e-6

    def Gx(x, y, xi, eta):
        return (G(x, y, xi + epsilon, eta) - G(x, y, xi - epsilon, eta)) / (2 * epsilon)

    def Gy(x, y, xi, eta):
        return (G(x, y, xi, eta + epsilon) - G(x, y, xi, eta - epsilon)) / (2 * epsilon)

    def Gx_reg(x, y, xi, eta):
        return (G_reg(x, y, xi + epsilon, eta) - G_reg(x, y, xi - epsilon, eta)) / (
            2 * epsilon
        )

    def Gy_reg(x, y, xi, eta):
        return (G_reg(x, y, xi, eta + epsilon) - G_reg(x, y, xi, eta - epsilon)) / (
            2 * epsilon
        )

    if abs(psi - theta) > 1e-10:
        return xi_p * Gy(x, y, xi, eta) - eta_p * Gx(x, y, xi, eta)
    else:
        r2prime = Xpp(theta) * Yp(theta) - Ypp(theta) * Xp(theta)
        Y0_n_w = r2prime / (4 * pi * w**2)
        G_reg_n_w = xi_p * Gy_reg(x, y, xi, eta) - eta_p * Gx_reg(x, y, xi, eta)
        return Y0_n_w + G_reg_n_w


def determinant(kb, a_i):
    k = kb / b
    S_2d = lattice.lattice_sums(2 * d, k, beta=0, M=200, Lh=20)
    G = lambda x, y, xi, eta: lattice.greens_dirichlet(
        x, y + b + a_i, xi, eta + b + a_i, S_2d, k, d
    )
    G_reg = lambda x, y, xi, eta: lattice.greens_dirichlet_reg(
        x, y + b + a_i, xi, eta + b + a_i, S_2d, k, d
    )

    Kw = np.zeros((M, M), dtype=np.complex128)
    for i in range(M):
        for j in range(M):
            Kw[i, j] = Gn_w(theta[i], theta[j], G, G_reg)

    A = np.identity(M) - 4 * np.pi / M * Kw
    return np.linalg.det(A)


def bisection(f, a_i, kb_left, kb_right, tol=1e-6, maxiter=50):
    fa, fb = np.real(f(kb_left, a_i)), np.real(f(kb_right, a_i))
    print("--------------------------------")
    print(fa, fb)
    print("--------------------------------")
    if fa * fb > 0:
        # raise ValueError("No sign change in interval, cannot use bisection.")
        return None

    for _ in range(maxiter):
        mid = 0.5 * (kb_left + kb_right)
        fm = np.real(f(mid, a_i))

        if abs(fm) < tol or (kb_right - kb_left) / 2 < tol:
            return mid

        if fa * fm < 0:
            kb_right, fb = mid, fm
        else:
            kb_left, fa = mid, fm
    return 0.5 * (kb_left + kb_right)


# ===========================================
M = 32  # circle divisions
# theta, to be used by function determinant
theta = (np.arange(1, M + 1) - 0.5) * 2 * np.pi / M

b = 1.0
d = 2 * b
beta = 0.1

epsilon_list = np.linspace(0.01, 0.2, 10)
kb_numeric_list = []
kb_analytic_list = []
sigma_analytic_list = []
error_list = []

for epsilon in epsilon_list:
    mu, nu = dip.dipole(beta, 0, 0, "x")
    Lambda1 = (pi/(2*b))**2
    Lambda2 = (pi/b)**2
    sigma_analytic = epsilon**2 * (np.pi**3 / b**3) * mu
    s_analytic = -2 * np.log10(sigma_analytic)
    k2_analytic = Lambda2 - sigma_analytic**2
    kb_analytic = np.sqrt(k2_analytic) * b
    a = 0

    kb_analytic_list.append(kb_analytic)
    sigma_analytic_list.append(sigma_analytic)

    kb_left = kb_analytic - 0.001 
    kb_right = kb_analytic + 0.001
    kb_test = np.linspace(kb_left, kb_right, 15)
    det_vals = [np.real(determinant(kb, a)) for kb in kb_test]

    sign_changes = False
    for i in range(len(det_vals) - 1):
        if det_vals[i] * det_vals[i + 1] <= 0:
            sign_changes = True
            kb_left = kb_test[i]
            kb_right = kb_test[i + 1]
            break

    if not sign_changes:
        kb_left = kb_analytic - 0.2
        kb_right = kb_analytic + 0.2
        kb_test = np.linspace(kb_left, kb_right, 30)
        det_vals = [np.real(determinant(kb, a)) for kb in kb_test]

        for i in range(len(det_vals) - 1):
            if det_vals[i] * det_vals[i + 1] <= 0:
                sign_changes = True
                kb_left = kb_test[i]
                kb_right = kb_test[i + 1]
                break

    f = lambda s, a: determinant(s, a)
    kb_root = bisection(f, a, kb_left, kb_right) if sign_changes else None
    print("kb_numeric", kb_root)

    error = np.abs((np.abs(kb_root) - np.abs(kb_analytic))/np.abs(kb_analytic))

    kb_numeric_list.append(kb_root)
    error_list.append(error)

# %%
print(kb_analytic_list)
print(kb_numeric_list)
print(sigma_analytic_list)
print(error_list)
# %%
plt.figure(figsize=(8, 5))

# --- First axis: σ_sol vs ε ---
fig, ax1 = plt.subplots(figsize=(8, 5))
color1 = 'tab:blue'
color3 = 'tab:green'

ax1.plot(epsilon_list, kb_numeric_list, 'o-', color=color1, label=r'$kb_{sol}$', markersize=6)
ax1.plot(epsilon_list, kb_analytic_list, 'o-', color=color3, label=r'$kb_{analytic}$', markersize=6)
ax1.set_xlabel(r"$\epsilon$", fontsize=12)
ax1.set_ylabel(r"$kb$", fontsize=12)
ax1.grid(True, which='both', linestyle='--', alpha=0.4)

ax1.legend(loc='upper left', frameon=True)

plt.tight_layout()
plt.savefig(f"./old_kb_comparison.png", dpi=200)
plt.show()
