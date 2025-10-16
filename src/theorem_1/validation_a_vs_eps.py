# %%
import math
import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice
import pandas as pd

pi = np.pi
cos, sin, sqrt = np.cos, np.sin, np.sqrt

def rho(t): return r0
def rho_p(t): return 0
def rho_pp(t): return 0

def Gn_w(psi, theta, G, G_reg):
    r_psi = rho(psi)
    x, y = r_psi*cos(psi), r_psi*sin(psi)

    r_theta = rho(theta)
    xi, eta = r_theta*cos(theta), r_theta*sin(theta)

    r_prime = rho_p(theta)
    xi_p = r_prime*cos(theta) - r_theta*sin(theta)
    eta_p = r_prime*sin(theta) + r_theta*cos(theta)
    w = sqrt(r_theta**2 + r_prime**2)

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
        r_2prime = rho_pp(theta)
        Y0_n_w = 1/(4*pi*w**2) * (r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
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
        raise ValueError("No sign change in interval, cannot use bisection.")

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

def find_sigma(epsilon, a):
    global r0
    r0 = epsilon*R0
    alpha = pi * a / b  # α = π a / b
    sigma_analytic = (epsilon**2 * pi**2) / (4 * b**3) * (pi * mu * np.sin(alpha/2)**2 - 0.5 * S * np.cos(alpha/2)**2)
    s_analytic = - 2*np.log10(sigma_analytic)

    try:
        convert = lambda s: b*np.sqrt(Lambda1 - 10**(-s))
        f = lambda s, a: determinant(convert(s), a)
        s_left, s_right = math.floor(s_analytic), math.ceil(s_analytic)+0.1
        s_root = bisection(f, a, s_left, s_right)
    except:
        s_root = 0

    return s_root

def bisection_epsilon(f, a_i, epsilon_init, step=0.1, tol=1e-6, maxiter=50):
    fa = np.real(f(epsilon_init, a_i))
    if fa != 0:
        return None

    epsilon_sol = epsilon_init-step
    for _ in range(maxiter):
        fm = np.real(f(epsilon_sol, a_i))

        if fm == 0:
            epsilon_sol -= step
        else:
            epsilon_sol += step/2
            step = step/2

        if step <= tol:
            return epsilon_sol

    return epsilon_sol

# ===========================================
M = 32 # circle divisions
# theta, to be used by function determinant
theta = (np.arange(1, M+1)-0.5)*2*np.pi/M

b = 1.0
d = 2*b
R0 = 1

Lambda1 = pi**2 / (4 * b**2)
S = pi * R0**2
mu = R0**2
a0 = (2*b/pi)*np.arctan(S/(2*pi*mu))

epsilon_list = []
a_list = np.linspace(0.5, 0.7, 5)
for a in a_list:  # 4 points: 0.4, 0.5, 0.6, 0.7
    epsilon_max = bisection_epsilon(find_sigma, a, 0.5)
    print('-----------------------------------')
    print(epsilon_max)
    epsilon_list.append(epsilon_max)

plt.figure(figsize=(8, 5))
plt.plot(a_list, epsilon_list, 'o-', color='tab:blue', markersize=6)

# Add labels for each point
for a, eps in zip(a_list, epsilon_list):
    plt.text(a, eps, f"{eps:.3f}", fontsize=10, ha='left', va='bottom')

plt.xlabel(r"$a$", fontsize=12)
plt.ylabel(r"$\varepsilon_{max}$", fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()