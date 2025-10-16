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
R0 = 1

Lambda1 = pi**2 / (4 * b**2)
S = pi * R0**2
mu = R0**2
a0 = (2*b/pi)*np.arctan(S/(2*pi*mu))

# Predefine lists
epsilon_list = np.linspace(0.5, 3.162e-4, 20)
s_list = []
error_list = []
sigma_analytic_list = []
s_analytic_list = []

for i, epsilon in enumerate(epsilon_list):
    r0 = epsilon*R0
    a = 0.6
    alpha = pi * a / b  # α = π a / b
    sigma_analytic = (epsilon**2 * pi**2) / (4 * b**3) * (pi * mu * np.sin(alpha/2)**2 - 0.5 * S * np.cos(alpha/2)**2)
    s_analytic = - 2*np.log10(sigma_analytic)
    k2_analytic = Lambda1 - sigma_analytic**2
    kbpi_analytic = np.sqrt(k2_analytic)*b/pi

    print(f"-----------------{i}-----------------")
    print(f'a0={a0}')
    print(f"sigma_analytic={sigma_analytic}")
    print(f"s_analytic={s_analytic}")
    print(f"k2_analytic={k2_analytic}")
    print(f"kb/pi analytic={kbpi_analytic}")

    try:
        convert = lambda s: b*np.sqrt(Lambda1 - 10**(-s))
        f = lambda s, a: determinant(convert(s), a)
        s_left, s_right = math.floor(s_analytic), math.ceil(s_analytic)+0.1
        s_root = bisection(f, a, s_left, s_right)
        s_numeric = 10**(-0.5*s_root)
        error = np.abs((np.abs(s_root) - np.abs(s_analytic))/np.abs(s_analytic))
    except:
        s_root = 0
        error = 1

    print("s_numeric", s_root)
    print("error", error)

    # Append results
    s_list.append(s_root)
    error_list.append(error)
    sigma_analytic_list.append(sigma_analytic)
    s_analytic_list.append(s_analytic)


# Save everything in a CSV
df = pd.DataFrame({
    "epsilon": epsilon_list,
    "sigma_analytic": sigma_analytic_list,
    "s_analytic": s_analytic_list,
    "s_root": s_list,
    "error": error_list,
})
df.to_csv("kevin.csv", index=False)

plt.figure(figsize=(8, 5))

# --- First axis: σ_sol vs ε ---
fig, ax1 = plt.subplots(figsize=(8, 5))

color1 = 'tab:blue'
ax1.plot(epsilon_list, s_list, 'o-', color=color1, label=r'$\sigma_{sol}$', markersize=6)
ax1.axhline(0, color='r', linewidth=1)
ax1.set_xlabel(r"$\epsilon$", fontsize=12)
ax1.set_ylabel(r"$\sigma_{sol}$", color=color1, fontsize=12)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, which='both', linestyle='--', alpha=0.4)

# --- Second axis: Error vs ε ---
ax2 = ax1.twinx()
color2 = 'tab:orange'
ax2.plot(epsilon_list, error_list, 's--', color=color2, label='Error', markersize=5, linewidth=2)
ax2.set_ylabel("Error", color=color2, fontsize=12)
ax2.tick_params(axis='y', labelcolor=color2)

# --- Combined Legend ---
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', frameon=True)

# --- Title & Styling ---
plt.title(r"Solution Stress $\sigma_{sol}$ and Error vs $\epsilon$", fontsize=14, pad=10)
plt.tight_layout()
plt.show()
# %%
# Se demostro que el limite para la a=0.6 es de alrededor de 0.295
# si el a disminuye el limite se presenta en un valor cada vez mayor
