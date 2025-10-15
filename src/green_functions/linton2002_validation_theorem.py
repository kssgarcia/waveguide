# %%
import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice

epsilon = 1e-1            # small parameter (scale)
r0 = 2 * epsilon          # obstacle nominal radius (circle)
b = 1.0
d = 2*b
pi = np.pi
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# BEM discretization
M = 32
theta = (np.arange(1, M+1)-0.5)*2*np.pi/M

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

M = 32
theta = (np.arange(1, M+1)-0.5)*2*np.pi/M

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

    # Store evaluation points
    eval_points = [kb_left, kb_right]
    eval_values = [fa, fb]

    if fa*fb > 0:
        # raise ValueError("No sign change in interval, cannot use bisection.")
        return None, None, None

    for _ in range(maxiter):
        mid = 0.5*(kb_left+kb_right)
        fm = np.real(f(mid, a_i))

        # Store the midpoint and its value
        eval_points.append(mid)
        eval_values.append(fm)

        if abs(fm) < tol or (kb_right-kb_left)/2 < tol:
            return mid, eval_points, eval_values

        if fa*fm < 0:
            kb_right, fb = mid, fm
        else:
            kb_left, fa = mid, fm

    return 0.5*(kb_left+kb_right), eval_points, eval_values

S = pi * r0**2          # area for a circle radius r0
mu = r0**2              # dipole vertical strength for circle (approx / exact)
a0_star = 2*b/pi * np.arctan(np.sqrt(S/(2*pi*mu)))   # formula (2.4)

print("Analytic a0* (from (2.4)) =", a0_star)

# choose three test heights: below, at, above
delta = 0.02 * b
a_tests = [max(-b+1e-6, a0_star - delta), a0_star, min(b-1e-6, a0_star + delta)]
print("Testing a values (below, analytic, above):", a_tests)

# kb search window: near pi/2 (kb = k*b), the first threshold sqrt(Lambda1) * b = pi/2
kb_center = pi/2
# we search slightly below kb_center (bound state has k < sqrt(Lambda1))
kb_min = pi * 0.45
kb_max = pi * 0.5

results = []
for a_i in a_tests:
    print("\n--- a_i =", a_i, "---")
    kb_root, kb_grid, det_vals = bisection(determinant, a_i, kb_min, kb_max)
    if kb_root is None:
        print("No root found in the interval for this a.")
        results.append((a_i, None))
    else:
        k_num = kb_root / b
        k2_num = k_num**2
        print("Found root kb =", kb_root, "=> k =", k_num, "=> k^2 =", k2_num)
        results.append((a_i, kb_root))

    # Plot determinant curve for diagnosis
    plt.figure(figsize=(6,3))
    plt.plot(np.array(kb_grid)/pi, np.real(np.array(det_vals)), '-o', markersize=3)
    plt.axhline(0, color='k', linestyle='--')
    if kb_root is not None:
        plt.axvline(kb_root/pi, color='g', linestyle='--', label=f'kb_root/pi={kb_root/pi:.6f}')
    plt.xlabel('kb / pi')
    plt.ylabel('Re[det A]')
    plt.title(f'a = {a_i:.5f}')
    plt.legend()
    plt.grid(True)
    plt.show()

# -------------------------
# analytic sigma and k^2 for the analytic a*
# -------------------------
alpha = pi * a0_star / b
sigma_analytic = (epsilon**2 * pi**2) / (4 * b**3) * (pi * mu * np.sin(alpha/2)**2 - 0.5 * S * np.cos(alpha/2)**2)
k2_analytic = (pi**2/(4*b**2)) - sigma_analytic**2

print("\nAnalytic quantities at a0*:")
print("sigma_analytic =", sigma_analytic)
print("k2_analytic =", k2_analytic)
print("sqrt(Lambda1) * b (kb threshold) = pi/2 =", pi/2)

# print summary of results
print("\nSummary of numeric root search:")
for a_i, root in results:
    if root is None:
        print(f"a={a_i:.6f} -> no root found in search interval")
    else:
        k_num = root / b
        print(f"a={a_i:.6f} -> kb_root = {root:.6f}, k^2_num = {k_num**2:.8f}, k^2_analytic = {k2_analytic:.8f}")
