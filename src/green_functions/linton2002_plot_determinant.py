# %%
import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice
import os

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
print(f'a0={a0}')

epsilon = 1e-2
r0 = epsilon*R0

a = 0.6

# a_list = np.arange(0,b,0.2)

# kb_sol = 0
# ai_sol = 0
# for a_i in a_list:
#     try:
#         kb_root = bisection(determinant, a_i, kb_left, kb_right)
#     except:
#         continue
#     if kb_root:
#         kb_sol = kb_root
#         ai_sol = a_i
#         break
#     print("a_i", a_i)
#     print("kb", kb_root)
#     print("-------------")

alpha = pi * a / b  # α = π a / b
sigma_analytic = (epsilon**2 * pi**2) / (4 * b**3) * (pi * mu * np.sin(alpha/2)**2 - 0.5 * S * np.cos(alpha/2)**2)
print(f"sigma_analytic={sigma_analytic}")
s_analytic = - 2*np.log10(sigma_analytic)
print(f"s_analytic={s_analytic}")
k2_analytic = Lambda1 - sigma_analytic**2
print(f"k2_analytic={k2_analytic}")
kbpi_analytic = np.sqrt(k2_analytic)*b/pi
print(f"kb/pi analytic={kbpi_analytic}")

# print('Numeric')
# print((kb_sol/b)**2)
# kb_left, kb_right = pi*0.499, pi*0.49999
# kb_test = np.linspace(kb_left, kb_right, 30)
# s = np.array([6.4, 7.0])

if False: # disable once determinant value have been computed
    s = np.linspace(6.4, 7.0, num=51)
    sigma2 = 10**(-s)
    kb_test = b*np.sqrt(Lambda1 - sigma2)
    det_vals = [np.real(determinant(kb, a)) for kb in kb_test]
    np.savez("determinant_data", s=s, kb_test=kb_test, det_vals=det_vals)
else:
    data = np.load("determinant_data.npz")
    s = data['s']
    det_vals = data['det_vals']

data1 = np.load("circle_k2_2.npz")
s_root = data1['s_root']
sigma_numeric = data1['sigma_numeric']

plt.figure()
# plt.plot(kb_test/pi, det_vals, 'o-')
plt.plot(s, det_vals, '-', linewidth=2)
plt.tick_params(axis='both', labelsize=14)
plt.axhline(0, color='k', linestyle='--')
plt.axvline(s_root, color='k', linestyle='--')
# plt.text(kb_sol/pi, 0, f"{kb_sol/pi:.6f}",
#          fontsize=10, ha='left', va='bottom', color='g')
plt.xlabel("$-\log_{10}(\Lambda_1 - k^2)$", fontsize=20)
plt.ylabel("Re(det $A$)", fontsize=20)
# plt.grid(True)
plt.annotate("Trapped mode $k_{num}$", (s_root, 0.75), (6.8, 0.75),
             horizontalalignment='right', verticalalignment='center',
             fontsize=18,
             arrowprops=dict(arrowstyle='->', facecolor='black'))
plt.xlim([6.4, 7.0])

waves_path = os.environ["WAVES_PATH"]
file_path = os.path.join(waves_path,"Portugal 2025. AGL","Portugal 2025. AGL",
                         "determinant.pdf")
plt.savefig(file_path, bbox_inches='tight')


plt.show()

# convert = lambda s: b*np.sqrt(Lambda1 - 10**(-s))
# f = lambda s, a: determinant(convert(s), a)
# s_left, s_right = 8.0, 9.0
# s_root = bisection(f, a, s_left, s_right)

# sigma_numeric = 10**(-0.5*s_root)

# np.savez("circle_k2_5", b=b, a=a, epsilon=epsilon, M=M,
#          sigma_analytic=sigma_analytic, s_analytic=s_analytic,
#          s_root=s_root, sigma_numeric=sigma_numeric)

         
