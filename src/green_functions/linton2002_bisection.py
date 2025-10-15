import numpy as np
import matplotlib.pyplot as plt
import lattice_sums as lattice

a = 0.2
b = 1.0
h = 0.7
d = 2*b
pi = np.pi
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Circle geometry
def rho(t): return a
def rho_p(t): return 0
def rho_pp(t): return 0

epsilon = 1e-6

def Gn_w(psi, theta, G, G_reg):
    r_psi = rho(psi)
    x, y = r_psi*cos(psi), r_psi*sin(psi)

    r_theta = rho(theta)
    xi, eta = r_theta*cos(theta), r_theta*sin(theta)

    r_prime = rho_p(theta)
    xi_p = r_prime*cos(theta) - r_theta*sin(theta)
    eta_p = r_prime*sin(theta) + r_theta*cos(theta)
    w = sqrt(r_theta**2 + r_prime**2)

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

# Function to compute determinant given kb
def determinant(kb):
    k = kb/b
    S_2d = lattice.lattice_sums(2*d, k, beta=0, M=200, Lh=20)
    G = lambda x,y,xi,eta: lattice.greens_dirichlet(x, y+b+h, xi, eta+b+h, S_2d, k, d)
    G_reg = lambda x,y,xi,eta: lattice.greens_dirichlet_reg(x, y+b+h, xi, eta+b+h, S_2d, k, d)

    Kw = np.zeros((M,M), dtype=np.complex128)
    for i in range(M):
        for j in range(M):
            Kw[i,j] = Gn_w(theta[i], theta[j], G, G_reg)

    A = np.identity(M) - 4*np.pi/M * Kw
    return np.linalg.det(A)

# Bisection method
def bisection(f, a, b, tol=1e-6, maxiter=50):
    fa, fb = np.real(f(a)), np.real(f(b))
    if fa*fb > 0:
        raise ValueError("No sign change in interval, cannot use bisection.")

    for _ in range(maxiter):
        mid = 0.5*(a+b)
        fm = np.real(f(mid))

        if abs(fm) < tol or (b-a)/2 < tol:
            return mid

        if fa*fm < 0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    return 0.5*(a+b)

# Example: find root between kb_left and kb_right
kb_left, kb_right = pi*0.496, pi*0.499  # choose interval where det changes sign
kb_root = bisection(determinant, kb_left, kb_right)
print("Root kb ≈", kb_root, " (normalized kb/pi =", kb_root/pi, ")")

# Optional: plot determinant around root
kb_test = np.linspace(kb_left, kb_right, 30)
det_vals = [np.real(determinant(kb)) for kb in kb_test]

plt.figure()
plt.plot(kb_test/pi, det_vals, 'o-')
plt.axhline(0, color='r')
plt.axvline(kb_root/pi, color='g', linestyle='--')
plt.text(kb_root/pi, 0, f"{kb_root/pi:.6f}",
         fontsize=10, ha='left', va='bottom', color='g')
plt.xlabel("$kb/\\pi$")
plt.ylabel("Re[det A]")
plt.grid(True)
plt.show()
