# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lattice_sums as lattice
import src.grid_generation.gen_mesh_center_y as gen

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
    print(kb)
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

# ===========================================
M = 32 # circle divisions
# theta, to be used by function determinant
theta = (np.arange(1, M+1)-0.5)*2*np.pi/M

b = 1.0
d = 2*b

# Predefine lists
r0 = 0.3
a = 0.0

kb = 0.95
f = lambda s, a: determinant(s, a)
kb_left, kb_right = pi*(kb-0.01), pi*(kb+0.01)
print(kb_left, kb_right)
kb_root = bisection(f, a, kb_left, kb_right)

print("kb_numeric", kb_root)
print("kb_numeric", kb_root/pi)

def plot_solution_mesh(node_sol):
    triangles = elements  # shape (num_triangles, 3)
    triang = tri.Triangulation(nodes[:,0], nodes[:,1], triangles)

    plt.figure(figsize=(10, 4))  # Wider aspect ratio like the reference

    # Create symmetric levels around 0
    node_sol_real = np.real(node_sol)
    vmax = np.max(np.abs(node_sol_real))
    levels = np.linspace(-vmax, vmax, 50)

    # Plot with diverging colormap
    plt.tricontourf(triang, node_sol_real, levels=levels, cmap='RdBu_r', extend='both')
    cbar = plt.colorbar(label='u(x,y)')
    cbar.ax.tick_params(labelsize=10)

    # Optional: add contour lines for clarity
    plt.tricontour(triang, node_sol_real, levels=10, colors='black', linewidths=0.5, alpha=0.3)

    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.axis('equal')
    plt.xlim(-2, 2)  # Adjust based on your domain
    plt.ylim(-1.5, 1.5)
    plt.tight_layout()
    plt.show()


def Gn_internal( x, y, theta):
    r_theta = rho(theta)
    xi = r_theta * np.cos(theta)
    eta = r_theta * np.sin(theta)

    r_prime = rho_p(theta)
    xi_p = r_prime * np.cos(theta) - r_theta * np.sin(theta)
    eta_p = r_prime * np.sin(theta) + r_theta * np.cos(theta)

    w = np.sqrt(r_theta**2 + r_prime**2)

    epsilon = 1e-6
    def Gx(x, y, xi, eta):
        return (G(x, y, xi+epsilon, eta) - G(x, y, xi-epsilon, eta))/(2*epsilon)
    def Gy(x, y, xi, eta):
        return (G(x, y, xi, eta+epsilon) - G(x, y, xi, eta-epsilon))/(2*epsilon)
    def Gx_reg(x, y, xi, eta):
        return (G_reg(x, y, xi+epsilon, eta) - G_reg(x, y, xi-epsilon, eta))/(2*epsilon)
    def Gy_reg(x, y, xi, eta):
        return (G_reg(x, y, xi, eta+epsilon) - G_reg(x, y, xi, eta-epsilon))/(2*epsilon)

    if np.sqrt((x - xi)**2 + (y - eta)**2) > 1e-2:
        return xi_p*Gy(x, y, xi, eta) - eta_p*Gx(x, y, xi, eta)
    else:
        r_2prime = rho_pp(theta)
        Y0_n_w = 1/(4*np.pi*w**2)*(r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
        G_reg_n_w = xi_p*Gy_reg(x, y, xi, eta) - eta_p*Gx_reg(x, y, xi, eta)
        return Y0_n_w + G_reg_n_w

L=2.0
xc=0.0
lc=0.05
nodes, elements = gen.mesh_with_obstacle_center(L=L, b=b, xc=xc, yc=a, r=r0, lc=lc)

k = kb_root/b
S_2d = lattice.lattice_sums(2*d, k, beta=0, M=200, Lh=20)
G = lambda x,y,xi,eta: lattice.greens_dirichlet(x, y+b+a, xi, eta+b+a, S_2d, k, d)
G_reg = lambda x,y,xi,eta: lattice.greens_dirichlet_reg(x, y+b+a, xi, eta+b+a, S_2d, k, d)
Kw = np.zeros((M,M), dtype=np.complex128)
for i in range(M):
    for j in range(M):
        Kw[i,j] = Gn_w(theta[i], theta[j], G, G_reg)
A = np.identity(M) - 4*np.pi/M * Kw
det_A = np.linalg.det(A)

eig_vals, eig_vecs = np.linalg.eig(A)
eig_val_idx = np.argmin(np.abs(eig_vals))
phi_x = eig_vecs[eig_val_idx]

# Compute the potential in the nodes of the mesh
# for vec in eig_vecs:
nodes_sol = np.zeros(nodes.shape[0])
for idx, node in enumerate(nodes):
    x_i = node[0]
    y_i = node[1]
    # Skip nodes inside the circular obstacle (avoid singular evaluation)
    if (x_i - xc)**2 + (y_i - a)**2 <= (r0 + 1e-6)**2:
        nodes_sol[idx] = 0.0
        continue

    sol = 0

    for j in range(M):
        try:
            sol += phi_x[j] * Gn_internal(x_i, y_i, theta[j])
        except:
            pass

    # Include quadrature weight for boundary integral (uniform in theta)
    nodes_sol[idx] = (2*np.pi/M) * sol
plot_solution_mesh(nodes_sol)