# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 14:31:37 2025

@author: agarz
"""
# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lattice_sums as lattice
import src.grid_generation.gen_mesh_center_y as gen

class Helmholtz:
    def __init__(self, k, b=1.0, r=0.2, yc=0.7, M_lattice=200, Lh=20):
        """
        Helmholtz solver for waveguides with circular obstacle.

        Parameters:
        k : float
            Wave number
        b : float
            Half-width of waveguide (-b <= y <= b)
        r : float
            Circle radius
        h : float
            Height of circle center above midline
        """
        self.k = k
        self.b = b
        self.kb = k/b
        self.r = r
        self.yc = yc
        self.d = 2*b
        self.M_lattice = M_lattice
        self.Lh = Lh
        self.pi = np.pi
        self.epsilon = 1e-6

        # Precompute 2D lattice sums
        self.S_2d = lattice.lattice_sums(2*self.d, k, beta=0, M=M_lattice, Lh=Lh)
        self.nodes, self.elements = self.initialize_mesh()

    def initialize_mesh(self):
        L=5.0
        a=0.0
        xc=0.0
        lc=0.1

        nodes, elements = gen.mesh_with_obstacle_center(L=L, b=self.b, a=a, xc=xc, yc=self.yc, r=self.r, lc=lc)
        return nodes, elements

    # --- Contour functions ---
    def rho(self, t):
        return self.r

    def rho_p(self, t):
        return 0

    def rho_pp(self, t):
        return 0

    # --- Green's functions ---
    def G(self, x, y, xi, eta):
        '''
        return lattice.greens_dirichlet(
            x, y + self.b + self.yc, xi, eta + self.b + self.yc,
            self.S_2d, self.k, self.d
        )
        whe remove the self.b translation because is moving the obstacle in y axes
        '''
        return lattice.greens_dirichlet(
            x, y + self.yc, xi, eta + self.yc,
            self.S_2d, self.k, self.d
        )

    def G_reg(self, x, y, xi, eta):
        '''
        return lattice.greens_dirichlet_reg(
            x, y + self.b + self.yc, xi, eta + self.b + self.yc,
            self.S_2d, self.k, self.d
        )
        whe remove the self.b translation because is moving the obstacle in y axes
        '''
        return lattice.greens_dirichlet_reg(
            x, y + self.yc, xi, eta + self.yc,
            self.S_2d, self.k, self.d
        )

    # --- Derivatives of Green's function ---
    def Gx(self, x, y, xi, eta):
        Gm = self.G(x, y, xi - self.epsilon, eta)
        Gp = self.G(x, y, xi + self.epsilon, eta)
        return (Gp - Gm) / (2 * self.epsilon)

    def Gy(self, x, y, xi, eta):
        Gm = self.G(x, y, xi, eta - self.epsilon)
        Gp = self.G(x, y, xi, eta + self.epsilon)
        return (Gp - Gm) / (2 * self.epsilon)

    def Gx_reg(self, x, y, xi, eta):
        Gm = self.G_reg(x, y, xi - self.epsilon, eta)
        Gp = self.G_reg(x, y, xi + self.epsilon, eta)
        return (Gp - Gm) / (2 * self.epsilon)

    def Gy_reg(self, x, y, xi, eta):
        Gm = self.G_reg(x, y, xi, eta - self.epsilon)
        Gp = self.G_reg(x, y, xi, eta + self.epsilon)
        return (Gp - Gm) / (2 * self.epsilon)

    # --- Kernel for boundary integral ---
    def Gn_w(self, psi, theta):
        r_psi = self.rho(psi)
        x = r_psi * np.cos(psi)
        y = r_psi * np.sin(psi)

        r_theta = self.rho(theta)
        xi = r_theta * np.cos(theta)
        eta = r_theta * np.sin(theta)

        r_prime = self.rho_p(theta)
        xi_p = r_prime * np.cos(theta) - r_theta * np.sin(theta)
        eta_p = r_prime * np.sin(theta) + r_theta * np.cos(theta)

        w = np.sqrt(r_theta**2 + r_prime**2)

        if abs(psi - theta) > 1e-10:
            return xi_p * self.Gy(x, y, xi, eta) - eta_p * self.Gx(x, y, xi, eta)
        else:
            r_2prime = self.rho_pp(theta)
            Y0_n_w = 1/(4*self.pi*w**2) * (r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
            G_reg_n_w = xi_p * self.Gy_reg(x, y, xi, eta) - eta_p * self.Gx_reg(x, y, xi, eta)
            return Y0_n_w + G_reg_n_w

    def Gn_internal(self, x, y, theta):
        r_theta = self.rho(theta)
        xi = r_theta * np.cos(theta)
        eta = r_theta * np.sin(theta)

        r_prime = self.rho_p(theta)
        xi_p = r_prime * np.cos(theta) - r_theta * np.sin(theta)
        eta_p = r_prime * np.sin(theta) + r_theta * np.cos(theta)

        w = np.sqrt(r_theta**2 + r_prime**2)

        if abs(np.sqrt((x - xi)**2 + (y - eta)**2)) > 1e-10:
            return xi_p*self.Gy(x, y, xi, eta) - eta_p*self.Gx(x, y, xi, eta)
        else:
            r_2prime = self.rho_pp(theta)
            Y0_n_w = 1/(4*np.pi*w**2)*(r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
            G_reg_n_w = xi_p*self.Gy_reg(x, y, xi, eta) - eta_p*self.Gx_reg(x, y, xi, eta)
            return Y0_n_w + G_reg_n_w

    def plot_solution_mesh(self, node_sol):
        triangles = self.elements  # shape (num_triangles, 3)
        triang = tri.Triangulation(self.nodes[:,0], self.nodes[:,1], triangles)

        plt.figure(figsize=(8,6))
        plt.tricontourf(triang, node_sol, cmap='viridis', levels=50)
        plt.colorbar(label='Potential')
        plt.triplot(triang, color='k', lw=0.5, alpha=0.3)  # optional mesh lines
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Helmholtz Solution')
        plt.axis('equal')
        plt.show()

    # --- Helmholtz solver ---
    def solve_helmholtz(self, M=32, plot=False):
        # Solve Helmholtz in the boundary
        theta = (np.arange(1, M+1) - 0.5) * 2 * self.pi / M
        A_list = []

        Kw = np.zeros((M, M), dtype=np.complex128)
        for i in range(M):
            for j in range(M):
                Kw[i, j] = self.Gn_w(theta[i], theta[j])

        A = np.identity(M) - 4 * self.pi / M * Kw
        detA = np.linalg.det(A)
        A_list.append(A)

        print(f'kb={self.kb}, det(A)={detA}')

        eig_vals, eig_vecs = np.linalg.eig(A)

        eig_val_idx = np.argmin(np.abs(eig_vals))
        phi_x = eig_vecs[eig_val_idx]

        # Compute the potential in the nodes of the mesh
        nodes_sol = np.zeros(self.nodes.shape[0])
        for idx, node in enumerate(self.nodes):
            x_i = node[0]
            y_i = node[1]
            sol = 0

            for j in range(M):
                try:
                    sol += phi_x[j] * self.Gn_internal(x_i, y_i, theta[j])
                except:
                    sol += 0

            nodes_sol[idx] = sol

        self.plot_solution_mesh(nodes_sol)

        if plot:
            # Plot eigenvalues
            plt.figure()
            plt.scatter(np.real(eig_vals), np.imag(eig_vals))
            plt.grid(True)
            plt.xlabel('Re(λ)')
            plt.ylabel('Im(λ)')
            plt.title('Eigenvalues of A')
            plt.show()

            # Plot first eigenvector
            plt.figure()
            plt.plot(np.real(eig_vecs[:, 0]), '.-')
            plt.grid(True)
            plt.xlabel('Index')
            plt.ylabel('Real part of eigenvector')
            plt.title('First Eigenvector')
            plt.show()

        return detA

helm = Helmholtz(k=1.564513141487717, b=1.0, r=0.2, yc=0.0)
detA = helm.solve_helmholtz(M=32)
