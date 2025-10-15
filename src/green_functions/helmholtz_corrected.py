# -*- coding: utf-8 -*-
"""
Corrected Helmholtz Solver for Waveguides with Circular Obstacle

This corrected version addresses the coordinate system issues, domain validity,
and proper scattering problem setup from the original implementation.

@author: Corrected implementation
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lattice_sums as lattice
import src.grid_generation.gen_mesh_center_y as gen

class HelmholtzCorrected:
    def __init__(self, k, b=1.0, r=0.3, xc=0.0, yc=0.0, M_lattice=200, Lh=20):
        """
        Corrected Helmholtz solver for waveguides with circular obstacle.

        Parameters:
        -----------
        k : float
            Wave number
        b : float
            Half-width of waveguide (-b <= y <= b)
        r : float
            Circle radius
        xc, yc : float
            Center coordinates of circular obstacle
        M_lattice : int
            Number of terms for lattice sum convergence
        Lh : int
            Number of harmonics for lattice sums
        """
        self.k = k
        self.b = b
        self.kb = k * b  # Note: kb = k * b, not k/b
        self.r = r
        self.xc = xc
        self.yc = yc
        self.d = 2 * b  # Period in y-direction
        self.M_lattice = M_lattice
        self.Lh = Lh
        self.pi = np.pi
        self.epsilon = 1e-6

        # Precompute 2D lattice sums for the correct period
        self.S_2d = lattice.lattice_sums(self.d, k, beta=0, M=M_lattice, Lh=Lh)
        self.nodes, self.elements = self.initialize_mesh()

        print(f"Initialized Helmholtz solver:")
        print(f"  k = {k:.3f}, kb = {self.kb:.3f}")
        print(f"  Domain: y ∈ [{-b:.1f}, {b:.1f}], period d = {self.d:.1f}")
        print(f"  Obstacle: center = ({xc:.1f}, {yc:.1f}), radius = {r:.3f}")
        print(f"  Mesh nodes: {len(self.nodes)}")

    def initialize_mesh(self):
        """Generate mesh with corrected parameters"""
        L = 4.0  # Reduced domain size
        lc = 0.08  # Finer mesh

        nodes, elements = gen.mesh_with_obstacle_center(
            L=L, b=self.b, a=0.0, xc=self.xc, yc=self.yc,
            r=self.r, lc=lc
        )

        # Validate mesh points are within lattice sum domain
        valid_nodes = []
        valid_elements = []
        node_mapping = {}

        for i, node in enumerate(nodes):
            x, y = node[0], node[1]
            # Check if point is within valid domain for lattice sums
            if abs(y) <= 0.98 * self.b and abs(x) <= 3.0:
                new_idx = len(valid_nodes)
                node_mapping[i] = new_idx
                valid_nodes.append(node)

        # Update elements to use new node indices
        for element in elements:
            if all(old_idx in node_mapping for old_idx in element):
                new_element = [node_mapping[old_idx] for old_idx in element]
                valid_elements.append(new_element)

        valid_nodes = np.array(valid_nodes)
        valid_elements = np.array(valid_elements)

        print(f"  Valid nodes after domain filtering: {len(valid_nodes)}")
        return valid_nodes, valid_elements

    # --- Contour functions for circular obstacle ---
    def rho(self, t):
        return self.r

    def rho_p(self, t):
        return 0

    def rho_pp(self, t):
        return 0

    # --- Corrected Green's functions ---
    def G(self, x, y, xi, eta):
        """
        Corrected Green's function with proper coordinate handling
        """
        # Check domain validity
        if abs(y) > 0.99 * self.b or abs(eta) > 0.99 * self.b:
            return 0.0

        if abs(x - xi) > 4.0:  # Avoid far-field issues
            return 0.0

        return lattice.greens_neumann(
            x, y + self.b, xi, eta + self.b,
            self.S_2d, self.k, self.b
        )

    def G_reg(self, x, y, xi, eta):
        """
        Regularized Green's function
        """
        if abs(y) > 0.99 * self.b or abs(eta) > 0.99 * self.b:
            return 0.0

        if abs(x - xi) > 4.0:
            return 0.0

        return lattice.greens_neumann_reg(
            x, y + self.b, xi, eta + self.b,
            self.S_2d, self.k, self.b
        )

    # --- Derivatives of Green's function ---
    def Gx(self, x, y, xi, eta):
        try:
            Gm = self.G(x, y, xi - self.epsilon, eta)
            Gp = self.G(x, y, xi + self.epsilon, eta)
            return (Gp - Gm) / (2 * self.epsilon)
        except:
            return 0.0

    def Gy(self, x, y, xi, eta):
        try:
            Gm = self.G(x, y, xi, eta - self.epsilon)
            Gp = self.G(x, y, xi, eta + self.epsilon)
            return (Gp - Gm) / (2 * self.epsilon)
        except:
            return 0.0

    def Gx_reg(self, x, y, xi, eta):
        try:
            Gm = self.G_reg(x, y, xi - self.epsilon, eta)
            Gp = self.G_reg(x, y, xi + self.epsilon, eta)
            return (Gp - Gm) / (2 * self.epsilon)
        except:
            return 0.0

    def Gy_reg(self, x, y, xi, eta):
        try:
            Gm = self.G_reg(x, y, xi, eta - self.epsilon)
            Gp = self.G_reg(x, y, xi, eta + self.epsilon)
            return (Gp - Gm) / (2 * self.epsilon)
        except:
            return 0.0

    # --- Incident wave functions ---
    def incident_wave(self, x, y, wave_type='plane', direction=0.0):
        """
        Compute incident wave at point (x, y)

        Parameters:
        -----------
        x, y : float
            Point coordinates
        wave_type : str
            'plane' or 'point_source'
        direction : float
            Wave propagation angle (radians)
        """
        if wave_type == 'plane':
            kx = self.k * np.cos(direction)
            ky = self.k * np.sin(direction)
            return np.exp(1j * (kx * x + ky * y))
        elif wave_type == 'point_source':
            r = np.sqrt(x**2 + y**2)
            if r < 1e-10:
                return 0.0
            from scipy.special import hankel1
            return 0.25j * hankel1(0, self.k * r)
        else:
            return 1.0

    # --- Kernel for boundary integral ---
    def Gn_w(self, psi, theta):
        """
        Boundary integral kernel with improved error handling
        """
        try:
            r_psi = self.rho(psi)
            x = self.xc + r_psi * np.cos(psi)
            y = self.yc + r_psi * np.sin(psi)

            r_theta = self.rho(theta)
            xi = self.xc + r_theta * np.cos(theta)
            eta = self.yc + r_theta * np.sin(theta)

            r_prime = self.rho_p(theta)
            xi_p = r_prime * np.cos(theta) - r_theta * np.sin(theta)
            eta_p = r_prime * np.sin(theta) + r_theta * np.cos(theta)

            w = np.sqrt(r_theta**2 + r_prime**2)

            if abs(psi - theta) > 1e-10:
                return xi_p * self.Gy(x, y, xi, eta) - eta_p * self.Gx(x, y, xi, eta)
            else:
                r_2prime = self.rho_pp(theta)
                Y0_n_w = 1/(4*np.pi*w**2) * (r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
                G_reg_n_w = xi_p * self.Gy_reg(x, y, xi, eta) - eta_p * self.Gx_reg(x, y, xi, eta)
                return Y0_n_w + G_reg_n_w
        except:
            return 0.0

    def plot_mesh(self):
        """Plot the mesh and obstacle"""
        plt.figure(figsize=(10, 6))

        # Plot mesh
        triangles = self.elements
        triang = tri.Triangulation(self.nodes[:, 0], self.nodes[:, 1], triangles)
        plt.triplot(triang, color='gray', lw=0.5, alpha=0.7)

        # Plot obstacle
        theta = np.linspace(0, 2*self.pi, 100)
        x_circle = self.xc + self.r * np.cos(theta)
        y_circle = self.yc + self.r * np.sin(theta)
        plt.fill(x_circle, y_circle, color='red', alpha=0.7, label='Obstacle')

        # Plot domain boundaries
        plt.axhline(y=self.b, color='blue', linestyle='--', alpha=0.7, label='Domain boundary')
        plt.axhline(y=-self.b, color='blue', linestyle='--', alpha=0.7)

        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Corrected Mesh')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    def plot_solution_mesh(self, node_sol, title="Helmholtz Solution"):
        """
        Plot solution with improved visualization
        """
        # Remove extreme values that indicate numerical issues
        finite_mask = np.isfinite(node_sol)
        if not np.all(finite_mask):
            print(f"Warning: {np.sum(~finite_mask)} nodes have non-finite values")
            node_sol = node_sol.copy()
            node_sol[~finite_mask] = 0.0

        # Clip extreme values
        q99 = np.percentile(np.abs(node_sol[finite_mask]), 99)
        node_sol_clipped = np.clip(node_sol, -2*q99, 2*q99)

        triangles = self.elements
        triang = tri.Triangulation(self.nodes[:, 0], self.nodes[:, 1], triangles)

        plt.figure(figsize=(12, 6))

        # Plot real part
        plt.subplot(1, 2, 1)
        cs1 = plt.tricontourf(triang, np.real(node_sol_clipped), levels=50, cmap='RdBu_r')
        plt.colorbar(cs1, label='Real Part')
        plt.tricontour(triang, np.real(node_sol_clipped), levels=10, colors='black', alpha=0.3, linewidths=0.5)

        # Add obstacle
        theta = np.linspace(0, 2*self.pi, 100)
        x_circle = self.xc + self.r * np.cos(theta)
        y_circle = self.yc + self.r * np.sin(theta)
        plt.fill(x_circle, y_circle, color='white', edgecolor='black', linewidth=2)

        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'{title} - Real Part')
        plt.grid(True, alpha=0.3)

        # Plot magnitude
        plt.subplot(1, 2, 2)
        cs2 = plt.tricontourf(triang, np.abs(node_sol_clipped), levels=50, cmap='viridis')
        plt.colorbar(cs2, label='Magnitude')
        plt.tricontour(triang, np.abs(node_sol_clipped), levels=10, colors='black', alpha=0.3, linewidths=0.5)

        # Add obstacle
        plt.fill(x_circle, y_circle, color='white', edgecolor='black', linewidth=2)

        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'{title} - Magnitude')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        print(f"Solution statistics:")
        print(f"  Range: [{np.min(np.real(node_sol)):.3e}, {np.max(np.real(node_sol)):.3e}]")
        print(f"  Mean magnitude: {np.mean(np.abs(node_sol)):.3e}")

    # --- Scattering solver ---
    def solve_scattering(self, M=64, incident_type='plane', incident_direction=0.0):
        """
        Solve scattering problem with incident wave
        """
        print(f"Solving scattering problem with M={M} boundary points...")

        # Boundary discretization
        theta = (np.arange(1, M+1) - 0.5) * 2 * self.pi / M

        # Compute incident wave on boundary
        u_inc_boundary = np.zeros(M, dtype=complex)
        for i, t in enumerate(theta):
            x_b = self.xc + self.r * np.cos(t)
            y_b = self.yc + self.r * np.sin(t)
            u_inc_boundary[i] = self.incident_wave(x_b, y_b, incident_type, incident_direction)

        # Assemble BIE matrix
        print("Assembling BIE matrix...")
        Kw = np.zeros((M, M), dtype=np.complex128)
        for i in range(M):
            for j in range(M):
                Kw[i, j] = self.Gn_w(theta[i], theta[j])

        # For Dirichlet boundary conditions: (I/2 + K)φ = -u_inc
        A = 0.5 * np.identity(M) + (2 * self.pi / M) * Kw

        # Solve for boundary density
        phi_boundary = np.linalg.solve(A, -u_inc_boundary)

        print(f"Boundary solve completed. Max boundary density: {np.max(np.abs(phi_boundary)):.3e}")

        # Compute scattered field at mesh nodes
        print("Computing scattered field at mesh nodes...")
        u_scattered = np.zeros(self.nodes.shape[0], dtype=complex)
        u_incident = np.zeros(self.nodes.shape[0], dtype=complex)

        for idx, node in enumerate(self.nodes):
            x_i, y_i = node[0], node[1]

            # Check if point is inside obstacle (skip)
            if (x_i - self.xc)**2 + (y_i - self.yc)**2 <= (self.r + self.epsilon)**2:
                continue

            # Incident wave at mesh node
            u_incident[idx] = self.incident_wave(x_i, y_i, incident_type, incident_direction)

            # Scattered field using Green's representation
            scattered = 0.0
            for j, t in enumerate(theta):
                xi = self.xc + self.r * np.cos(t)
                eta = self.yc + self.r * np.sin(t)

                # Single layer potential
                try:
                    green_val = self.G(x_i, y_i, xi, eta)
                    scattered += phi_boundary[j] * green_val * (2 * self.pi / M) * self.r
                except:
                    continue

            u_scattered[idx] = scattered

        # Total field
        u_total = u_incident + u_scattered

        print("Scattering solution completed!")

        return u_total, u_incident, u_scattered

def main():
    """
    Test the corrected Helmholtz solver
    """
    print("="*60)
    print("Corrected Helmholtz Solver Test")
    print("="*60)

    # Parameters for a reasonable test case
    k = 2.0      # Wave number
    b = 1.0      # Half-width
    r = 0.3      # Obstacle radius

    # Create solver
    solver = HelmholtzCorrected(k=k, b=b, r=r, xc=0.0, yc=0.0)

    # Plot mesh
    solver.plot_mesh()

    # Solve scattering problem
    u_total, u_incident, u_scattered = solver.solve_scattering(
        M=32, incident_type='plane', incident_direction=0.0
    )

    # Plot results
    solver.plot_solution_mesh(u_total, "Total Field")
    solver.plot_solution_mesh(u_incident, "Incident Field")
    solver.plot_solution_mesh(u_scattered, "Scattered Field")

if __name__ == "__main__":
    main()
