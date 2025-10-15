# -*- coding: utf-8 -*-
"""
Trapped Modes Solver for Waveguides with Circular Obstacle

This solver finds trapped modes (bound states) around circular obstacles
in waveguides. Trapped modes are localized wave solutions that decay
exponentially away from the obstacle and exist at discrete eigenfrequencies
below the radiation threshold.

@author: Trapped modes implementation
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lattice_sums as lattice
import src.grid_generation.gen_mesh_center_y as gen

class TrappedModesSolver:
    def __init__(self, b=1.0, r=0.3, xc=0.0, yc=0.0, M_lattice=200, Lh=20):
        """
        Trapped modes solver for circular obstacle in waveguide.

        Parameters:
        -----------
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
        self.b = b
        self.r = r
        self.xc = xc
        self.yc = yc
        self.d = 2 * b  # Period in y-direction
        self.M_lattice = M_lattice
        self.Lh = Lh
        self.pi = np.pi
        self.epsilon = 1e-6

        # Cutoff frequency for first propagating mode
        self.k_cutoff = np.pi / (2 * b)

        print(f"Trapped Modes Solver initialized:")
        print(f"  Waveguide: y ∈ [{-b:.1f}, {b:.1f}]")
        print(f"  Obstacle: center = ({xc:.1f}, {yc:.1f}), radius = {r:.3f}")
        print(f"  Cutoff frequency: k_c = π/(2b) = {self.k_cutoff:.4f}")
        print(f"  Trapped modes exist for k < {self.k_cutoff:.4f}")

    def initialize_mesh(self, k, L=4.0, lc=0.06):
        """Generate mesh for given frequency"""
        nodes, elements = gen.mesh_with_obstacle_center(
            L=L, b=self.b, a=0.0, xc=self.xc, yc=self.yc,
            r=self.r, lc=lc
        )

        # Filter nodes to avoid lattice sum domain issues
        valid_nodes = []
        valid_elements = []
        node_mapping = {}

        for i, node in enumerate(nodes):
            x, y = node[0], node[1]
            # Keep nodes within reasonable bounds for lattice sums
            if abs(y) <= 0.95 * self.b and abs(x) <= L:
                new_idx = len(valid_nodes)
                node_mapping[i] = new_idx
                valid_nodes.append(node)

        # Update elements
        for element in elements:
            if all(old_idx in node_mapping for old_idx in element):
                new_element = [node_mapping[old_idx] for old_idx in element]
                valid_elements.append(new_element)

        return np.array(valid_nodes), np.array(valid_elements)

    # --- Contour functions for circular obstacle ---
    def rho(self, t):
        return self.r

    def rho_p(self, t):
        return 0

    def rho_pp(self, t):
        return 0

    # --- Green's functions ---
    def setup_greens_functions(self, k):
        """Setup Green's functions for given frequency"""
        self.k = k
        self.S_2d = lattice.lattice_sums(self.d, k, beta=0, M=self.M_lattice, Lh=self.Lh)

    def G(self, x, y, xi, eta):
        """Green's function with domain checking"""
        if abs(y) > 0.98 * self.b or abs(eta) > 0.98 * self.b:
            return 0.0
        if abs(x - xi) > 6.0:  # Avoid far-field issues
            return 0.0

        try:
            return lattice.greens_neumann(
                x, y + self.b, xi, eta + self.b,
                self.S_2d, self.k, self.b
            )
        except:
            return 0.0

    def G_reg(self, x, y, xi, eta):
        """Regularized Green's function"""
        if abs(y) > 0.98 * self.b or abs(eta) > 0.98 * self.b:
            return 0.0
        if abs(x - xi) > 6.0:
            return 0.0

        try:
            return lattice.greens_neumann_reg(
                x, y + self.b, xi, eta + self.b,
                self.S_2d, self.k, self.b
            )
        except:
            return 0.0

    # --- Derivatives ---
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

    # --- Boundary integral kernel ---
    def Gn_w(self, psi, theta):
        """Boundary integral kernel for trapped modes"""
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
                # Diagonal term with regularization
                r_2prime = self.rho_pp(theta)
                Y0_n_w = 1/(4*np.pi*w**2) * (r_theta*r_2prime - r_theta**2 - 2*r_prime**2)
                G_reg_n_w = xi_p * self.Gy_reg(x, y, xi, eta) - eta_p * self.Gx_reg(x, y, xi, eta)
                return Y0_n_w + G_reg_n_w
        except:
            return 0.0

    def find_trapped_modes(self, k_min=None, k_max=None, M=64, num_k=50):
        """
        Find trapped modes by scanning frequency range

        Parameters:
        -----------
        k_min, k_max : float
            Frequency range to scan. If None, uses reasonable defaults
        M : int
            Number of boundary collocation points
        num_k : int
            Number of frequency points to scan

        Returns:
        --------
        k_values : ndarray
            Frequency values scanned
        determinants : ndarray
            Determinant values (trapped modes occur where |det| is small)
        eigenvalues : list
            Eigenvalues for each frequency
        """
        if k_min is None:
            k_min = 0.1 * self.k_cutoff
        if k_max is None:
            k_max = 0.98 * self.k_cutoff

        print(f"Scanning for trapped modes:")
        print(f"  Frequency range: [{k_min:.4f}, {k_max:.4f}]")
        print(f"  Cutoff frequency: {self.k_cutoff:.4f}")
        print(f"  Boundary points: {M}")

        k_values = np.linspace(k_min, k_max, num_k)
        determinants = np.zeros(num_k, dtype=complex)
        eigenvalues = []

        # Boundary discretization
        theta = (np.arange(1, M+1) - 0.5) * 2 * self.pi / M

        for i, k in enumerate(k_values):
            if i % 10 == 0:
                print(f"  Processing k = {k:.4f} ({i+1}/{num_k})")

            # Setup Green's functions for this frequency
            self.setup_greens_functions(k)

            # Assemble BIE matrix
            Kw = np.zeros((M, M), dtype=np.complex128)
            for j in range(M):
                for l in range(M):
                    Kw[j, l] = self.Gn_w(theta[j], theta[l])

            # For Dirichlet BC: A = I - (2π/M) * K
            A = np.identity(M) - (2 * self.pi / M) * Kw

            # Compute determinant and eigenvalues
            det_A = np.linalg.det(A)
            eig_vals = np.linalg.eigvals(A)

            determinants[i] = det_A
            eigenvalues.append(eig_vals)

        return k_values, determinants, eigenvalues

    def compute_trapped_mode_field(self, k, M=64):
        """
        Compute trapped mode field for a specific frequency

        Parameters:
        -----------
        k : float
            Frequency (should be close to an eigenfrequency)
        M : int
            Number of boundary points

        Returns:
        --------
        nodes : ndarray
            Mesh node coordinates
        field : ndarray
            Trapped mode field at mesh nodes
        """
        print(f"Computing trapped mode field at k = {k:.6f}")

        # Check frequency is below cutoff
        if k >= self.k_cutoff:
            print(f"Warning: k = {k:.4f} >= k_cutoff = {self.k_cutoff:.4f}")
            print("This frequency may not support trapped modes!")

        # Setup
        self.setup_greens_functions(k)
        nodes, elements = self.initialize_mesh(k)

        # Boundary discretization
        theta = (np.arange(1, M+1) - 0.5) * 2 * self.pi / M

        # Assemble BIE matrix
        print("Assembling BIE matrix...")
        Kw = np.zeros((M, M), dtype=np.complex128)
        for i in range(M):
            for j in range(M):
                Kw[i, j] = self.Gn_w(theta[i], theta[j])

        A = np.identity(M) - (2 * self.pi / M) * Kw

        # Find eigenmode
        eig_vals, eig_vecs = np.linalg.eig(A)
        eig_val_idx = np.argmin(np.abs(eig_vals))
        phi_boundary = eig_vecs[:, eig_val_idx]

        print(f"Smallest eigenvalue: {eig_vals[eig_val_idx]:.6e}")

        # Compute field at mesh nodes using Green's representation
        print("Computing field at mesh nodes...")
        field = np.zeros(len(nodes), dtype=complex)

        for idx, node in enumerate(nodes):
            x_i, y_i = node[0], node[1]

            # Skip if inside obstacle
            if (x_i - self.xc)**2 + (y_i - self.yc)**2 <= (self.r + 1e-6)**2:
                field[idx] = 0.0
                continue

            # Compute field using single-layer potential
            field_val = 0.0
            for j, t in enumerate(theta):
                xi = self.xc + self.r * np.cos(t)
                eta = self.yc + self.r * np.sin(t)

                try:
                    green_val = self.G(x_i, y_i, xi, eta)
                    field_val += phi_boundary[j] * green_val * (2 * self.pi / M) * self.r
                except:
                    continue

            field[idx] = field_val

        return nodes, elements, field

    def plot_determinant_scan(self, k_values, determinants):
        """Plot determinant vs frequency to identify trapped modes"""
        plt.figure(figsize=(12, 4))

        plt.subplot(1, 2, 1)
        plt.plot(k_values, np.abs(determinants), 'b-', linewidth=2)
        plt.axvline(self.k_cutoff, color='r', linestyle='--',
                   label=f'Cutoff k_c = {self.k_cutoff:.3f}')
        plt.yscale('log')
        plt.xlabel('Frequency k')
        plt.ylabel('|det(A)|')
        plt.title('Determinant Magnitude')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.plot(np.real(determinants), np.imag(determinants), 'b.-')
        plt.plot(0, 0, 'ro', markersize=8, label='Origin')
        plt.xlabel('Re(det)')
        plt.ylabel('Im(det)')
        plt.title('Determinant in Complex Plane')
        plt.grid(True)
        plt.legend()
        plt.axis('equal')

        plt.tight_layout()
        plt.show()

        # Find potential trapped mode frequencies
        min_indices = []
        threshold = 0.1 * np.max(np.abs(determinants))

        for i in range(1, len(determinants)-1):
            if (np.abs(determinants[i]) < threshold and
                np.abs(determinants[i]) < np.abs(determinants[i-1]) and
                np.abs(determinants[i]) < np.abs(determinants[i+1])):
                min_indices.append(i)

        if min_indices:
            print("Potential trapped mode frequencies:")
            for idx in min_indices:
                print(f"  k = {k_values[idx]:.6f}, |det| = {np.abs(determinants[idx]):.6e}")
        else:
            print("No clear trapped mode frequencies found in this range")

        return min_indices

    def plot_trapped_mode(self, nodes, elements, field, k, title_suffix=""):
        """Plot trapped mode field"""
        # Filter out extreme values
        finite_mask = np.isfinite(field)
        if not np.all(finite_mask):
            print(f"Warning: {np.sum(~finite_mask)} nodes have non-finite values")
            field = field.copy()
            field[~finite_mask] = 0.0

        # Create triangulation
        triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements)

        plt.figure(figsize=(15, 5))

        # Real part
        plt.subplot(1, 3, 1)
        cs1 = plt.tricontourf(triang, np.real(field), levels=50, cmap='RdBu_r')
        plt.colorbar(cs1, label='Real Part')
        self._add_obstacle_to_plot()
        plt.axis('equal')
        plt.title(f'Trapped Mode - Real Part\nk = {k:.4f}{title_suffix}')
        plt.xlabel('x')
        plt.ylabel('y')

        # Imaginary part
        plt.subplot(1, 3, 2)
        cs2 = plt.tricontourf(triang, np.imag(field), levels=50, cmap='RdBu_r')
        plt.colorbar(cs2, label='Imaginary Part')
        self._add_obstacle_to_plot()
        plt.axis('equal')
        plt.title(f'Trapped Mode - Imaginary Part\nk = {k:.4f}{title_suffix}')
        plt.xlabel('x')
        plt.ylabel('y')

        # Magnitude
        plt.subplot(1, 3, 3)
        cs3 = plt.tricontourf(triang, np.abs(field), levels=50, cmap='viridis')
        plt.colorbar(cs3, label='Magnitude')
        self._add_obstacle_to_plot()
        plt.axis('equal')
        plt.title(f'Trapped Mode - Magnitude\nk = {k:.4f}{title_suffix}')
        plt.xlabel('x')
        plt.ylabel('y')

        plt.tight_layout()
        plt.show()

        print(f"Field statistics:")
        print(f"  Range (real): [{np.min(np.real(field)):.3e}, {np.max(np.real(field)):.3e}]")
        print(f"  Mean magnitude: {np.mean(np.abs(field)):.3e}")

    def _add_obstacle_to_plot(self):
        """Add obstacle circle to current plot"""
        theta = np.linspace(0, 2*np.pi, 100)
        x_circle = self.xc + self.r * np.cos(theta)
        y_circle = self.yc + self.r * np.sin(theta)
        plt.fill(x_circle, y_circle, color='white', edgecolor='black', linewidth=2)

def main():
    """
    Example usage of trapped modes solver
    """
    print("="*60)
    print("Trapped Modes Solver")
    print("="*60)

    # Create solver
    solver = TrappedModesSolver(b=1.0, r=0.3, xc=0.0, yc=0.0)

    # Scan for trapped modes
    k_values, determinants, eigenvalues = solver.find_trapped_modes(
        k_min=0.2, k_max=1.5, M=32, num_k=100
    )

    # Plot results
    min_indices = solver.plot_determinant_scan(k_values, determinants)

    # Compute field at a promising frequency
    if min_indices:
        k_trapped = k_values[min_indices[0]]
        print(f"\nComputing trapped mode field at k = {k_trapped:.6f}")

        nodes, elements, field = solver.compute_trapped_mode_field(k_trapped, M=32)
        solver.plot_trapped_mode(nodes, elements, field, k_trapped)
    else:
        # Try a frequency close to cutoff
        k_test = 0.9 * solver.k_cutoff
        print(f"\nTesting frequency k = {k_test:.6f}")

        nodes, elements, field = solver.compute_trapped_mode_field(k_test, M=32)
        solver.plot_trapped_mode(nodes, elements, field, k_test, " (test)")

if __name__ == "__main__":
    main()
