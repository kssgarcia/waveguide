# -*- coding: utf-8 -*-
"""
Computation of Green's function by lattice sums.

Original version:
    Created on Thu Aug 28 12:21:05 2025
    @author: agarz

Updated for complex-k evaluation:
    This version preserves the original real-k formulas, but also supports
    complex k values needed by contour nonlinear-eigenvalue methods such as
    Beyn's method.

IMPORTANT ANALYTICITY ASSUMPTION
--------------------------------
For complex k, the contour must remain inside a single modal regime and must
not cross / enclose a Rayleigh threshold |beta_m|.  In the current trapped-mode
experiment this is exactly the situation:

    beta = 0,
    0 < Re(k) < pi/(2 b),

so the m=0 lattice channel is the only open channel and all m != 0 channels
remain closed throughout the Beyn contour.

The longitudinal factor gamma_m is continued analytically using

    open channel:   gamma_m = -i sqrt(k^2 - beta_m^2)
    closed channel: gamma_m =     sqrt(beta_m^2 - k^2)

with principal square roots, provided the channel classification is constant
along the contour.

The old theta_m construction used real comparisons and inverse trigonometric
branch choices.  It has been replaced by the algebraic identities

    exp(-i theta_m) = i (gamma_m - beta_m) / k
    exp(+i theta_m) = i (gamma_m + beta_m) / k,

which are equivalent on the original physical branch and are much better
suited to analytic continuation in complex k.
"""

from __future__ import annotations

from math import factorial

import numpy as np
from scipy.special import hankel1, jv, zeta
from sympy import bernoulli


pi = np.pi
sqrt = np.sqrt
exp = np.exp
cos = np.cos

H0 = lambda z: hankel1(0, z)
J0 = lambda z: jv(0, z)

_COMPLEX_TOL = 1.0e-13


def _as_complex_scalar(z) -> complex:
    """Return z as a Python complex scalar."""
    return complex(np.asarray(z).item())


def _is_effectively_real(z: complex, tol: float = _COMPLEX_TOL) -> bool:
    return abs(z.imag) <= tol * max(1.0, abs(z.real))


def lattice_sums(d, k, beta, M, Lh):
    """
    Compute lattice sums.

    Parameters
    ----------
    d : float
        Period.
    k : float or complex
        Wavenumber. Complex values are supported for contour methods.
    beta : float
        Bloch / quasi-periodic parameter. The current waveguide calculations
        use beta = 0.
    M : int
        Number of lattice terms.
    Lh : int
        Highest harmonic pair index.

    Notes
    -----
    For real k, this function reproduces the original physical branch:
        gamma_m > 0 for evanescent channels,
        Im(gamma_m) < 0 for propagating channels.

    For complex k, the contour must not cross a modal threshold.  The open /
    closed channel classification therefore has to remain constant along the
    contour.  This is satisfied by the current Beyn contour below Lambda_1.
    """
    d = float(d)
    k = _as_complex_scalar(k)
    beta_c = _as_complex_scalar(beta)

    if not _is_effectively_real(beta_c):
        raise ValueError(
            "This lattice-sum implementation expects real beta. "
            f"Received beta={beta_c!r}."
        )
    beta = float(beta_c.real)

    if abs(k) == 0.0:
        raise ValueError("k must be non-zero because the lattice sums contain log(k) and 1/k.")

    # Preserve the original admissibility test on the real axis.
    # Do NOT replace this with abs(beta) < abs(k): |k| is non-analytic.
    if _is_effectively_real(k):
        if not abs(beta) < k.real:
            raise ValueError(
                f"abs(beta) should be less than k on the physical real axis; "
                f"got beta={beta}, k={k.real}."
            )
    else:
        # For the complex continuation used here, the contour is assumed to be
        # a deformation of a positive-real-k interval in the same modal regime.
        if k.real <= abs(beta):
            raise ValueError(
                "Complex-k continuation requires Re(k) > abs(beta) for the "
                "current physical branch. Move the contour or revise the branch "
                f"construction. Got beta={beta}, k={k}."
            )

    p = 2.0 * pi / d

    def gamma_m(beta_m: float) -> complex:
        """
        Longitudinal lattice factor on the physical branch.

        Classification by Re(k) is legitimate only because the Beyn contour is
        required not to cross a threshold |beta_m|.  Within such a contour the
        classification is constant, while the formula itself remains analytic
        in k.
        """
        beta_m = float(beta_m)

        if abs(beta_m) < k.real:
            # Propagating/open channel on the reference real interval.
            # For real k this gives -i*sqrt(k^2-beta_m^2), exactly as before.
            q = np.sqrt(k * k - beta_m * beta_m + 0.0j)
            return -1j * q

        # Evanescent/closed channel on the reference real interval.
        # For real k this is the positive real square root.
        return np.sqrt(beta_m * beta_m - k * k + 0.0j)

    def exp_minus_i_theta(beta_m: float, gamma: complex) -> complex:
        # e^{-i theta} = cos(theta) - i sin(theta)
        # with sin(theta)=beta_m/k and cos(theta)=i*gamma/k.
        return 1j * (gamma - beta_m) / k

    def exp_plus_i_theta(beta_m: float, gamma: complex) -> complex:
        # e^{+i theta} = cos(theta) + i sin(theta)
        return 1j * (gamma + beta_m) / k

    # ============ lattice sums computation ============
    S = np.zeros(2 * Lh + 1, dtype=np.complex128)

    # ---------- compute S0 ----------
    gamma_0 = gamma_m(beta)
    zeta3 = float(zeta(3))
    C = np.euler_gamma

    # np.log is the complex principal logarithm.  The Beyn contour used for
    # Theorem 2.1 stays in Re(k)>0 and does not encircle k=0, hence this branch
    # remains consistent on the contour.
    S0 = (
        -1
        - 2j / pi * (C + np.log(k / (2.0 * p)))
        - 2j / (gamma_0 * d)
        - 2j * (k**2 + 2.0 * beta**2) * zeta3 / (p**3 * d)
    )

    lattice_sum = 0.0j
    for m in range(1, M + 1):
        beta_mp = beta + m * p
        gamma_mp = gamma_m(beta_mp)

        beta_mm = beta - m * p
        gamma_mm = gamma_m(beta_mm)

        common = -1.0 / (p * m) - (k**2 + 2.0 * beta**2) / (2.0 * p**3 * m**3)
        term_p = 1.0 / gamma_mp + common
        term_m = 1.0 / gamma_mm + common

        lattice_sum += term_p + term_m

    S0 += -2j / d * lattice_sum
    S[0] = S0

    # Instead of theta_0 = arcsin(beta/k), use the analytic exponential
    # directly. This avoids non-analytic abs(...) branch tests.
    e_minus_theta_0 = exp_minus_i_theta(beta, gamma_0)

    for l in range(1, Lh + 1):
        # ---------------- Compute S_{2l-1} ----------------
        odd_order = 2 * l - 1

        S2lm1 = (
            2j * e_minus_theta_0**odd_order / (gamma_0 * d)
            + 2.0
            * (-1) ** l
            * beta
            * d
            * l
            / (pi**2)
            * (k / (2.0 * p)) ** odd_order
            * zeta(2 * l + 1)
        )

        sum1 = 0.0j
        for m in range(1, M + 1):
            beta_mp = beta + m * p
            gamma_mp = gamma_m(beta_mp)

            beta_mm = beta - m * p
            gamma_mm = gamma_m(beta_mm)

            e_minus_mp = exp_minus_i_theta(beta_mp, gamma_mp)
            e_plus_mm = exp_plus_i_theta(beta_mm, gamma_mm)

            term1 = e_minus_mp**odd_order / (gamma_mp * d)
            term2 = -(e_plus_mm**odd_order) / (gamma_mm * d)
            term3 = (
                1j
                * (-1) ** l
                * beta
                * d
                * l
                / (m**2 * pi**2)
                * (k / (2.0 * m * p)) ** odd_order
            )

            sum1 += term1 + term2 + term3

        S2lm1 += 2j * sum1

        sum2 = 0.0j
        for m in range(l):
            num = (-1) ** m * 2 ** (2 * m) * factorial(l + m - 1)
            deno = factorial(2 * m + 1) * factorial(l - m - 1)
            fac = (p / k) ** (2 * m + 1)

            # Convert the SymPy number explicitly so complex NumPy arithmetic
            # never falls back to object dtype.
            B2m1 = complex(bernoulli(2 * m + 1, beta / p))

            sum2 += num / deno * fac * B2m1

        S2lm1 += -2.0 / pi * sum2
        S[2 * l - 1] = complex(S2lm1)

        # ---------------- Compute S_{2l} ----------------
        even_order = 2 * l

        S2l = (
            -2j * e_minus_theta_0**even_order / (gamma_0 * d)
            - 2j
            * (-1) ** l
            / pi
            * (k / (2.0 * p)) ** even_order
            * zeta(2 * l + 1)
            + 1j / (l * pi)
        )

        sum1 = 0.0j
        for m in range(1, M + 1):
            beta_mp = beta + m * p
            gamma_mp = gamma_m(beta_mp)

            beta_mm = beta - m * p
            gamma_mm = gamma_m(beta_mm)

            e_minus_mp = exp_minus_i_theta(beta_mp, gamma_mp)
            e_plus_mm = exp_plus_i_theta(beta_mm, gamma_mm)

            term1 = e_minus_mp**even_order / (gamma_mp * d)
            term2 = e_plus_mm**even_order / (gamma_mm * d)
            term3 = (
                -(-1) ** l
                / (m * pi)
                * (k / (2.0 * m * p)) ** even_order
            )

            sum1 += term1 + term2 + term3

        S2l += -2j * sum1

        sum2 = 0.0j
        for m in range(1, l + 1):
            num = (-1) ** m * 2 ** (2 * m) * factorial(l + m - 1)
            deno = factorial(2 * m) * factorial(l - m)
            fac = (p / k) ** (2 * m)
            B2m = complex(bernoulli(2 * m, beta / p))

            sum2 += num / deno * fac * B2m

        S2l += 1j / pi * sum2
        S[2 * l] = complex(S2l)

    return S


def greens_periodic(X, Y, S, k, d):
    """
    Compute the periodic Green's function.

    k may be real or complex. scipy.special.hankel1 and jv both support
    complex arguments.
    """
    r = sqrt(X**2 + Y**2)

    if r > 0.99 * d:
        raise ValueError(f"greens_periodic requires r <= 0.99*d; got r={r}, d={d}.")

    series_sum = H0(k * r)
    theta = np.arctan2(Y, X)
    L = len(S)

    for l in range(L):
        epsilon = 1 if l == 0 else 2
        series_sum += epsilon * S[l] * jv(l, k * r) * cos(l * (pi / 2.0 - theta))

    return -1j / 4.0 * series_sum


def greens_periodic_reg(X, Y, S, k, d):
    """
    Compute the periodic Green's function without the singular H0 term.

    k may be real or complex.
    """
    r = sqrt(X**2 + Y**2)

    if r > 0.99 * d:
        raise ValueError(
            f"greens_periodic_reg requires r <= 0.99*d; got r={r}, d={d}."
        )

    series_sum = J0(k * r)  # not adding the singular i*Y0(k*r) contribution
    theta = np.arctan2(Y, X)
    L = len(S)

    for l in range(L):
        epsilon = 1 if l == 0 else 2
        series_sum += epsilon * S[l] * jv(l, k * r) * cos(l * (pi / 2.0 - theta))

    return -1j / 4.0 * series_sum


def greens_neumann(x, y, xi, eta, S_2d, k, d):
    """
    Compute Green's function satisfying Neumann boundary conditions.
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_2d, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2 * d)
    return term1 + term2


def greens_neumann_reg(x, y, xi, eta, S_2d, k, d):
    """
    Compute the Neumann Green's function without the singular 0.25*Y0(k*r) term.
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_2d, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2 * d)
    return term1 + term2


def greens_dirichlet(x, y, xi, eta, S_2d, k, d):
    """
    Compute Green's function satisfying Dirichlet boundary conditions.
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_2d, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2 * d)
    return term1 - term2


def greens_dirichlet_reg(x, y, xi, eta, S_2d, k, d):
    """
    Compute the Dirichlet Green's function without the singular 0.25*Y0(k*r) term.
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_2d, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_2d, k, 2 * d)
    return term1 - term2


def greens_dir_neu(x, y, xi, eta, S_dn, k, d):
    """
    Compute Green's function satisfying Dirichlet at y=0 and Neumann at y=d.
    """
    X = x - xi
    term1 = greens_periodic(X, y - eta, S_dn, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_dn, k, 2 * d)
    return term1 - term2


def greens_dir_neu_reg(x, y, xi, eta, S_dn, k, d):
    """
    Compute the mixed Dirichlet-Neumann Green's function without the singular
    0.25*Y0(k*r) term.
    """
    X = x - xi
    term1 = greens_periodic_reg(X, y - eta, S_dn, k, 2 * d)
    term2 = greens_periodic(X, y + eta, S_dn, k, 2 * d)
    return term1 - term2


if __name__ == "__main__":
    # Original real-k check: values associated with Table 2 of
    # linton1998greens.pdf in the original script.
    d = 1.0
    k = 2.0
    beta = np.sqrt(2.0)
    M = 80
    Lh = 2

    S = lattice_sums(d, k, beta, M, Lh)
    X = 0.0
    Y = 0.01
    G = greens_periodic(X, Y, S, k, d)
    print(f"G_lattice(real k) = {G}")

    # Complex-k smoke test. This is not a reference-value test; it checks that
    # the analytic continuation path evaluates to finite complex values.
    k_complex = 2.0 + 1.0e-2j
    S_complex = lattice_sums(d, k_complex, beta, M, Lh)
    G_complex = greens_periodic(X, Y, S_complex, k_complex, d)

    finite = np.all(np.isfinite(S_complex)) and np.isfinite(G_complex)
    print(f"complex-k smoke test finite = {finite}")
    print(f"G_lattice(complex k) = {G_complex}")
