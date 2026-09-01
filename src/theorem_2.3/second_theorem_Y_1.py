"""
Focused numerical validation of Theorem 2.3(iv): one y-symmetric BIC.

For sufficiently small epsilon, Theorem 2.3(iv) states that if the obstacle
Gamma is symmetric with respect to the y-axis (X odd, Y even), then nu=0.
A unique embedded eigenvalue exists in

    Lambda_1 <= k^2 < Lambda_2

provided the vertical position is tuned as

    a = epsilon * a_1 + O(epsilon^2),

and the eigenvalue satisfies

    k^2 = Lambda_2 - sigma^2,
    sigma = epsilon^2 * pi^3 / b^3 * mu + O(epsilon^3 log epsilon).

For the explicit small-beta obstacle used in the paper examples,

    X0(t) = sin(t) - beta/2 sin(2t),
    Y0(t) = -cos(t) + beta/2 cos(2t),

we use the leading placement formula

    a_1 = -beta/12 + O(beta^2),

so the numerical experiment sets a = epsilon*(-beta/12). This validates the
leading-order explicit y-symmetric example. If an exact/numerically evaluated
a_1 from the theorem integral is available, replace a1_leading() below.

The numerical questions are separated:

1. GEOMETRIC/SYMMETRY CONDITIONS
   - X(-t) = -X(t),
   - Y(-t) =  Y(t),
   - nu is numerically zero,
   - a follows epsilon*a_1 at leading order.

2. EXPECTED BIC BRANCH
   Does the near-singular minimum predicted by the asymptotic formula persist
   and converge as the BEM boundary discretization M is refined?

3. BIC / NON-RADIATION TEST
   Since Lambda_1 < k^2 < Lambda_2, the first transverse channel is open.
   A BIC must have a numerically vanishing projection onto that propagating
   channel on far cross-sections.

4. UNIQUENESS SCREEN
   Does a whole-band scan on Lambda_1 < k^2 < Lambda_2 reveal any additional
   resolved non-radiating near-singular minima away from the expected BIC?

5. ASYMPTOTIC ACCURACY
   How close is sigma_BEM = sqrt(Lambda_2-k_BEM^2) to the leading asymptotic
   prediction sigma_asym?

Important:
- sv_min(A) denotes the smallest singular value of the BEM matrix.
- sigma denotes the spectral distance in Theorem 2.3; they are different.
- The 5% accuracy threshold and other tolerances are numerical reporting
  criteria, not constants appearing in the theorem.
- This is a numerical validation, not a mathematical proof of uniqueness.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from collections.abc import Callable
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

# Same project import convention used by the previous scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

import dipole_theorem as dip
import lattice_sums as lattice

PI = np.pi


@dataclass(frozen=True)
class Config:
    # Guide geometry.
    b: float = 1.0

    # y-axis-symmetric reference shape:
    #   X0(t) = sin(t) - beta/2 sin(2t)       (odd)
    #   Y0(t) = -cos(t) + beta/2 cos(2t)      (even)
    # Physical obstacle: epsilon*(X0,Y0), centered at (0,a(epsilon)).
    shape_beta: float = 0.10

    # Same first experiment used for the x-symmetric case.
    epsilon_values: tuple[float, ...] = (
        0.01, 0.03111111, 0.05222222, 0.07333333, 0.09444444,
        0.11555556, 0.13666667, 0.15777778, 0.17888889, 0.2,
    )

    # BEM / Green function.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Mesh refinement for the expected BIC branch.
    refinement_M: tuple[int, ...] = (16, 24, 32, 40, 48)

    # Expected-mode search window in delta_2 = Lambda_2-k^2 = sigma^2.
    expected_delta_lower_factor: float = 0.20
    expected_delta_upper_factor: float = 5.00
    upper_cutoff_delta_floor: float = 1.0e-9
    minimizer_xatol: float = 2.0e-10

    # Near-singularity diagnostics.
    minimum_drop_factor: float = 100.0
    near_singular_tolerance: float = 1.0e-4
    mesh_sigma_relative_tolerance: float = 0.03

    # BIC diagnostics.
    propagating_ratio_tolerance: float = 1.0e-2
    projection_x_over_b: float = 2.0
    projection_quadrature_points: int = 64

    # Geometry / dipole / placement checks.
    geometry_symmetry_tolerance: float = 1.0e-12
    nu_tolerance: float = 1.0e-8
    placement_tolerance: float = 1.0e-14

    # Leading asymptotic approximation accuracy, reported separately.
    relative_sigma_error_tolerance: float = 0.05

    # Whole-band uniqueness screen.
    uniqueness_scan_M: int = 16
    uniqueness_refine_M: int = 32
    uniqueness_linear_points: int = 36
    uniqueness_log_points: int = 36
    lower_cutoff_kb_margin: float = 2.0e-4
    upper_cutoff_kb_margin: float = 1.0e-7

    # Persistence of additional candidates under two BEM resolutions.
    additional_candidate_kb_tolerance: float = 2.0e-3

    output_directory: str = "theorem_2_3_y_symmetry_validation"


@dataclass
class RefinementRow:
    epsilon: float
    a: float
    M: int
    kb: float
    sigma_bem: float
    sv_min: float
    left_value: float
    right_value: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_previous: float
    yaxis_even_residual: float
    yaxis_odd_residual: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    first_mode_coefficient_left: float
    first_mode_coefficient_right: float
    field_norm_left: float
    field_norm_right: float


@dataclass
class AdditionalCandidate:
    epsilon: float
    a: float
    bracket_left: float
    bracket_right: float
    kb_coarse: float
    sv_min_coarse: float
    kb_fine: float
    sv_min_fine: float
    fine_drop_factor: float
    kb_shift: float
    yaxis_even_residual: float
    yaxis_odd_residual: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    persistent_nonradiating: bool


@dataclass
class ValidationResult:
    epsilon: float
    beta: float
    a1_leading: float
    a: float
    placement_residual: float
    mu: float
    nu: float
    lambda_1: float
    lambda_2: float
    kb_cutoff_1: float
    kb_cutoff_2: float
    kb_asymptotic: float
    kb_numerical: float
    sigma_asymptotic: float
    sigma_numerical: float
    sv_min_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    yaxis_even_residual: float
    yaxis_odd_residual: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    relative_error_kb: float
    relative_error_sigma: float
    symmetry_conditions_verified: bool
    placement_condition_verified: bool
    expected_bic_resolved: bool
    nonradiating_bic_verified: bool
    additional_resolved_bics: int
    uniqueness_screen_passed: bool
    unique_bic_verified: bool
    asymptotic_agreement_verified: bool


CONFIG = Config()


# ---------------------------------------------------------------------------
# Spectral and asymptotic quantities
# ---------------------------------------------------------------------------


def lambda_1(config: Config) -> float:
    return (PI / (2.0 * config.b)) ** 2


def lambda_2(config: Config) -> float:
    return (PI / config.b) ** 2


def kb_cutoff_1(config: Config) -> float:
    return math.sqrt(lambda_1(config)) * config.b


def kb_cutoff_2(config: Config) -> float:
    return math.sqrt(lambda_2(config)) * config.b


@cache
def dipole_moments(shape_beta: float) -> tuple[float, float]:
    """Return (mu,nu) for the unscaled y-symmetric reference obstacle."""
    mu, nu = dip.dipole(shape_beta, 0.0, 0.0, "y")
    return float(mu), float(nu)


def a1_leading(config: Config) -> float:
    """
    Leading explicit formula used in the paper example:

        a_1 = -beta/12 + O(beta^2).

    Replace this function if you later evaluate the full a_1 integral from
    Theorem 2.3(iv).
    """
    return -config.shape_beta / 12.0


def obstacle_position(epsilon: float, config: Config) -> float:
    """Leading placement a = epsilon*a_1."""
    return epsilon * a1_leading(config)


def placement_residual(epsilon: float, config: Config) -> float:
    target = epsilon * a1_leading(config)
    return abs(obstacle_position(epsilon, config) - target)


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    """Return (kb_asymptotic, sigma_asymptotic) for Theorem 2.3(iv)."""
    mu, _ = dipole_moments(config.shape_beta)

    sigma = epsilon**2 * PI**3 / config.b**3 * mu
    if sigma <= 0.0:
        raise ValueError("The leading-order sigma must be positive.")

    k_squared = lambda_2(config) - sigma**2
    if not (lambda_1(config) <= k_squared < lambda_2(config)):
        raise ValueError(
            "The asymptotic prediction lies outside [Lambda_1,Lambda_2). "
            "Use a smaller epsilon or inspect the dipole moment."
        )

    return math.sqrt(k_squared) * config.b, sigma


def sigma_from_kb(kb: float, config: Config) -> float:
    """Theorem sigma = sqrt(Lambda_2-k^2)."""
    k = kb / config.b
    gap = lambda_2(config) - k**2
    return math.sqrt(max(gap, 0.0))


def kb_from_delta_2(delta_2: float, config: Config) -> float:
    """delta_2 = Lambda_2-k^2 = sigma^2; return dimensionless kb."""
    lam1 = lambda_1(config)
    lam2 = lambda_2(config)
    max_gap = lam2 - lam1
    delta_2 = min(
        max(delta_2, config.upper_cutoff_delta_floor),
        max_gap * (1.0 - 1.0e-10),
    )
    return config.b * math.sqrt(lam2 - delta_2)


# ---------------------------------------------------------------------------
# Geometry and symmetry
# ---------------------------------------------------------------------------


def obstacle_geometry(
    t: np.ndarray | float,
    epsilon: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    y-axis-symmetric reference obstacle in obstacle-centred coordinates:

      X(t) = eps [sin t - beta/2 sin(2t)]       (odd)
      Y(t) = eps [-cos t + beta/2 cos(2t)]      (even)

    The physical vertical coordinate is Y(t)+a(epsilon).
    """
    t = np.asarray(t)
    beta = config.shape_beta

    X = epsilon * (np.sin(t) - 0.5 * beta * np.sin(2.0 * t))
    Y = epsilon * (-np.cos(t) + 0.5 * beta * np.cos(2.0 * t))

    Xp = epsilon * (np.cos(t) - beta * np.cos(2.0 * t))
    Yp = epsilon * (np.sin(t) - beta * np.sin(2.0 * t))

    Xpp = epsilon * (-np.sin(t) + 2.0 * beta * np.sin(2.0 * t))
    Ypp = epsilon * (np.cos(t) - 2.0 * beta * np.cos(2.0 * t))

    return X, Y, Xp, Yp, Xpp, Ypp


def geometry_symmetry_residual(config: Config) -> tuple[float, float]:
    """Return max residuals for X(-t)=-X(t) and Y(-t)=Y(t)."""
    t = np.linspace(-PI, PI, 257)
    Xp, Yp, *_ = obstacle_geometry(t, 1.0, config)
    Xm, Ym, *_ = obstacle_geometry(-t, 1.0, config)

    x_odd = float(np.max(np.abs(Xm + Xp)))
    y_even = float(np.max(np.abs(Ym - Yp)))
    return x_odd, y_even


def obstacle_fits_guide(epsilon: float, config: Config) -> bool:
    """Check that the vertically shifted obstacle remains inside y in (-b,b)."""
    t = np.linspace(0.0, 2.0 * PI, 1025)
    _, Y, *_ = obstacle_geometry(t, epsilon, config)
    a = obstacle_position(epsilon, config)
    Y_global = Y + a
    return bool(np.max(np.abs(Y_global)) < config.b)


# ---------------------------------------------------------------------------
# BEM assembly
# ---------------------------------------------------------------------------


def boundary_nodes(M: int) -> np.ndarray:
    return (np.arange(M) + 0.5) * (2.0 * PI / M)


def make_green_functions(
    kb: float,
    epsilon: float,
    config: Config,
) -> tuple[Callable[..., complex], Callable[..., complex]]:
    b = config.b
    d = 2.0 * b
    k = kb / b
    a = obstacle_position(epsilon, config)

    lattice_coefficients = lattice.lattice_sums(
        2.0 * d,
        k,
        beta=0.0,
        M=config.lattice_terms,
        Lh=config.harmonic_order,
    )

    def green(x: float, y: float, xi: float, eta: float) -> complex:
        return lattice.greens_dirichlet(
            x,
            y + b + a,
            xi,
            eta + b + a,
            lattice_coefficients,
            k,
            d,
        )

    def green_regularized(x: float, y: float, xi: float, eta: float) -> complex:
        return lattice.greens_dirichlet_reg(
            x,
            y + b + a,
            xi,
            eta + b + a,
            lattice_coefficients,
            k,
            d,
        )

    return green, green_regularized


def source_derivatives(
    G: Callable[..., complex],
    x: float,
    y: float,
    xi: float,
    eta: float,
    h: float,
) -> tuple[complex, complex]:
    dG_dxi = (G(x, y, xi + h, eta) - G(x, y, xi - h, eta)) / (2.0 * h)
    dG_deta = (G(x, y, xi, eta + h) - G(x, y, xi, eta - h)) / (2.0 * h)
    return dG_dxi, dG_deta


def weighted_normal_kernel(
    psi: float,
    theta: float,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
    G_regularized: Callable[..., complex],
) -> complex:
    x, y, _, _, _, _ = obstacle_geometry(psi, epsilon, config)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = obstacle_geometry(theta, epsilon, config)

    w = float(np.hypot(xi_p, eta_p))
    h = config.finite_difference_step

    if abs(psi - theta) > 1.0e-12:
        G_xi, G_eta = source_derivatives(G, x, y, xi, eta, h)
        return xi_p * G_eta - eta_p * G_xi

    G_xi_reg, G_eta_reg = source_derivatives(G_regularized, x, y, xi, eta, h)
    geometric_term = (xi_pp * eta_p - eta_pp * xi_p) / (4.0 * PI * w**2)
    regularized_term = xi_p * G_eta_reg - eta_p * G_xi_reg
    return geometric_term + regularized_term


def assemble_matrix(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, Callable[..., complex]]:
    G, G_regularized = make_green_functions(kb, epsilon, config)
    theta = boundary_nodes(M)

    K_weighted = np.empty((M, M), dtype=np.complex128)
    for i, psi in enumerate(theta):
        for j, source_theta in enumerate(theta):
            K_weighted[i, j] = weighted_normal_kernel(
                float(psi),
                float(source_theta),
                epsilon,
                config,
                G,
                G_regularized,
            )

    # Full contour [0,2pi]:
    #   1/2 u_i = (2pi/M) sum_j K^w_ij u_j
    # => A = I - (4pi/M) K^w.
    A = np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted
    return A, theta, G


def smallest_singular_pair(A: np.ndarray) -> tuple[float, np.ndarray]:
    _, singular_values, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= max(np.linalg.norm(vector), 1.0e-30)

    # Fix an arbitrary global phase for stable diagnostics/plots.
    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))

    return float(singular_values[-1]), vector


def smallest_singular_value(kb: float, epsilon: float, M: int, config: Config) -> float:
    A, _, _ = assemble_matrix(kb, epsilon, M, config)
    singular_values = np.linalg.svd(A, compute_uv=False)
    return float(singular_values[-1])


def spectral_state(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, np.ndarray, np.ndarray, Callable[..., complex]]:
    A, theta, G = assemble_matrix(kb, epsilon, M, config)
    sv_min, vector = smallest_singular_pair(A)
    return sv_min, vector, theta, G


# ---------------------------------------------------------------------------
# BIC diagnostics: symmetry and cancellation of the open channel
# ---------------------------------------------------------------------------


def boundary_parity_residuals(boundary_vector: np.ndarray) -> tuple[float, float]:
    """
    Reflection across the y-axis maps (x,y)->(-x,y), which corresponds to
    t -> -t (mod 2pi). With midpoint nodes this maps j -> M-1-j.

    The theorem does not require a specific eigenfunction parity under this
    reflection, so both even and odd residuals are reported as diagnostics;
    neither is used as a BIC acceptance condition.
    """
    reflected = boundary_vector[::-1]
    norm = max(float(np.linalg.norm(boundary_vector)), 1.0e-30)

    even_residual = float(np.linalg.norm(boundary_vector - reflected) / norm)
    odd_residual = float(np.linalg.norm(boundary_vector + reflected) / norm)
    return even_residual, odd_residual


def weighted_kernel_at_field_point(
    x: float,
    y: float,
    theta: float,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
) -> complex:
    xi, eta, xi_p, eta_p, _, _ = obstacle_geometry(theta, epsilon, config)
    G_xi, G_eta = source_derivatives(
        G,
        x,
        y,
        float(xi),
        float(eta),
        config.finite_difference_step,
    )
    return float(xi_p) * G_eta - float(eta_p) * G_xi


def reconstruct_field(
    x: float,
    y_values: np.ndarray,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
) -> np.ndarray:
    delta_theta = 2.0 * PI / len(theta)
    field = np.empty(len(y_values), dtype=np.complex128)

    for i, y in enumerate(y_values):
        kernel_values = np.array(
            [
                weighted_kernel_at_field_point(
                    x,
                    float(y),
                    float(t),
                    epsilon,
                    config,
                    G,
                )
                for t in theta
            ],
            dtype=np.complex128,
        )
        field[i] = delta_theta * np.dot(kernel_values, boundary_vector)

    return field


def first_propagating_component(
    epsilon: float,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    G: Callable[..., complex],
    config: Config,
) -> dict[str, float]:
    """
    Project the reconstructed field onto the only open transverse mode in
    Lambda_1 < k^2 < Lambda_2.

    Geometry is stored in obstacle-centred y coordinates. The physical guide
    coordinate is y_global = y + a(epsilon), so

        phi_1 = cos(pi*y_global/(2b)).

    For a BIC this propagating-channel coefficient must vanish numerically.
    """
    b = config.b
    a = obstacle_position(epsilon, config)
    x_far = config.projection_x_over_b * b

    nodes, weights = leggauss(config.projection_quadrature_points)

    # In obstacle-centred coordinates the guide walls are at
    # y=-b-a and y=b-a.
    y_lower = -b - a
    y_upper = b - a
    y_values = 0.5 * (y_upper - y_lower) * nodes + 0.5 * (y_upper + y_lower)
    physical_weights = 0.5 * (y_upper - y_lower) * weights

    y_global = y_values + a
    phi_1 = np.cos(PI * y_global / (2.0 * b))
    norm_phi_squared = float(np.sum(physical_weights * np.abs(phi_1) ** 2))

    diagnostics: dict[str, float] = {}
    for label, x in (("left", -x_far), ("right", x_far)):
        field = reconstruct_field(
            x,
            y_values,
            boundary_vector,
            theta,
            epsilon,
            config,
            G,
        )

        coefficient = (
            np.sum(physical_weights * field * np.conjugate(phi_1)) / norm_phi_squared
        )

        projected_norm = abs(coefficient) * math.sqrt(norm_phi_squared)
        field_norm = float(math.sqrt(np.sum(physical_weights * np.abs(field) ** 2)))
        ratio = float(projected_norm / max(field_norm, 1.0e-30))

        diagnostics[f"ratio_{label}"] = ratio
        diagnostics[f"coefficient_{label}"] = float(abs(coefficient))
        diagnostics[f"field_norm_{label}"] = field_norm

    return diagnostics


# ---------------------------------------------------------------------------
# Expected BIC branch: targeted mesh refinement
# ---------------------------------------------------------------------------


def expected_mode_bracket(epsilon: float, config: Config) -> tuple[float, float]:
    """Bracket the expected BIC in delta_2 = Lambda_2-k^2 coordinates."""
    _, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2
    max_gap = lambda_2(config) - lambda_1(config)

    delta_far = min(delta_asym * config.expected_delta_upper_factor, max_gap * 0.95)
    delta_near = max(
        delta_asym * config.expected_delta_lower_factor,
        config.upper_cutoff_delta_floor,
    )

    # Larger delta means smaller k.
    left = kb_from_delta_2(delta_far, config)
    right = kb_from_delta_2(delta_near, config)

    if not left < right:
        raise RuntimeError("Invalid expected-BIC search bracket.")
    return left, right


def refine_expected_mode_for_M(
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    bool,
    float,
    float,
    dict[str, float],
]:
    """
    Return
      kb, sv_min, left_value, right_value, drop_factor, interior,
      yaxis_even_residual, yaxis_odd_residual, propagation_diagnostics.
    """
    left, right = expected_mode_bracket(epsilon, config)
    tiny = np.finfo(float).tiny

    left_value = smallest_singular_value(left, epsilon, M, config)
    right_value = smallest_singular_value(right, epsilon, M, config)

    result = minimize_scalar(
        lambda kb: math.log10(
            max(smallest_singular_value(float(kb), epsilon, M, config), tiny)
        ),
        bounds=(left, right),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    if not result.success:
        raise RuntimeError(
            f"Expected-BIC minimization failed for M={M}: {result.message}"
        )

    kb = float(result.x)
    sv_min, boundary_vector, theta, G = spectral_state(kb, epsilon, M, config)

    endpoint_level = min(left_value, right_value)
    drop_factor = endpoint_level / max(sv_min, tiny)

    width = right - left
    edge_margin = 0.01 * width
    interior = (kb > left + edge_margin) and (kb < right - edge_margin)

    even_residual, odd_residual = boundary_parity_residuals(boundary_vector)
    propagation = first_propagating_component(
        epsilon,
        boundary_vector,
        theta,
        G,
        config,
    )

    return (
        kb,
        sv_min,
        left_value,
        right_value,
        drop_factor,
        interior,
        even_residual,
        odd_residual,
        propagation,
    )


def run_expected_branch_refinement(
    epsilon: float,
    config: Config,
) -> list[RefinementRow]:
    rows: list[RefinementRow] = []
    previous_sigma = math.nan
    a = obstacle_position(epsilon, config)

    print("  expected-BIC mesh refinement:")
    for M in config.refinement_M:
        (
            kb,
            sv_min,
            left_value,
            right_value,
            drop_factor,
            interior,
            even_residual,
            odd_residual,
            propagation,
        ) = refine_expected_mode_for_M(epsilon, M, config)

        sigma_bem = sigma_from_kb(kb, config)
        change = (
            abs(sigma_bem - previous_sigma) / max(abs(sigma_bem), 1.0e-30)
            if np.isfinite(previous_sigma)
            else math.nan
        )

        row = RefinementRow(
            epsilon=epsilon,
            a=a,
            M=M,
            kb=kb,
            sigma_bem=sigma_bem,
            sv_min=sv_min,
            left_value=left_value,
            right_value=right_value,
            drop_factor=drop_factor,
            minimum_is_interior=interior,
            relative_sigma_change_from_previous=change,
            yaxis_even_residual=even_residual,
            yaxis_odd_residual=odd_residual,
            propagating_ratio_left=propagation["ratio_left"],
            propagating_ratio_right=propagation["ratio_right"],
            first_mode_coefficient_left=propagation["coefficient_left"],
            first_mode_coefficient_right=propagation["coefficient_right"],
            field_norm_left=propagation["field_norm_left"],
            field_norm_right=propagation["field_norm_right"],
        )
        rows.append(row)
        previous_sigma = sigma_bem

        change_text = f"{change:.3%}" if np.isfinite(change) else "--"
        max_prop = max(row.propagating_ratio_left, row.propagating_ratio_right)
        parity_label = "even" if even_residual <= odd_residual else "odd"
        parity_res = min(even_residual, odd_residual)
        print(
            f"    M={M:>2}: kb={kb:.12f}, "
            f"sigma_BEM={sigma_bem:.8e}, sv_min={sv_min:.3e}, "
            f"drop={drop_factor:.2e}, mesh_change={change_text}, "
            f"y-axis parity={parity_label}({parity_res:.3e}), "
            f"Rprop={max_prop:.3e}, "
            f"interior={'yes' if interior else 'no'}"
        )

    return rows


def expected_bic_is_resolved(
    rows: list[RefinementRow], config: Config
) -> tuple[bool, bool]:
    """
    Return (spectral_mode_resolved, nonradiating_bic_verified).

    Unlike statement (iii), no specific y-parity cancellation is imposed here.
    For statement (iv), the tuned vertical placement is essential and the BIC
    is verified directly by near-singularity, mesh persistence, and vanishing
    projection onto the open first channel.
    """
    if len(rows) < 2:
        return False, False

    final = rows[-1]
    previous = rows[-2]

    final_change = abs(final.sigma_bem - previous.sigma_bem) / max(
        abs(final.sigma_bem), 1.0e-30
    )

    near_singular = final.sv_min <= config.near_singular_tolerance
    deep_minimum = final.drop_factor >= config.minimum_drop_factor
    converged = final_change <= config.mesh_sigma_relative_tolerance

    spectral_resolved = bool(
        final.minimum_is_interior and near_singular and deep_minimum and converged
    )

    propagation_ok = (
        final.propagating_ratio_left <= config.propagating_ratio_tolerance
        and final.propagating_ratio_right <= config.propagating_ratio_tolerance
    )

    nonradiating_bic = bool(spectral_resolved and propagation_ok)
    return spectral_resolved, nonradiating_bic


# ---------------------------------------------------------------------------
# Whole-band uniqueness screen on Lambda_1 < k^2 < Lambda_2
# ---------------------------------------------------------------------------


def build_uniqueness_grid(epsilon: float, config: Config) -> np.ndarray:
    kb1 = kb_cutoff_1(config)
    kb2 = kb_cutoff_2(config)
    lam1 = lambda_1(config)
    lam2 = lambda_2(config)

    _, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2

    linear = np.linspace(
        kb1 + config.lower_cutoff_kb_margin,
        kb2 - config.upper_cutoff_kb_margin,
        config.uniqueness_linear_points,
    )

    # Logarithmic resolution in delta_2 close to Lambda_2, where the BIC is
    # asymptotically located.
    delta_max = (lam2 - lam1) * (1.0 - 1.0e-6)
    deltas = np.geomspace(
        config.upper_cutoff_delta_floor,
        delta_max,
        config.uniqueness_log_points,
    )
    logarithmic = config.b * np.sqrt(lam2 - deltas)

    local_factors = np.array([0.15, 0.25, 0.5, 0.75, 1.0, 1.5, 2.5, 5.0])
    local_delta = delta_asym * local_factors
    local_delta = local_delta[
        (local_delta > config.upper_cutoff_delta_floor) & (local_delta < (lam2 - lam1))
    ]
    local = config.b * np.sqrt(lam2 - local_delta)

    grid = np.concatenate([linear, logarithmic, local])
    grid = grid[(grid > kb1 + 0.5 * config.lower_cutoff_kb_margin) & (grid < kb2)]
    return np.unique(np.sort(grid))


def sampled_local_minimum_brackets(
    scan: np.ndarray,
    values: np.ndarray,
) -> list[tuple[float, float]]:
    brackets: list[tuple[float, float]] = []
    for i in range(1, len(scan) - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            brackets.append((float(scan[i - 1]), float(scan[i + 1])))
    return brackets


def bracket_overlaps_expected_mode(
    left: float,
    right: float,
    epsilon: float,
    config: Config,
) -> bool:
    expected_left, expected_right = expected_mode_bracket(epsilon, config)
    return max(left, expected_left) <= min(right, expected_right)


def refine_bracket_once(
    left: float,
    right: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, np.ndarray, np.ndarray, Callable[..., complex]]:
    tiny = np.finfo(float).tiny
    left_value = smallest_singular_value(left, epsilon, M, config)
    right_value = smallest_singular_value(right, epsilon, M, config)

    result = minimize_scalar(
        lambda kb: math.log10(
            max(smallest_singular_value(float(kb), epsilon, M, config), tiny)
        ),
        bounds=(left, right),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    if not result.success:
        return math.nan, math.nan, math.nan, np.array([]), np.array([]), None  # type: ignore[return-value]

    kb = float(result.x)
    sv_min, vector, theta, G = spectral_state(kb, epsilon, M, config)
    drop_factor = min(left_value, right_value) / max(sv_min, tiny)
    return kb, sv_min, drop_factor, vector, theta, G


def whole_band_uniqueness_screen(
    epsilon: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, list[AdditionalCandidate]]:
    scan = build_uniqueness_grid(epsilon, config)
    values = np.empty(len(scan), dtype=float)
    a = obstacle_position(epsilon, config)

    print(
        f"  whole-band uniqueness screen: {len(scan)} points "
        f"(M={config.uniqueness_scan_M})"
    )
    for i, kb in enumerate(scan):
        values[i] = smallest_singular_value(
            float(kb), epsilon, config.uniqueness_scan_M, config
        )
        if (i + 1) % 16 == 0 or i + 1 == len(scan):
            print(f"    {i + 1:>3}/{len(scan)}")

    brackets = sampled_local_minimum_brackets(scan, values)
    other_brackets = [
        bracket
        for bracket in brackets
        if not bracket_overlaps_expected_mode(bracket[0], bracket[1], epsilon, config)
    ]

    print(
        f"  sampled local minima: {len(brackets)} total; "
        f"{len(other_brackets)} outside expected-BIC window"
    )

    additional: list[AdditionalCandidate] = []
    for j, (left, right) in enumerate(other_brackets, start=1):
        # Cheap rejection: if the sampled center is nowhere near singular and
        # does not form a visible drop, avoid expensive field reconstruction.
        print(f"    checking additional minimum {j}/{len(other_brackets)}")

        kb_c, sv_c, _, _, _, _ = refine_bracket_once(
            left, right, epsilon, config.uniqueness_scan_M, config
        )
        kb_f, sv_f, drop_f, vector_f, theta_f, G_f = refine_bracket_once(
            left, right, epsilon, config.uniqueness_refine_M, config
        )

        kb_shift = (
            abs(kb_f - kb_c) if np.isfinite(kb_c) and np.isfinite(kb_f) else math.nan
        )

        if np.isfinite(kb_f) and len(vector_f) > 0:
            even_residual, odd_residual = boundary_parity_residuals(vector_f)
            propagation = first_propagating_component(
                epsilon, vector_f, theta_f, G_f, config
            )
            ratio_left = propagation["ratio_left"]
            ratio_right = propagation["ratio_right"]
        else:
            even_residual = math.nan
            odd_residual = math.nan
            ratio_left = math.nan
            ratio_right = math.nan

        persistent_nonradiating = bool(
            np.isfinite(kb_f)
            and np.isfinite(sv_f)
            and sv_f <= config.near_singular_tolerance
            and drop_f >= config.minimum_drop_factor
            and np.isfinite(kb_shift)
            and kb_shift <= config.additional_candidate_kb_tolerance
            and np.isfinite(ratio_left)
            and np.isfinite(ratio_right)
            and ratio_left <= config.propagating_ratio_tolerance
            and ratio_right <= config.propagating_ratio_tolerance
        )

        additional.append(
            AdditionalCandidate(
                epsilon=epsilon,
                a=a,
                bracket_left=left,
                bracket_right=right,
                kb_coarse=kb_c,
                sv_min_coarse=sv_c,
                kb_fine=kb_f,
                sv_min_fine=sv_f,
                fine_drop_factor=drop_f,
                kb_shift=kb_shift,
                yaxis_even_residual=even_residual,
                yaxis_odd_residual=odd_residual,
                propagating_ratio_left=ratio_left,
                propagating_ratio_right=ratio_right,
                persistent_nonradiating=persistent_nonradiating,
            )
        )

        max_prop = max(ratio_left, ratio_right) if np.isfinite(ratio_left) else math.nan
        parity_res = (
            min(even_residual, odd_residual) if np.isfinite(even_residual) else math.nan
        )
        parity_label = (
            "even"
            if np.isfinite(even_residual) and even_residual <= odd_residual
            else "odd"
        )
        print(f"      kb(M={config.uniqueness_scan_M})={kb_c:.12f}, sv_min={sv_c:.3e}")
        print(
            f"      kb(M={config.uniqueness_refine_M})={kb_f:.12f}, "
            f"sv_min={sv_f:.3e}, drop={drop_f:.2e}, "
            f"shift={kb_shift:.3e}, y-axis parity={parity_label}({parity_res:.3e}), "
            f"Rprop={max_prop:.3e}, "
            f"persistent BIC={'YES' if persistent_nonradiating else 'no'}"
        )

    return scan, values, additional


# ---------------------------------------------------------------------------
# Validation and reporting
# ---------------------------------------------------------------------------


def symmetry_conditions_verified(
    config: Config,
) -> tuple[bool, float, float, float, float]:
    x_odd_residual, y_even_residual = geometry_symmetry_residual(config)
    mu, nu = dipole_moments(config.shape_beta)

    passed = bool(
        x_odd_residual <= config.geometry_symmetry_tolerance
        and y_even_residual <= config.geometry_symmetry_tolerance
        and abs(nu) <= config.nu_tolerance
    )
    return passed, x_odd_residual, y_even_residual, mu, nu


def placement_condition_verified(epsilon: float, config: Config) -> tuple[bool, float]:
    residual = placement_residual(epsilon, config)
    passed = bool(
        residual <= config.placement_tolerance and obstacle_fits_guide(epsilon, config)
    )
    return passed, residual


def validate_epsilon(
    epsilon: float,
    config: Config,
) -> tuple[
    ValidationResult,
    list[RefinementRow],
    np.ndarray,
    np.ndarray,
    list[AdditionalCandidate],
]:
    symmetry_ok, _, _, mu, nu = symmetry_conditions_verified(config)
    placement_ok, placement_err = placement_condition_verified(epsilon, config)
    a = obstacle_position(epsilon, config)
    a1 = a1_leading(config)
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)

    refinement = run_expected_branch_refinement(epsilon, config)
    spectral_resolved, nonradiating_bic = expected_bic_is_resolved(refinement, config)

    final = refinement[-1]
    kb_num = final.kb
    sigma_num = final.sigma_bem
    final_mesh_change = final.relative_sigma_change_from_previous

    error_kb = abs(kb_num - kb_asym) / abs(kb_asym)
    error_sigma = abs(sigma_num - sigma_asym) / abs(sigma_asym)
    asymptotic_ok = error_sigma <= config.relative_sigma_error_tolerance

    scan, scan_values, additional = whole_band_uniqueness_screen(epsilon, config)
    additional_resolved = [
        candidate for candidate in additional if candidate.persistent_nonradiating
    ]
    uniqueness_screen_passed = len(additional_resolved) == 0

    unique_bic_verified = bool(
        symmetry_ok and placement_ok and nonradiating_bic and uniqueness_screen_passed
    )

    result = ValidationResult(
        epsilon=epsilon,
        beta=config.shape_beta,
        a1_leading=a1,
        a=a,
        placement_residual=placement_err,
        mu=mu,
        nu=nu,
        lambda_1=lambda_1(config),
        lambda_2=lambda_2(config),
        kb_cutoff_1=kb_cutoff_1(config),
        kb_cutoff_2=kb_cutoff_2(config),
        kb_asymptotic=kb_asym,
        kb_numerical=kb_num,
        sigma_asymptotic=sigma_asym,
        sigma_numerical=sigma_num,
        sv_min_final=final.sv_min,
        final_drop_factor=final.drop_factor,
        final_relative_mesh_change=final_mesh_change,
        yaxis_even_residual=final.yaxis_even_residual,
        yaxis_odd_residual=final.yaxis_odd_residual,
        propagating_ratio_left=final.propagating_ratio_left,
        propagating_ratio_right=final.propagating_ratio_right,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        symmetry_conditions_verified=symmetry_ok,
        placement_condition_verified=placement_ok,
        expected_bic_resolved=spectral_resolved,
        nonradiating_bic_verified=nonradiating_bic,
        additional_resolved_bics=len(additional_resolved),
        uniqueness_screen_passed=uniqueness_screen_passed,
        unique_bic_verified=unique_bic_verified,
        asymptotic_agreement_verified=asymptotic_ok,
    )

    return result, refinement, scan, scan_values, additional


def write_dataclass_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def plot_full_band_scan(
    epsilon: float,
    scan: np.ndarray,
    values: np.ndarray,
    result: ValidationResult,
    output_directory: Path,
) -> None:
    plt.figure(figsize=(8, 5))
    plt.semilogy(scan, values, "o-", markersize=3, label=r"coarse $s_{\min}(A)$")
    plt.axvline(result.kb_asymptotic, linestyle="--", label=r"$kb_{\mathrm{asym}}$")
    plt.axvline(result.kb_numerical, linestyle=":", label=r"$kb_{\mathrm{BEM}}$")
    plt.xlabel(r"$kb$")
    plt.ylabel(r"$s_{\min}(A(k))$")
    plt.title(rf"Theorem 2.3(iv), $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"full_band_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_refinement(
    epsilon: float,
    rows: list[RefinementRow],
    sigma_asym: float,
    output_directory: Path,
    config: Config,
) -> None:
    M = np.array([row.M for row in rows])
    sigma_bem = np.array([row.sigma_bem for row in rows])
    sv_min = np.array([row.sv_min for row in rows])
    best_parity = np.array(
        [min(row.yaxis_even_residual, row.yaxis_odd_residual) for row in rows]
    )
    prop = np.array(
        [max(row.propagating_ratio_left, row.propagating_ratio_right) for row in rows]
    )

    plt.figure(figsize=(7, 4.5))
    plt.plot(M, sigma_bem, "o-", label=r"$\sigma_{\mathrm{BEM}}$")
    plt.axhline(sigma_asym, linestyle="--", label=r"$\sigma_{\mathrm{asym}}$")
    plt.xlabel(r"$M$")
    plt.ylabel(r"spectral distance $\sigma$")
    plt.title(rf"BIC mesh refinement, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        output_directory / f"refinement_sigma_epsilon_{epsilon:.2f}.png", dpi=180
    )
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(M, sv_min, "o-")
    plt.axhline(
        config.near_singular_tolerance, linestyle="--", label="diagnostic tolerance"
    )
    plt.xlabel(r"$M$")
    plt.ylabel(r"$s_{\min}(A(k_*))$")
    plt.title(rf"Near-singularity, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"refinement_svd_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(M, prop, "o-")
    plt.axhline(
        config.propagating_ratio_tolerance,
        linestyle="--",
        label="BIC radiation tolerance",
    )
    plt.xlabel(r"$M$")
    plt.ylabel(r"first propagating-mode fraction")
    plt.title(rf"Open-channel cancellation, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        output_directory / f"refinement_propagation_epsilon_{epsilon:.2f}.png", dpi=180
    )
    plt.close()

    # Pure diagnostic: for a simple eigenstate of the y-reflection-symmetric
    # geometry one parity residual is usually small, but Theorem 2.3(iv) does
    # not use a prescribed parity as the BIC mechanism.
    plt.figure(figsize=(7, 4.5))
    plt.semilogy(M, best_parity, "o-")
    plt.xlabel(r"$M$")
    plt.ylabel(r"best y-axis reflection-parity residual")
    plt.title(rf"y-axis reflection diagnostic, $\varepsilon={epsilon:.2f}$")
    plt.tight_layout()
    plt.savefig(
        output_directory / f"refinement_parity_epsilon_{epsilon:.2f}.png", dpi=180
    )
    plt.close()


def plot_summary(
    results: list[ValidationResult],
    output_directory: Path,
    config: Config,
) -> None:
    eps = np.array([row.epsilon for row in results])
    sigma_asym = np.array([row.sigma_asymptotic for row in results])
    sigma_num = np.array([row.sigma_numerical for row in results])
    kb_asym = np.array([row.kb_asymptotic for row in results])
    kb_num = np.array([row.kb_numerical for row in results])
    error = np.array([row.relative_error_sigma for row in results])
    prop = np.array(
        [
            max(row.propagating_ratio_left, row.propagating_ratio_right)
            for row in results
        ]
    )

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, kb_num, "o-", label=r"$kb_{\mathrm{BEM}}$")
    plt.plot(eps, kb_asym, "s--", label=r"$kb_{\mathrm{asym}}$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$kb$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "kb_vs_epsilon.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, sigma_num, "o-", label=r"$\sigma_{\mathrm{BEM}}$")
    plt.plot(eps, sigma_asym, "s--", label=r"$\sigma_{\mathrm{asym}}$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$\sigma$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "sigma_bem_vs_asymptotic.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, error, "o-")
    plt.axhline(
        config.relative_sigma_error_tolerance,
        linestyle="--",
        label="5% tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"relative error in $\sigma$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "relative_error_sigma.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(eps, prop, "o-")
    plt.axhline(
        config.propagating_ratio_tolerance,
        linestyle="--",
        label="BIC radiation tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"first propagating-mode fraction")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "propagating_fraction.png", dpi=180)
    plt.close()


def print_result(result: ValidationResult, config: Config) -> None:
    max_prop = max(result.propagating_ratio_left, result.propagating_ratio_right)
    if result.yaxis_even_residual <= result.yaxis_odd_residual:
        parity_label = "even"
        parity_res = result.yaxis_even_residual
    else:
        parity_label = "odd"
        parity_res = result.yaxis_odd_residual

    print("\n  --- validation result ---")
    print(
        f"  shape symmetry (X odd, Y even, nu~0): "
        f"{'PASS' if result.symmetry_conditions_verified else 'FAIL'}"
    )
    print(
        f"  placement a=epsilon*a1 (leading order): "
        f"{'PASS' if result.placement_condition_verified else 'FAIL'} "
        f"[a1={result.a1_leading:.8e}, a={result.a:.8e}]"
    )
    print(
        f"  expected spectral minimum resolved: "
        f"{'PASS' if result.expected_bic_resolved else 'FAIL'}"
    )
    print(
        f"  non-radiating BIC diagnostics: "
        f"{'PASS' if result.nonradiating_bic_verified else 'FAIL'}"
    )
    print(f"  additional resolved BICs: {result.additional_resolved_bics}")
    print(
        f"  uniqueness screen: {'PASS' if result.uniqueness_screen_passed else 'FAIL'}"
    )
    print(
        f"  exactly one resolved BIC in [Lambda_1,Lambda_2): "
        f"{'PASS' if result.unique_bic_verified else 'FAIL'}"
    )

    print(f"  kb asymptotic = {result.kb_asymptotic:.12f}")
    print(f"  kb BEM        = {result.kb_numerical:.12f}")
    print(f"  sigma asym    = {result.sigma_asymptotic:.8e}")
    print(f"  sigma BEM     = {result.sigma_numerical:.8e}")
    print(f"  final sv_min(A) = {result.sv_min_final:.3e}")
    print(f"  final minimum drop = {result.final_drop_factor:.2e}")
    print(f"  final mesh change in sigma = {result.final_relative_mesh_change:.3%}")
    print(f"  y-axis parity diagnostic = {parity_label} residual {parity_res:.3e}")
    print(f"  max first propagating-mode fraction = {max_prop:.3e}")
    print(f"  relative sigma error = {result.relative_error_sigma:.3%}")
    print(
        f"  asymptotic accuracy <= {config.relative_sigma_error_tolerance:.1%}: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )

    print(
        f"  EXISTENCE / BIC / UNIQUENESS CHECK: "
        f"{'PASS' if result.unique_bic_verified else 'FAIL'}"
    )
    print(
        f"  ASYMPTOTIC APPROXIMATION CHECK: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    symmetry_ok, x_odd_res, y_even_res, mu, nu = symmetry_conditions_verified(config)

    print("=== Focused validation of Theorem 2.3(iv): y-symmetric BIC ===")
    print(f"b = {config.b}")
    print(f"beta = {config.shape_beta}")
    print(f"leading a1 = -beta/12 = {a1_leading(config):.12e}")
    print(f"mu = {mu:.12e}")
    print(f"nu = {nu:.12e}")
    print(f"X-odd residual  = {x_odd_res:.3e}")
    print(f"Y-even residual = {y_even_res:.3e}")
    print(f"Theorem shape symmetry conditions: {'PASS' if symmetry_ok else 'FAIL'}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"Lambda_2 = {lambda_2(config):.12f}")
    print(f"sqrt(Lambda_1) b = {kb_cutoff_1(config):.12f}")
    print(f"sqrt(Lambda_2) b = {kb_cutoff_2(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")
    print(f"refinement M = {config.refinement_M}")

    if not symmetry_ok:
        raise RuntimeError(
            "The selected obstacle does not satisfy the y-axis symmetry / nu=0 checks."
        )

    summary: list[ValidationResult] = []
    all_refinement: list[RefinementRow] = []
    all_additional: list[AdditionalCandidate] = []

    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.3f} ===")
        a = obstacle_position(epsilon, config)
        placement_ok, placement_err = placement_condition_verified(epsilon, config)
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)

        print(f"  a = epsilon*a1 = {a:.12e}")
        print(
            f"  placement condition: {'PASS' if placement_ok else 'FAIL'} "
            f"(residual={placement_err:.3e})"
        )
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted sigma = {sigma_asym:.8e}")
        print(
            f"  predicted kb gap to sqrt(Lambda_2) = {kb_cutoff_2(config) - kb_asym:.3e}"
        )

        result, refinement, scan, scan_values, additional = validate_epsilon(
            epsilon, config
        )

        summary.append(result)
        all_refinement.extend(refinement)
        all_additional.extend(additional)

        print_result(result, config)
        plot_full_band_scan(epsilon, scan, scan_values, result, output_directory)
        plot_refinement(epsilon, refinement, sigma_asym, output_directory, config)

    write_dataclass_csv(output_directory / "summary.csv", summary)
    write_dataclass_csv(output_directory / "mesh_refinement.csv", all_refinement)
    write_dataclass_csv(output_directory / "additional_candidates.csv", all_additional)
    plot_summary(summary, output_directory, config)

    unique_passed = sum(row.unique_bic_verified for row in summary)
    asymptotic_passed = sum(row.asymptotic_agreement_verified for row in summary)

    print("\n=== FINAL SUMMARY ===")
    print(
        "Interval screened: Lambda_1 < k^2 < Lambda_2 "
        "(excluding tiny numerical layers at the cutoffs)."
    )
    print(f"Unique non-radiating BIC checks passed: {unique_passed}/{len(summary)}")
    print(
        f"Leading asymptotic sigma within "
        f"{config.relative_sigma_error_tolerance:.1%}: "
        f"{asymptotic_passed}/{len(summary)}"
    )

    for row in summary:
        max_prop = max(row.propagating_ratio_left, row.propagating_ratio_right)
        parity_res = min(row.yaxis_even_residual, row.yaxis_odd_residual)
        parity_label = (
            "even" if row.yaxis_even_residual <= row.yaxis_odd_residual else "odd"
        )
        print(
            f"  epsilon={row.epsilon:.3f}: a={row.a:.6e}, "
            f"BIC={'PASS' if row.unique_bic_verified else 'FAIL'}, "
            f"asymptotic={'PASS' if row.asymptotic_agreement_verified else 'FAIL'}, "
            f"sv_min={row.sv_min_final:.3e}, "
            f"mesh_change={row.final_relative_mesh_change:.3%}, "
            f"yparity={parity_label}({parity_res:.3e}), "
            f"Rprop={max_prop:.3e}, "
            f"sigma_error={row.relative_error_sigma:.3%}"
        )

    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
