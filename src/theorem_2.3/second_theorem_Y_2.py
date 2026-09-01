"""
Theorem 2.3(iv): y-symmetric BIC with numerical tuning of the vertical position.

This is a second version of the y-symmetric validation.

The theorem states, for sufficiently small epsilon, that a y-axis-symmetric
obstacle (X odd, Y even, hence nu=0) supports a unique embedded eigenvalue in

    Lambda_1 <= k^2 < Lambda_2

when its vertical position satisfies

    a = epsilon*a_1 + O(epsilon^2),

and

    k^2 = Lambda_2 - sigma^2,
    sigma = epsilon^2*pi^3*mu/b^3 + O(epsilon^3 log epsilon).

For the explicit small-beta example used here,

    a_1 = -beta/12 + O(beta^2).

VERSION 2 DOES NOT SET a = epsilon*a_1 AND STOP THERE.

For every epsilon it uses

    a_asym = epsilon*(-beta/12)

as the centre of an O(epsilon^2) search window

    a in [a_asym - C_a epsilon^2, a_asym + C_a epsilon^2].

For each trial a:
  1. the BEM matrix is assembled,
  2. the spectral minimum in k near Lambda_2 is located,
  3. the corresponding boundary state is reconstructed,
  4. its projection onto the first open waveguide channel is measured.

Then we solve numerically

    a_BIC = argmin_a R_prop(a),

where R_prop is the maximum left/right fraction of the field in the first
propagating transverse mode.  Finally we report

    Delta a = a_BIC - epsilon*a_1,
    Delta a / epsilon^2,

which is the natural quantity for checking the O(epsilon^2) placement
correction predicted by the theorem.

After a_BIC is found we:
  - validate the spectral branch under mesh refinement,
  - check that the open-channel projection is small,
  - screen the complete interval Lambda_1 < k^2 < Lambda_2 for additional
    resolved non-radiating BIC candidates,
  - compare sigma_BEM with the leading asymptotic sigma.

Important terminology:
  sigma             = sqrt(Lambda_2-k^2), spectral quantity in the theorem.
  sv_min(A)          = smallest singular value of the BEM matrix.
  R_prop             = fraction in the first open propagating channel.
  a_asym             = epsilon*(-beta/12), leading placement prediction.
  a_BIC              = numerically tuned placement minimizing R_prop.

The numerical tolerances below are reporting/diagnostic criteria; they are
not constants appearing in Theorem 2.3.

IMPORTANT ABOUT THE O(epsilon^2) TEST:
The cleanest test uses the full theorem value of a_1.  If we use only

    a_1 = -beta/12 + O(beta^2),

then Delta a also contains the error coming from O(beta^2) in a_1.  Therefore
Delta a/epsilon^2 should be interpreted cautiously until the full a_1 integral
is supplied through Config.a1_override.
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

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

import dipole_theorem as dip
import lattice_sums as lattice

PI = np.pi


@dataclass(frozen=True)
class Config:
    # ------------------------------------------------------------------
    # Physical / geometric parameters
    # ------------------------------------------------------------------
    b: float = 1.0
    shape_beta: float = 0.10

    # If you later evaluate the full a_1 integral in Theorem 2.3(iv), put the
    # numerical value here.  Leaving it as None uses the small-beta formula
    # a_1 = -beta/12 + O(beta^2).
    a1_override: float | None = None

    # Keep this small at first because the placement optimization is nested:
    # each a evaluation contains a spectral minimization in k.
    epsilon_values: tuple[float, ...] = (0.05, 0.07, 0.09, 0.11)

    # ------------------------------------------------------------------
    # BEM / Green function
    # ------------------------------------------------------------------
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # ------------------------------------------------------------------
    # Placement tuning a = epsilon*a1 + O(epsilon^2)
    # ------------------------------------------------------------------
    # Search half-width = a_window_C * epsilon^2.
    # If the optimized a lands near an edge, increase this value.
    a_window_C: float = 1.0

    # Coarse samples before bounded refinement in a.
    a_coarse_points: int = 7

    # We optimize a at these BEM resolutions. The last one defines a_BIC.
    placement_optimization_M: tuple[int, ...] = (16, 24, 32)

    # After a_BIC is fixed, independently refine the spectral branch with M.
    final_refinement_M: tuple[int, ...] = (16, 24, 32, 40, 48)

    a_minimizer_xatol: float = 1.0e-8

    # For M>first placement mesh, search locally around the previous a_BIC.
    # local half-width = factor * epsilon^2, clipped to the full O(eps^2) box.
    a_follow_halfwidth_factor: float = 0.35

    # ------------------------------------------------------------------
    # Spectral search near Lambda_2
    # delta_2 = Lambda_2-k^2 = sigma^2
    # ------------------------------------------------------------------
    expected_delta_lower_factor: float = 0.15
    expected_delta_upper_factor: float = 6.0
    upper_cutoff_delta_floor: float = 1.0e-9
    k_minimizer_xatol: float = 2.0e-10

    # ------------------------------------------------------------------
    # Spectral/BIC diagnostics
    # ------------------------------------------------------------------
    minimum_drop_factor: float = 100.0
    near_singular_tolerance: float = 1.0e-4
    mesh_sigma_relative_tolerance: float = 0.03

    # The point of this version is to reduce radiation below the leading-a
    # calculation. 1e-4 is deliberately stricter than the first script's 1e-2.
    propagating_ratio_tolerance: float = 1.0e-4
    projection_x_over_b: float = 2.0
    projection_quadrature_points: int = 64

    geometry_symmetry_tolerance: float = 1.0e-12
    nu_tolerance: float = 1.0e-8

    # Leading asymptotic sigma accuracy is reported separately.
    relative_sigma_error_tolerance: float = 0.05

    # ------------------------------------------------------------------
    # Whole-band uniqueness screen at the tuned a_BIC
    # ------------------------------------------------------------------
    uniqueness_scan_M: int = 16
    uniqueness_refine_M: int = 32
    uniqueness_linear_points: int = 32
    uniqueness_log_points: int = 32
    lower_cutoff_kb_margin: float = 2.0e-4
    upper_cutoff_kb_margin: float = 1.0e-7
    additional_candidate_kb_tolerance: float = 2.0e-3
    expected_candidate_match_kb_tolerance: float = 5.0e-4

    output_directory: str = "theorem_2_3_y_symmetry_tuned_a_validation"


@dataclass
class BranchState:
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
    yaxis_even_residual: float
    yaxis_odd_residual: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    first_mode_coefficient_left: float
    first_mode_coefficient_right: float
    field_norm_left: float
    field_norm_right: float


@dataclass
class PlacementOptimizationRow:
    epsilon: float
    M: int
    a_asymptotic: float
    search_left: float
    search_right: float
    a_bic: float
    delta_a: float
    delta_a_over_epsilon2: float
    kb_bic: float
    sigma_bem: float
    sv_min: float
    drop_factor: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    max_propagating_ratio: float
    optimum_is_interior: bool


@dataclass
class MeshRefinementRow:
    epsilon: float
    a_bic: float
    M: int
    kb: float
    sigma_bem: float
    sv_min: float
    drop_factor: float
    relative_sigma_change_from_previous: float
    yaxis_even_residual: float
    yaxis_odd_residual: float
    propagating_ratio_left: float
    propagating_ratio_right: float


@dataclass
class AdditionalCandidate:
    epsilon: float
    a_bic: float
    bracket_left: float
    bracket_right: float
    kb_coarse: float
    sv_min_coarse: float
    kb_fine: float
    sv_min_fine: float
    fine_drop_factor: float
    kb_shift: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    persistent_nonradiating: bool
    matches_expected_bic: bool


@dataclass
class ValidationResult:
    epsilon: float
    beta: float
    mu: float
    nu: float
    a1_leading: float
    a_asymptotic: float
    a_bic: float
    delta_a: float
    delta_a_over_epsilon2: float
    a_search_halfwidth: float
    a_optimum_is_interior: bool
    leading_a_propagating_ratio: float
    tuned_a_propagating_ratio: float
    radiation_reduction_factor: float
    lambda_1: float
    lambda_2: float
    kb_asymptotic: float
    kb_numerical: float
    sigma_asymptotic: float
    sigma_numerical: float
    sv_min_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    yaxis_even_residual: float
    yaxis_odd_residual: float
    relative_error_kb: float
    relative_error_sigma: float
    symmetry_conditions_verified: bool
    tuned_placement_consistent_with_O_epsilon2: bool
    expected_bic_resolved: bool
    nonradiating_bic_verified: bool
    additional_resolved_bics: int
    uniqueness_screen_passed: bool
    unique_bic_verified: bool
    asymptotic_agreement_verified: bool


CONFIG = Config()


# ===========================================================================
# Spectral / asymptotic quantities
# ===========================================================================


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
    mu, nu = dip.dipole(shape_beta, 0.0, 0.0, "y")
    return float(mu), float(nu)


def a1_leading(config: Config) -> float:
    if config.a1_override is not None:
        return float(config.a1_override)
    return -config.shape_beta / 12.0


def a_asymptotic(epsilon: float, config: Config) -> float:
    return epsilon * a1_leading(config)


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    mu, _ = dipole_moments(config.shape_beta)
    sigma = epsilon**2 * PI**3 / config.b**3 * mu

    k_squared = lambda_2(config) - sigma**2
    if not (lambda_1(config) <= k_squared < lambda_2(config)):
        raise ValueError(
            "Leading asymptotic prediction is outside [Lambda_1,Lambda_2)."
        )
    return config.b * math.sqrt(k_squared), float(sigma)


def sigma_from_kb(kb: float, config: Config) -> float:
    k = kb / config.b
    return math.sqrt(max(lambda_2(config) - k**2, 0.0))


def kb_from_delta_2(delta_2: float, config: Config) -> float:
    lam1 = lambda_1(config)
    lam2 = lambda_2(config)
    max_gap = lam2 - lam1
    delta_2 = min(
        max(delta_2, config.upper_cutoff_delta_floor),
        max_gap * (1.0 - 1.0e-10),
    )
    return config.b * math.sqrt(lam2 - delta_2)


# ===========================================================================
# Geometry
# ===========================================================================


def obstacle_geometry(
    t: np.ndarray | float,
    epsilon: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Obstacle coordinates relative to its centre."""
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
    t = np.linspace(-PI, PI, 257)
    Xp, Yp, *_ = obstacle_geometry(t, 1.0, config)
    Xm, Ym, *_ = obstacle_geometry(-t, 1.0, config)
    return (
        float(np.max(np.abs(Xm + Xp))),
        float(np.max(np.abs(Ym - Yp))),
    )


def symmetry_conditions_verified(
    config: Config,
) -> tuple[bool, float, float, float, float]:
    x_odd, y_even = geometry_symmetry_residual(config)
    mu, nu = dipole_moments(config.shape_beta)
    passed = bool(
        x_odd <= config.geometry_symmetry_tolerance
        and y_even <= config.geometry_symmetry_tolerance
        and abs(nu) <= config.nu_tolerance
    )
    return passed, x_odd, y_even, mu, nu


def obstacle_fits_guide(epsilon: float, a: float, config: Config) -> bool:
    t = np.linspace(0.0, 2.0 * PI, 1025)
    _, Y, *_ = obstacle_geometry(t, epsilon, config)
    return bool(np.max(np.abs(Y + a)) < config.b)


def placement_bounds(epsilon: float, config: Config) -> tuple[float, float]:
    centre = a_asymptotic(epsilon, config)
    halfwidth = config.a_window_C * epsilon**2
    return centre - halfwidth, centre + halfwidth


# ===========================================================================
# BEM assembly for arbitrary a
# ===========================================================================


def boundary_nodes(M: int) -> np.ndarray:
    return (np.arange(M) + 0.5) * (2.0 * PI / M)


def make_green_functions(
    kb: float,
    a: float,
    config: Config,
) -> tuple[Callable[..., complex], Callable[..., complex]]:
    b = config.b
    d = 2.0 * b
    k = kb / b

    coefficients = lattice.lattice_sums(
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
            coefficients,
            k,
            d,
        )

    def green_regularized(x: float, y: float, xi: float, eta: float) -> complex:
        return lattice.greens_dirichlet_reg(
            x,
            y + b + a,
            xi,
            eta + b + a,
            coefficients,
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
    a: float,
    M: int,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, Callable[..., complex]]:
    G, G_regularized = make_green_functions(kb, a, config)
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

    A = np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted
    return A, theta, G


def smallest_singular_pair(A: np.ndarray) -> tuple[float, np.ndarray]:
    _, singular_values, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= max(float(np.linalg.norm(vector)), 1.0e-30)

    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))
    return float(singular_values[-1]), vector


def smallest_singular_value(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> float:
    A, _, _ = assemble_matrix(kb, epsilon, a, M, config)
    return float(np.linalg.svd(A, compute_uv=False)[-1])


# ===========================================================================
# Field reconstruction / BIC diagnostics
# ===========================================================================


def boundary_parity_residuals(vector: np.ndarray) -> tuple[float, float]:
    reflected = vector[::-1]
    norm = max(float(np.linalg.norm(vector)), 1.0e-30)
    even = float(np.linalg.norm(vector - reflected) / norm)
    odd = float(np.linalg.norm(vector + reflected) / norm)
    return even, odd


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
        G, x, y, float(xi), float(eta), config.finite_difference_step
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
        kernel = np.array(
            [
                weighted_kernel_at_field_point(
                    x, float(y), float(t), epsilon, config, G
                )
                for t in theta
            ],
            dtype=np.complex128,
        )
        field[i] = delta_theta * np.dot(kernel, boundary_vector)
    return field


def first_propagating_component(
    epsilon: float,
    a: float,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    G: Callable[..., complex],
    config: Config,
) -> dict[str, float]:
    """Projection onto the only open transverse mode in (Lambda_1,Lambda_2)."""
    b = config.b
    x_far = config.projection_x_over_b * b

    nodes, weights = leggauss(config.projection_quadrature_points)
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
            x, y_values, boundary_vector, theta, epsilon, config, G
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


# ===========================================================================
# Spectral branch at a prescribed vertical position a
# ===========================================================================


def expected_k_bracket(epsilon: float, config: Config) -> tuple[float, float]:
    _, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2
    max_gap = lambda_2(config) - lambda_1(config)

    delta_far = min(
        delta_asym * config.expected_delta_upper_factor,
        max_gap * 0.95,
    )
    delta_near = max(
        delta_asym * config.expected_delta_lower_factor,
        config.upper_cutoff_delta_floor,
    )

    left = kb_from_delta_2(delta_far, config)
    right = kb_from_delta_2(delta_near, config)
    if not left < right:
        raise RuntimeError("Invalid expected-BIC k bracket.")
    return left, right


def solve_branch_at_a(
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> BranchState:
    """
    For a fixed a, locate the near-singular spectral state in k and measure
    its propagating-channel content.
    """
    if not obstacle_fits_guide(epsilon, a, config):
        raise ValueError("Trial obstacle position lies outside the guide.")

    left, right = expected_k_bracket(epsilon, config)
    tiny = np.finfo(float).tiny

    left_value = smallest_singular_value(left, epsilon, a, M, config)
    right_value = smallest_singular_value(right, epsilon, a, M, config)

    result = minimize_scalar(
        lambda kb: math.log10(
            max(
                smallest_singular_value(float(kb), epsilon, a, M, config),
                tiny,
            )
        ),
        bounds=(left, right),
        method="bounded",
        options={"xatol": config.k_minimizer_xatol},
    )
    if not result.success:
        raise RuntimeError(f"k minimization failed at a={a:.8e}, M={M}")

    kb = float(result.x)
    A, theta, G = assemble_matrix(kb, epsilon, a, M, config)
    sv_min, vector = smallest_singular_pair(A)

    drop = min(left_value, right_value) / max(sv_min, tiny)
    width = right - left
    edge = 0.01 * width
    interior = (kb > left + edge) and (kb < right - edge)

    even_res, odd_res = boundary_parity_residuals(vector)
    prop = first_propagating_component(epsilon, a, vector, theta, G, config)

    return BranchState(
        epsilon=epsilon,
        a=a,
        M=M,
        kb=kb,
        sigma_bem=sigma_from_kb(kb, config),
        sv_min=sv_min,
        left_value=left_value,
        right_value=right_value,
        drop_factor=drop,
        minimum_is_interior=interior,
        yaxis_even_residual=even_res,
        yaxis_odd_residual=odd_res,
        propagating_ratio_left=prop["ratio_left"],
        propagating_ratio_right=prop["ratio_right"],
        first_mode_coefficient_left=prop["coefficient_left"],
        first_mode_coefficient_right=prop["coefficient_right"],
        field_norm_left=prop["field_norm_left"],
        field_norm_right=prop["field_norm_right"],
    )


def max_rprop(state: BranchState) -> float:
    return max(state.propagating_ratio_left, state.propagating_ratio_right)


def branch_is_spectrally_resolved(state: BranchState, config: Config) -> bool:
    return bool(
        state.minimum_is_interior
        and state.sv_min <= config.near_singular_tolerance
        and state.drop_factor >= config.minimum_drop_factor
    )


# ===========================================================================
# Optimize a_BIC = argmin_a R_prop(a)
# ===========================================================================


def placement_objective(state: BranchState, config: Config) -> float:
    """
    Main quantity is R_prop. A large penalty prevents the outer optimizer from
    preferring a point where the tracked spectral branch is not actually
    resolved.
    """
    rprop = max_rprop(state)
    if branch_is_spectrally_resolved(state, config):
        return rprop

    penalty = 1.0
    if state.sv_min > config.near_singular_tolerance:
        penalty += min(state.sv_min / config.near_singular_tolerance, 100.0)
    if state.drop_factor < config.minimum_drop_factor:
        penalty += min(
            config.minimum_drop_factor / max(state.drop_factor, 1.0e-30),
            100.0,
        )
    if not state.minimum_is_interior:
        penalty += 100.0
    return penalty + rprop


def optimize_a_for_M(
    epsilon: float,
    M: int,
    config: Config,
    previous_a: float | None = None,
) -> tuple[PlacementOptimizationRow, BranchState, np.ndarray, np.ndarray]:
    """
    Coarse-to-fine minimization of R_prop(a).

    Returns the optimization row, final BranchState, coarse a samples, and
    coarse R_prop values. The latter are useful for diagnostic plots.
    """
    full_left, full_right = placement_bounds(epsilon, config)
    a0 = a_asymptotic(epsilon, config)

    if previous_a is None:
        left, right = full_left, full_right
    else:
        half = config.a_follow_halfwidth_factor * epsilon**2
        left = max(full_left, previous_a - half)
        right = min(full_right, previous_a + half)

    if not left < right:
        raise RuntimeError("Invalid a optimization interval.")

    cache: dict[float, BranchState] = {}

    def evaluate(a_value: float) -> BranchState:
        key = round(float(a_value), 14)
        if key not in cache:
            cache[key] = solve_branch_at_a(epsilon, float(a_value), M, config)
        return cache[key]

    # Coarse profile protects the bounded optimizer from a poor global box.
    coarse_a = np.linspace(left, right, config.a_coarse_points)
    coarse_states = [evaluate(float(a)) for a in coarse_a]
    coarse_obj = np.array(
        [placement_objective(state, config) for state in coarse_states], dtype=float
    )
    best_i = int(np.argmin(coarse_obj))

    local_left = coarse_a[max(best_i - 1, 0)]
    local_right = coarse_a[min(best_i + 1, len(coarse_a) - 1)]
    if local_left == local_right:
        local_left, local_right = left, right

    result = minimize_scalar(
        lambda a_value: placement_objective(evaluate(float(a_value)), config),
        bounds=(float(local_left), float(local_right)),
        method="bounded",
        options={"xatol": config.a_minimizer_xatol},
    )
    if not result.success:
        raise RuntimeError(f"a minimization failed for epsilon={epsilon}, M={M}")

    a_bic = float(result.x)
    state = evaluate(a_bic)

    full_width = full_right - full_left
    edge_margin = 0.01 * full_width
    optimum_interior = (
        a_bic > full_left + edge_margin and a_bic < full_right - edge_margin
    )

    delta_a = a_bic - a0
    row = PlacementOptimizationRow(
        epsilon=epsilon,
        M=M,
        a_asymptotic=a0,
        search_left=left,
        search_right=right,
        a_bic=a_bic,
        delta_a=delta_a,
        delta_a_over_epsilon2=delta_a / epsilon**2,
        kb_bic=state.kb,
        sigma_bem=state.sigma_bem,
        sv_min=state.sv_min,
        drop_factor=state.drop_factor,
        propagating_ratio_left=state.propagating_ratio_left,
        propagating_ratio_right=state.propagating_ratio_right,
        max_propagating_ratio=max_rprop(state),
        optimum_is_interior=optimum_interior,
    )

    coarse_rprop = np.array([max_rprop(s) for s in coarse_states], dtype=float)
    return row, state, coarse_a, coarse_rprop


def run_placement_optimization(
    epsilon: float,
    config: Config,
) -> tuple[
    list[PlacementOptimizationRow], float, list[tuple[int, np.ndarray, np.ndarray]]
]:
    history: list[PlacementOptimizationRow] = []
    profiles: list[tuple[int, np.ndarray, np.ndarray]] = []
    previous_a: float | None = None

    print("  tuning a_BIC by minimizing the open-channel projection:")
    for M in config.placement_optimization_M:
        row, _, coarse_a, coarse_rprop = optimize_a_for_M(
            epsilon, M, config, previous_a
        )
        history.append(row)
        profiles.append((M, coarse_a, coarse_rprop))
        previous_a = row.a_bic

        print(
            f"    M={M:>2}: a_BIC={row.a_bic:+.10e}, "
            f"Delta a={row.delta_a:+.3e}, "
            f"Delta a/eps^2={row.delta_a_over_epsilon2:+.6f}, "
            f"kb={row.kb_bic:.12f}, "
            f"Rprop={row.max_propagating_ratio:.3e}, "
            f"sv_min={row.sv_min:.3e}, drop={row.drop_factor:.2e}, "
            f"a-interior={'yes' if row.optimum_is_interior else 'NO'}"
        )

    return history, history[-1].a_bic, profiles


# ===========================================================================
# Fixed-a mesh refinement after tuning
# ===========================================================================


def run_fixed_a_mesh_refinement(
    epsilon: float,
    a_bic: float,
    config: Config,
) -> list[MeshRefinementRow]:
    rows: list[MeshRefinementRow] = []
    previous_sigma = math.nan

    print("  fixed-a_BIC spectral mesh refinement:")
    for M in config.final_refinement_M:
        state = solve_branch_at_a(epsilon, a_bic, M, config)
        change = (
            abs(state.sigma_bem - previous_sigma) / max(abs(state.sigma_bem), 1.0e-30)
            if np.isfinite(previous_sigma)
            else math.nan
        )
        row = MeshRefinementRow(
            epsilon=epsilon,
            a_bic=a_bic,
            M=M,
            kb=state.kb,
            sigma_bem=state.sigma_bem,
            sv_min=state.sv_min,
            drop_factor=state.drop_factor,
            relative_sigma_change_from_previous=change,
            yaxis_even_residual=state.yaxis_even_residual,
            yaxis_odd_residual=state.yaxis_odd_residual,
            propagating_ratio_left=state.propagating_ratio_left,
            propagating_ratio_right=state.propagating_ratio_right,
        )
        rows.append(row)
        previous_sigma = state.sigma_bem

        change_text = f"{change:.3%}" if np.isfinite(change) else "--"
        print(
            f"    M={M:>2}: kb={state.kb:.12f}, "
            f"sigma_BEM={state.sigma_bem:.8e}, "
            f"sv_min={state.sv_min:.3e}, drop={state.drop_factor:.2e}, "
            f"mesh_change={change_text}, Rprop={max_rprop(state):.3e}"
        )

    return rows


def fixed_a_bic_is_resolved(
    rows: list[MeshRefinementRow],
    config: Config,
) -> tuple[bool, bool]:
    if len(rows) < 2:
        return False, False

    final = rows[-1]
    previous = rows[-2]
    change = abs(final.sigma_bem - previous.sigma_bem) / max(
        abs(final.sigma_bem), 1.0e-30
    )

    spectral = bool(
        final.sv_min <= config.near_singular_tolerance
        and final.drop_factor >= config.minimum_drop_factor
        and change <= config.mesh_sigma_relative_tolerance
    )
    radiation = bool(
        final.propagating_ratio_left <= config.propagating_ratio_tolerance
        and final.propagating_ratio_right <= config.propagating_ratio_tolerance
    )
    return spectral, bool(spectral and radiation)


# ===========================================================================
# Whole-band uniqueness screen at a_BIC
# ===========================================================================


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

    deltas = np.geomspace(
        config.upper_cutoff_delta_floor,
        (lam2 - lam1) * (1.0 - 1.0e-6),
        config.uniqueness_log_points,
    )
    logarithmic = config.b * np.sqrt(lam2 - deltas)

    factors = np.array([0.15, 0.25, 0.5, 0.75, 1.0, 1.5, 2.5, 5.0])
    local_delta = delta_asym * factors
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


def expected_window_overlap(
    left: float,
    right: float,
    epsilon: float,
    config: Config,
) -> bool:
    eleft, eright = expected_k_bracket(epsilon, config)
    return max(left, eleft) <= min(right, eright)


def refine_k_bracket(
    left: float,
    right: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[float, BranchState]:
    tiny = np.finfo(float).tiny
    left_value = smallest_singular_value(left, epsilon, a, M, config)
    right_value = smallest_singular_value(right, epsilon, a, M, config)

    result = minimize_scalar(
        lambda kb: math.log10(
            max(smallest_singular_value(float(kb), epsilon, a, M, config), tiny)
        ),
        bounds=(left, right),
        method="bounded",
        options={"xatol": config.k_minimizer_xatol},
    )
    if not result.success:
        raise RuntimeError("Additional candidate k refinement failed.")

    kb = float(result.x)
    A, theta, G = assemble_matrix(kb, epsilon, a, M, config)
    sv_min, vector = smallest_singular_pair(A)
    drop = min(left_value, right_value) / max(sv_min, tiny)

    even_res, odd_res = boundary_parity_residuals(vector)
    prop = first_propagating_component(epsilon, a, vector, theta, G, config)

    state = BranchState(
        epsilon=epsilon,
        a=a,
        M=M,
        kb=kb,
        sigma_bem=sigma_from_kb(kb, config),
        sv_min=sv_min,
        left_value=left_value,
        right_value=right_value,
        drop_factor=drop,
        minimum_is_interior=True,
        yaxis_even_residual=even_res,
        yaxis_odd_residual=odd_res,
        propagating_ratio_left=prop["ratio_left"],
        propagating_ratio_right=prop["ratio_right"],
        first_mode_coefficient_left=prop["coefficient_left"],
        first_mode_coefficient_right=prop["coefficient_right"],
        field_norm_left=prop["field_norm_left"],
        field_norm_right=prop["field_norm_right"],
    )
    return kb, state


def whole_band_uniqueness_screen(
    epsilon: float,
    a_bic: float,
    expected_kb: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, list[AdditionalCandidate]]:
    """
    Screen the complete band at the tuned placement.

    Importantly, every sampled local minimum is refined.  We do NOT discard
    the whole asymptotic search window, because doing so could hide a second
    eigenvalue close to the expected one.  After refinement we simply identify
    the minimum that matches the already-resolved BIC by proximity in kb; all
    other persistent non-radiating minima count as additional BICs.
    """
    scan = build_uniqueness_grid(epsilon, config)
    values = np.empty(len(scan), dtype=float)

    print(
        f"  whole-band uniqueness screen at tuned a_BIC: {len(scan)} points "
        f"(M={config.uniqueness_scan_M})"
    )
    for i, kb in enumerate(scan):
        values[i] = smallest_singular_value(
            float(kb), epsilon, a_bic, config.uniqueness_scan_M, config
        )
        if (i + 1) % 16 == 0 or i + 1 == len(scan):
            print(f"    {i + 1:>3}/{len(scan)}")

    brackets = sampled_local_minimum_brackets(scan, values)
    print(f"  sampled local minima: {len(brackets)} total; refining all of them")

    candidates: list[AdditionalCandidate] = []
    for j, (left, right) in enumerate(brackets, start=1):
        print(f"    checking minimum {j}/{len(brackets)}")
        kb_c, coarse = refine_k_bracket(
            left, right, epsilon, a_bic, config.uniqueness_scan_M, config
        )
        kb_f, fine = refine_k_bracket(
            left, right, epsilon, a_bic, config.uniqueness_refine_M, config
        )
        shift = abs(kb_f - kb_c)
        matches_expected = (
            abs(kb_f - expected_kb) <= config.expected_candidate_match_kb_tolerance
        )

        persistent = bool(
            fine.sv_min <= config.near_singular_tolerance
            and fine.drop_factor >= config.minimum_drop_factor
            and shift <= config.additional_candidate_kb_tolerance
            and max_rprop(fine) <= config.propagating_ratio_tolerance
        )

        candidates.append(
            AdditionalCandidate(
                epsilon=epsilon,
                a_bic=a_bic,
                bracket_left=left,
                bracket_right=right,
                kb_coarse=kb_c,
                sv_min_coarse=coarse.sv_min,
                kb_fine=kb_f,
                sv_min_fine=fine.sv_min,
                fine_drop_factor=fine.drop_factor,
                kb_shift=shift,
                propagating_ratio_left=fine.propagating_ratio_left,
                propagating_ratio_right=fine.propagating_ratio_right,
                persistent_nonradiating=persistent,
                matches_expected_bic=matches_expected,
            )
        )

        label = "expected BIC" if matches_expected else "other"
        print(
            f"      kb(M={config.uniqueness_refine_M})={kb_f:.12f}, "
            f"sv_min={fine.sv_min:.3e}, drop={fine.drop_factor:.2e}, "
            f"shift={shift:.3e}, Rprop={max_rprop(fine):.3e}, "
            f"persistent={'YES' if persistent else 'no'}, class={label}"
        )

    return scan, values, candidates


# ===========================================================================
# Per-epsilon validation
# ===========================================================================


def validate_epsilon(
    epsilon: float,
    config: Config,
) -> tuple[
    ValidationResult,
    list[PlacementOptimizationRow],
    list[MeshRefinementRow],
    np.ndarray,
    np.ndarray,
    list[AdditionalCandidate],
    list[tuple[int, np.ndarray, np.ndarray]],
]:
    symmetry_ok, _, _, mu, nu = symmetry_conditions_verified(config)
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
    a0 = a_asymptotic(epsilon, config)

    # Baseline radiation using only the leading placement at the finest
    # placement-optimization M.
    baseline_M = config.final_refinement_M[-1]
    leading_state = solve_branch_at_a(epsilon, a0, baseline_M, config)
    leading_rprop = max_rprop(leading_state)

    placement_history, a_bic, profiles = run_placement_optimization(epsilon, config)
    final_placement = placement_history[-1]

    mesh_rows = run_fixed_a_mesh_refinement(epsilon, a_bic, config)
    spectral_ok, nonradiating_ok = fixed_a_bic_is_resolved(mesh_rows, config)
    final_mesh = mesh_rows[-1]

    tuned_rprop = max(
        final_mesh.propagating_ratio_left, final_mesh.propagating_ratio_right
    )
    reduction = leading_rprop / max(tuned_rprop, 1.0e-30)

    scan, scan_values, additional = whole_band_uniqueness_screen(
        epsilon, a_bic, final_mesh.kb, config
    )
    additional_resolved = [
        c
        for c in additional
        if c.persistent_nonradiating and not c.matches_expected_bic
    ]
    uniqueness_ok = len(additional_resolved) == 0

    delta_a = a_bic - a0
    halfwidth = config.a_window_C * epsilon**2
    placement_scaling_ok = bool(
        final_placement.optimum_is_interior
        and abs(delta_a) <= halfwidth
        and obstacle_fits_guide(epsilon, a_bic, config)
    )

    error_kb = abs(final_mesh.kb - kb_asym) / abs(kb_asym)
    error_sigma = abs(final_mesh.sigma_bem - sigma_asym) / abs(sigma_asym)
    asymptotic_ok = error_sigma <= config.relative_sigma_error_tolerance

    unique_bic = bool(
        symmetry_ok and placement_scaling_ok and nonradiating_ok and uniqueness_ok
    )

    result = ValidationResult(
        epsilon=epsilon,
        beta=config.shape_beta,
        mu=mu,
        nu=nu,
        a1_leading=a1_leading(config),
        a_asymptotic=a0,
        a_bic=a_bic,
        delta_a=delta_a,
        delta_a_over_epsilon2=delta_a / epsilon**2,
        a_search_halfwidth=halfwidth,
        a_optimum_is_interior=final_placement.optimum_is_interior,
        leading_a_propagating_ratio=leading_rprop,
        tuned_a_propagating_ratio=tuned_rprop,
        radiation_reduction_factor=reduction,
        lambda_1=lambda_1(config),
        lambda_2=lambda_2(config),
        kb_asymptotic=kb_asym,
        kb_numerical=final_mesh.kb,
        sigma_asymptotic=sigma_asym,
        sigma_numerical=final_mesh.sigma_bem,
        sv_min_final=final_mesh.sv_min,
        final_drop_factor=final_mesh.drop_factor,
        final_relative_mesh_change=final_mesh.relative_sigma_change_from_previous,
        yaxis_even_residual=final_mesh.yaxis_even_residual,
        yaxis_odd_residual=final_mesh.yaxis_odd_residual,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        symmetry_conditions_verified=symmetry_ok,
        tuned_placement_consistent_with_O_epsilon2=placement_scaling_ok,
        expected_bic_resolved=spectral_ok,
        nonradiating_bic_verified=nonradiating_ok,
        additional_resolved_bics=len(additional_resolved),
        uniqueness_screen_passed=uniqueness_ok,
        unique_bic_verified=unique_bic,
        asymptotic_agreement_verified=asymptotic_ok,
    )

    return (
        result,
        placement_history,
        mesh_rows,
        scan,
        scan_values,
        additional,
        profiles,
    )


# ===========================================================================
# Output helpers
# ===========================================================================


def write_dataclass_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def plot_a_profile(
    epsilon: float,
    profiles: list[tuple[int, np.ndarray, np.ndarray]],
    a0: float,
    a_bic: float,
    output_directory: Path,
) -> None:
    for M, a_values, rprop in profiles:
        plt.figure(figsize=(7, 4.5))
        plt.semilogy(a_values, rprop, "o-")
        plt.axvline(a0, linestyle="--", label=r"$a_{\mathrm{asym}}$")
        plt.axvline(a_bic, linestyle=":", label=r"$a_{\mathrm{BIC}}$")
        plt.xlabel(r"$a$")
        plt.ylabel(r"$R_{\mathrm{prop}}$")
        plt.title(rf"Placement tuning, $\varepsilon={epsilon:.2f}$, $M={M}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            output_directory / f"a_profile_epsilon_{epsilon:.2f}_M{M}.png",
            dpi=180,
        )
        plt.close()


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
    plt.title(rf"Theorem 2.3(iv), tuned $a$, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"full_band_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_summary(
    results: list[ValidationResult],
    output_directory: Path,
    config: Config,
) -> None:
    eps = np.array([r.epsilon for r in results])
    a0 = np.array([r.a_asymptotic for r in results])
    abic = np.array([r.a_bic for r in results])
    scaled = np.array([r.delta_a_over_epsilon2 for r in results])
    r0 = np.array([r.leading_a_propagating_ratio for r in results])
    rbic = np.array([r.tuned_a_propagating_ratio for r in results])
    kb_asym = np.array([r.kb_asymptotic for r in results])
    kb_num = np.array([r.kb_numerical for r in results])
    sigma_error = np.array([r.relative_error_sigma for r in results])

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
    plt.plot(eps, a0, "s--", label=r"$a_{\mathrm{asym}}=\varepsilon a_1$")
    plt.plot(eps, abic, "o-", label=r"$a_{\mathrm{BIC}}$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$a$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "a_bic_vs_leading.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, scaled, "o-")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$(a_{\mathrm{BIC}}-\varepsilon a_1)/\varepsilon^2$")
    plt.tight_layout()
    plt.savefig(output_directory / "scaled_a_correction.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(eps, r0, "s--", label=r"leading $a=\varepsilon a_1$")
    plt.semilogy(eps, rbic, "o-", label=r"tuned $a_{\mathrm{BIC}}$")
    plt.axhline(
        config.propagating_ratio_tolerance,
        linestyle=":",
        label="BIC radiation tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$R_{\mathrm{prop}}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "radiation_before_after_tuning.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, sigma_error, "o-")
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


def print_result(result: ValidationResult, config: Config) -> None:
    print("\n  --- tuned-placement validation result ---")
    print(f"  a_asym = epsilon*a1      = {result.a_asymptotic:+.10e}")
    print(f"  a_BIC (numerical)        = {result.a_bic:+.10e}")
    print(f"  Delta a                  = {result.delta_a:+.6e}")
    print(f"  Delta a / epsilon^2      = {result.delta_a_over_epsilon2:+.8f}")
    print(
        f"  O(epsilon^2) search consistency: "
        f"{'PASS' if result.tuned_placement_consistent_with_O_epsilon2 else 'FAIL'}"
    )
    print(f"  Rprop at leading a       = {result.leading_a_propagating_ratio:.3e}")
    print(f"  Rprop at tuned a_BIC     = {result.tuned_a_propagating_ratio:.3e}")
    print(f"  radiation reduction      = {result.radiation_reduction_factor:.2e}x")
    print(
        f"  spectral BIC resolved    = "
        f"{'PASS' if result.expected_bic_resolved else 'FAIL'}"
    )
    print(
        f"  non-radiating diagnostic = "
        f"{'PASS' if result.nonradiating_bic_verified else 'FAIL'}"
    )
    print(f"  additional resolved BICs = {result.additional_resolved_bics}")
    print(
        f"  uniqueness screen        = "
        f"{'PASS' if result.uniqueness_screen_passed else 'FAIL'}"
    )
    print(
        f"  unique tuned BIC         = "
        f"{'PASS' if result.unique_bic_verified else 'FAIL'}"
    )
    print(f"  kb asymptotic            = {result.kb_asymptotic:.12f}")
    print(f"  kb BEM                   = {result.kb_numerical:.12f}")
    print(f"  sigma asymptotic         = {result.sigma_asymptotic:.8e}")
    print(f"  sigma BEM                = {result.sigma_numerical:.8e}")
    print(f"  final sv_min(A)          = {result.sv_min_final:.3e}")
    print(f"  final minimum drop       = {result.final_drop_factor:.2e}")
    print(f"  relative sigma error     = {result.relative_error_sigma:.3%}")
    print(
        f"  asymptotic accuracy <= {config.relative_sigma_error_tolerance:.1%}: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )


# ===========================================================================
# Main
# ===========================================================================


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    symmetry_ok, x_odd, y_even, mu, nu = symmetry_conditions_verified(config)

    print("=== Theorem 2.3(iv), version 2: tune a numerically ===")
    print(f"b = {config.b}")
    print(f"beta = {config.shape_beta}")
    if config.a1_override is None:
        print(
            f"a1 used = -beta/12 = {a1_leading(config):.12e} [small-beta approximation]"
        )
    else:
        print(f"a1 used = {a1_leading(config):.12e} [full/override value]")
    print(f"mu = {mu:.12e}")
    print(f"nu = {nu:.12e}")
    print(f"X-odd residual  = {x_odd:.3e}")
    print(f"Y-even residual = {y_even:.3e}")
    print(f"symmetry conditions: {'PASS' if symmetry_ok else 'FAIL'}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"Lambda_2 = {lambda_2(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")
    print(f"a search half-width = {config.a_window_C} * epsilon^2")
    print(f"placement optimization M = {config.placement_optimization_M}")
    print(f"final spectral refinement M = {config.final_refinement_M}")

    if not symmetry_ok:
        raise RuntimeError("Selected shape does not satisfy Theorem 2.3(iv) symmetry.")

    summary: list[ValidationResult] = []
    all_placement: list[PlacementOptimizationRow] = []
    all_mesh: list[MeshRefinementRow] = []
    all_additional: list[AdditionalCandidate] = []

    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.3f} ===")
        a0 = a_asymptotic(epsilon, config)
        left, right = placement_bounds(epsilon, config)
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)

        print(f"  leading a_asym = {a0:+.10e}")
        print(f"  O(eps^2) search window = [{left:+.10e}, {right:+.10e}]")
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted sigma = {sigma_asym:.8e}")

        (
            result,
            placement_history,
            mesh_rows,
            scan,
            scan_values,
            additional,
            profiles,
        ) = validate_epsilon(epsilon, config)

        summary.append(result)
        all_placement.extend(placement_history)
        all_mesh.extend(mesh_rows)
        all_additional.extend(additional)

        print_result(result, config)
        plot_a_profile(
            epsilon,
            profiles,
            result.a_asymptotic,
            result.a_bic,
            output_directory,
        )
        plot_full_band_scan(epsilon, scan, scan_values, result, output_directory)

    write_dataclass_csv(output_directory / "summary.csv", summary)
    write_dataclass_csv(output_directory / "placement_optimization.csv", all_placement)
    write_dataclass_csv(output_directory / "mesh_refinement.csv", all_mesh)
    write_dataclass_csv(output_directory / "additional_candidates.csv", all_additional)
    plot_summary(summary, output_directory, config)

    print("\n=== FINAL SUMMARY ===")
    for row in summary:
        print(
            f"  epsilon={row.epsilon:.3f}: "
            f"a_asym={row.a_asymptotic:+.6e}, "
            f"a_BIC={row.a_bic:+.6e}, "
            f"Delta a/eps^2={row.delta_a_over_epsilon2:+.6f}, "
            f"Rprop {row.leading_a_propagating_ratio:.3e} -> "
            f"{row.tuned_a_propagating_ratio:.3e} "
            f"({row.radiation_reduction_factor:.2e}x), "
            f"BIC={'PASS' if row.unique_bic_verified else 'FAIL'}, "
            f"sigma_error={row.relative_error_sigma:.3%}"
        )

    print("\nKey quantity for Theorem 2.3(iv):")
    print("  (a_BIC - epsilon*a1) / epsilon^2")
    print("If the O(epsilon^2) placement correction is resolved, this quantity should")
    print("remain O(1) and ideally approach a finite limit as epsilon -> 0.")
    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
