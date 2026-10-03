from __future__ import annotations

"""
Numerical validation of Theorem 2.3(iv) from

    P. Zhevandrov et al.,
    "Discrete and embedded trapped modes in a plane quantum waveguide
     with a small obstacle: exact solutions" (March 4, 2025).

Target statement: Theorem 2.3(iv)
---------------------------------
For a small obstacle symmetric with respect to the y-axis,

    X(-t) = -X(t),     Y(-t) = Y(t),     nu = 0,

there exists, for sufficiently small epsilon, a unique embedded eigenvalue in

    Lambda_1 <= k^2 < Lambda_2

provided the vertical displacement is tuned as

    a(epsilon) = epsilon*a_1 + O(epsilon^2),

where

    a_1 = 1/[2(2S + pi*mu)] *
          integral_{-pi}^{pi}
          (Y X' - 3 X Y') (Y - Psi|_Gamma) dt,

and the eigenvalue satisfies

    k^2 = Lambda_2 - sigma^2,
    sigma = epsilon^2 * pi^3 * mu / b^3
            + O(epsilon^3 log epsilon).

Important methodological difference from Theorem 2.3(iii)
---------------------------------------------------------
The y-axis symmetry is a symmetry in x.  As Remark 2.5 of the paper notes,
restricting to even/odd functions in x does NOT remove the continuous spectrum
on [Lambda_1,Lambda_2).  Therefore this script does not replace the full BEM
operator by an A_odd sector.

Instead it:

  1. builds the FULL BEM matrix A(k,a) with fixed M=24;
  2. computes mu, nu, Psi and a_1 independently from the inflated obstacle
     using an exterior-Laplace Method of Fundamental Solutions (MFS);
  3. starts from the theorem predictions
         a ~ epsilon*a_1,
         sigma ~ epsilon^2*pi^3*mu/b^3,
     and uses a coarse SVD minimization ONLY to obtain a seed;
  4. solves the two real equations
         Re(lambda_0(A(k,a))) = 0,
         Im(lambda_0(A(k,a))) = 0,
     where lambda_0 is the eigenvalue of the full BEM matrix closest to zero;
     continuation in epsilon is used to remain on the same branch;
  5. uses fixed-M SVD only as an independent local certification of that root;
  6. reconstructs the FULL field and verifies:
         - wall Dirichlet values,
         - off-grid boundary integral equation,
         - evenness under x -> -x (expected by the proof),
         - suppression of the first open transverse channel,
         - exponential decay with rate sigma;
  7. optionally scans the entire real interval Lambda_1 < k^2 < Lambda_2 at
     the tuned placement to look for additional non-radiating singular minima;
  8. checks both asymptotic laws
         sigma/epsilon^2 -> pi^3*mu/b^3
     and
         (a - epsilon*a_1)/epsilon^2 = O(1);
  9. keeps M fixed at 24 everywhere.  The numerical sensitivity study varies
     only finite-difference step, lattice terms and harmonic order.

The paper also gives, in Example 5.4 for a slightly perturbed circle,

    a_1 = -beta/12 + O(beta^2).

This script does NOT hard-code that approximation.  It evaluates the theorem's
formula (2.7) numerically from Psi.  The value -beta/12 is written to the output
only as a separate paper-example diagnostic.

This is a numerical validation, not a proof of theorem-wide existence,
uniqueness, analyticity, or Big-O bounds.
"""

import csv
import math
import os
import sys
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

# Avoid nested BLAS oversubscription on macOS / worker environments.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLAS_NUM_THREADS",
):
    os.environ.setdefault(_thread_env, "1")
os.environ.setdefault("OMP_DYNAMIC", "FALSE")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from numpy.polynomial.legendre import leggauss  # noqa: E402
from scipy.optimize import least_squares, minimize, minimize_scalar  # noqa: E402

# Same project convention as the existing theorem-validation scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice  # noqa: E402

PI = np.pi


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Config:
    # Waveguide.
    b: float = 1.0

    # Reference y-axis-symmetric obstacle:
    #   X(t) = sin(t) - beta/2 sin(2t)      (odd)
    #   Y(t) = -cos(t) + beta/2 cos(2t)     (even)
    #
    # The physical obstacle is
    #   x = epsilon X(t),
    #   y = a(epsilon) + epsilon Y(t).
    shape_beta: float = 0.50

    epsilon_values: tuple[float, ...] = (
        0.02, 0.04, 0.06, 0.08, 0.10,
        0.12, 0.14, 0.16, 0.18, 0.20,
    )
    publication_asymptotic_epsilon_max: float = 0.10
    publication_min_asymptotic_points: int = 4

    # Single BEM order everywhere.
    bem_M: int = 24

    # Green function / BEM.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Independent MFS calculation of Psi, mu, nu and a1.
    mfs_orders: tuple[int, ...] = (60, 80, 120, 160)
    mfs_source_scale: float = 0.65
    mfs_rcond: float = 1.0e-12
    a1_quadrature_points: int = 4096
    mfs_circle_mu_tolerance: float = 1.0e-6
    mfs_circle_a1_tolerance: float = 1.0e-8
    mfs_convergence_relative_tolerance: float = 5.0e-6
    run_small_beta_a1_diagnostic: bool = True
    small_beta_a1_values: tuple[float, ...] = (0.02, 0.05, 0.10)
    small_beta_a1_order: int = 120

    # Search in the natural spectral variable
    #   sigma = sqrt(Lambda_2-k^2).
    sigma_search_factor: float = 5.0
    sigma_floor_fraction_of_band: float = 1.0e-8

    # The theorem gives a = epsilon*a1 + O(epsilon^2).
    # We optimize the dimensionless correction c in
    #
    #   a = epsilon*a1 + c*epsilon^2.
    #
    # If the optimum lands too close to the boundary, progressively larger
    # O(epsilon^2) windows are attempted.
    placement_correction_factors: tuple[float, ...] = (2.0, 4.0, 8.0)
    placement_absolute_halfwidth_cap: float = 0.20
    placement_boundary_fraction: float = 0.92

    # Stage 1: coarse SVD minimization used ONLY as a root seed.
    coarse_seed_maxiter: int = 70
    coarse_seed_xtol: float = 1.0e-6
    coarse_seed_ftol: float = 1.0e-8

    # Stage 2: solve Re(lambda_0)=Im(lambda_0)=0 for the eigenvalue of A
    # closest to zero. Variables are c and log(sigma/sigma_asym), where
    # a=epsilon*a1+c*epsilon^2.
    root_max_nfev: int = 80
    root_xtol: float = 1.0e-11
    root_ftol: float = 1.0e-11
    root_gtol: float = 1.0e-11
    root_residual_tolerance: float = 2.0e-8
    root_multistart_c_offsets: tuple[float, ...] = (0.0, -0.5, 0.5)
    root_multistart_log_sigma_offsets: tuple[float, ...] = (0.0, -0.08, 0.08)

    # Continuation in epsilon: the previous BIC supplies the first seed for
    # the next epsilon, expressed in the scaled variables c and sigma/sigma_asym.
    use_epsilon_continuation: bool = True

    # Scalar sigma polishing remains available for diagnostic scans and the
    # implementation-sensitivity study, but is NOT used to move the final BIC
    # away from the complex-eigenvalue root.
    polish_sigma_factor: float = 1.35
    minimizer_xatol: float = 1.0e-12

    # Spectral certification at fixed M.
    relative_near_singular_tolerance: float = 1.0e-4
    minimum_drop_factor: float = 100.0
    drop_probe_relative_sigma: float = 0.08

    # Physical diagnostics.
    run_physical_diagnostics: bool = True
    physical_quadrature_points: int = 96
    physical_x_over_b: tuple[float, ...] = (1.5, 2.0, 2.5, 3.0)
    physical_boundary_residual_samples: int = 128
    # Oversample ONLY the diagnostic boundary quadrature by Fourier interpolation.
    # The eigenproblem itself remains M=24.
    physical_boundary_source_oversample_factor: int = 4
    physical_channel_phase_consistency_tolerance: float = 5.0e-2
    physical_wall_relative_tolerance: float = 1.0e-6
    physical_x_even_relative_tolerance: float = 1.0e-5
    physical_open_channel_relative_tolerance: float = 1.0e-5
    physical_boundary_residual_tolerance: float = 5.0e-3
    physical_decay_relative_tolerance: float = 0.20

    # Whole-real-band uniqueness screen.  This is intentionally restricted to
    # the small-epsilon publication window by default because Theorem 2.3 is a
    # small-obstacle theorem and the scan is expensive.
    run_uniqueness_screen: bool = True
    uniqueness_epsilon_max: float = 0.10
    uniqueness_linear_points: int = 44
    uniqueness_log_points: int = 52
    uniqueness_lower_k_margin: float = 2.0e-4
    uniqueness_upper_k_margin: float = 1.0e-7
    uniqueness_exclusion_k_radius: float = 2.0e-3
    uniqueness_max_candidates: int = 10

    # Geometry/symmetry checks.
    geometry_symmetry_tolerance: float = 1.0e-12
    nu_tolerance: float = 1.0e-8

    # Representative implementation sensitivity study: no M sweep.
    run_internal_convergence_study: bool = True
    internal_convergence_epsilon: float = 0.10
    finite_difference_steps_test: tuple[float, ...] = (
        1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4,
    )
    lattice_terms_test: tuple[int, ...] = (25, 50, 100, 200)
    harmonic_orders_test: tuple[int, ...] = (5, 10, 15, 20, 30, 40)

    plot_dpi: int = 220
    output_directory: str = "theorem_2_3_iv_embedded_validation"


CONFIG = Config()


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class ShapeDiagnostics:
    beta: float
    area: float
    max_abs_X: float
    max_abs_Y: float
    min_reference_speed: float
    x_odd_defect: float
    y_even_defect: float


@dataclass
class MFSRow:
    beta: float
    order: int
    source_scale: float
    area: float
    mu: float
    nu: float
    a1_formula_2_7: float
    a1_formula_5_47: float
    a1_internal_relative_difference: float
    a1_paper_example_leading: float
    boundary_relative_residual: float
    zero_total_source_error: float


@dataclass
class TunedMode:
    epsilon: float
    a1: float
    a_leading: float
    a_numerical: float
    scaled_a_correction: float
    sigma_asymptotic: float
    sigma_numerical: float
    kb_numerical: float
    sigma_min: float
    sigma_max: float
    relative_singular_value: float
    drop_factor: float
    search_correction_bound: float
    search_minimum_is_interior: bool
    resolved: bool
    leading_point_relative_singular_value: float
    root_success: bool
    root_residual_norm: float
    root_eigenvalue_real: float
    root_eigenvalue_imag: float
    root_eigenvalue_abs: float
    root_nfev: int
    root_seed_source: str
    continuation_seed_used: bool


@dataclass
class PhysicalDiagnostics:
    epsilon: float
    a: float
    kb: float
    full_M: int
    boundary_x_even_residual: float
    field_x_even_residual: float
    first_open_channel_relative_amplitude: float
    open_channel_left_real: float
    open_channel_left_imag: float
    open_channel_right_real: float
    open_channel_right_imag: float
    open_channel_phase_consistency: float
    wall_relative_residual: float
    boundary_integral_relative_residual: float
    decay_rate_left: float
    decay_rate_right: float
    expected_decay_rate: float
    decay_relative_error_left: float
    decay_relative_error_right: float
    monotone_decay_left: bool
    monotone_decay_right: bool
    boundary_x_even_verified: bool
    field_x_even_verified: bool
    open_channel_suppressed: bool
    open_channel_phase_consistent: bool
    walls_verified: bool
    boundary_integral_verified: bool
    decay_verified: bool


@dataclass
class AdditionalCandidate:
    epsilon: float
    a: float
    kb: float
    sigma: float
    relative_singular_value: float
    drop_factor: float
    boundary_even_residual: float
    open_channel_relative_amplitude: float
    nonradiating_resolved_candidate: bool


@dataclass
class ValidationResult:
    epsilon: float
    beta: float
    area: float
    mu: float
    nu: float
    a1_formula_2_7: float
    a1_paper_example_leading: float
    a_leading: float
    a_numerical: float
    scaled_a_correction: float
    placement_search_bound: float
    placement_minimum_is_interior: bool
    lambda_1: float
    lambda_2: float
    kb_asymptotic: float
    kb_numerical: float
    sigma_asymptotic: float
    sigma_numerical: float
    sigma_over_epsilon_squared: float
    asymptotic_coefficient: float
    scaled_sigma_remainder: float
    relative_error_sigma: float
    leading_point_relative_singular_value: float
    sigma_min_final: float
    relative_singular_value_final: float
    final_drop_factor: float
    root_success: bool
    root_residual_norm: float
    root_eigenvalue_abs: float
    root_seed_source: str
    continuation_seed_used: bool
    expected_bic_resolved: bool
    physical_bic_verified: bool | None
    uniqueness_screen_run: bool
    additional_resolved_bics: int
    uniqueness_screen_passed: bool | None


@dataclass
class InternalConvergenceRow:
    epsilon: float
    parameter: str
    value: float
    kb: float
    sigma_bem: float
    relative_singular_value: float
    relative_kb_shift_from_reference: float
    relative_sigma_shift_from_reference: float


@dataclass
class PublicationClaimRow:
    claim: str
    status: str
    evidence: str
    limitation: str


# ---------------------------------------------------------------------------
# Reference obstacle, theorem quantities and MFS
# ---------------------------------------------------------------------------


def reference_obstacle_geometry(
    t: np.ndarray | float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Inflated y-axis-symmetric obstacle:

        X(t) = sin t - beta/2 sin(2t)      odd,
        Y(t) = -cos t + beta/2 cos(2t)     even.
    """
    t = np.asarray(t)
    beta = float(config.shape_beta)

    X = np.sin(t) - 0.5 * beta * np.sin(2.0 * t)
    Y = -np.cos(t) + 0.5 * beta * np.cos(2.0 * t)

    Xp = np.cos(t) - beta * np.cos(2.0 * t)
    Yp = np.sin(t) - beta * np.sin(2.0 * t)

    Xpp = -np.sin(t) + 2.0 * beta * np.sin(2.0 * t)
    Ypp = np.cos(t) - 2.0 * beta * np.cos(2.0 * t)

    return X, Y, Xp, Yp, Xpp, Ypp


def obstacle_geometry(
    t: np.ndarray | float,
    epsilon: float,
    a: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Physical obstacle in global waveguide coordinates."""
    X, Y, Xp, Yp, Xpp, Ypp = reference_obstacle_geometry(t, config)
    eps = float(epsilon)
    return (
        eps * X,
        float(a) + eps * Y,
        eps * Xp,
        eps * Yp,
        eps * Xpp,
        eps * Ypp,
    )


def lambda_1(config: Config) -> float:
    return PI**2 / (4.0 * config.b**2)


def lambda_2(config: Config) -> float:
    return PI**2 / config.b**2


def kb_1(config: Config) -> float:
    return math.sqrt(lambda_1(config)) * config.b


def kb_2(config: Config) -> float:
    return math.sqrt(lambda_2(config)) * config.b


def asymptotic_coefficient(mu: float, config: Config) -> float:
    return PI**3 * float(mu) / config.b**3


def sigma_from_kb(kb: float, config: Config) -> float:
    cutoff = kb_2(config)
    if kb >= cutoff:
        return 0.0
    gap_sq_dimensionless = max((cutoff - kb) * (cutoff + kb), 0.0)
    return math.sqrt(gap_sq_dimensionless) / config.b


def kb_from_sigma(sigma: float, config: Config) -> float:
    sigma_b = float(sigma) * config.b
    cutoff = kb_2(config)
    radicand = cutoff**2 - sigma_b**2
    if radicand <= 0.0:
        return 0.0
    root = math.sqrt(radicand)
    # Stable subtraction close to the cutoff.
    gap = sigma_b**2 / max(cutoff + root, np.finfo(float).tiny)
    return cutoff - gap


def sigma_band(config: Config) -> tuple[float, float]:
    """Open interval in sigma corresponding to Lambda_1 < k^2 < Lambda_2."""
    k_right = kb_2(config) - config.uniqueness_upper_k_margin
    k_left = kb_1(config) + config.uniqueness_lower_k_margin
    sigma_min = sigma_from_kb(k_right, config)
    sigma_max = sigma_from_kb(k_left, config)
    sigma_min = max(
        sigma_min,
        config.sigma_floor_fraction_of_band * max(sigma_max, 1.0),
    )
    return float(sigma_min), float(sigma_max)


def asymptotic_prediction(
    epsilon: float,
    mu: float,
    config: Config,
) -> tuple[float, float]:
    sigma = asymptotic_coefficient(mu, config) * float(epsilon) ** 2
    k2 = lambda_2(config) - sigma**2
    if not (lambda_1(config) < k2 < lambda_2(config)):
        raise ValueError(
            "The leading Theorem 2.3(iv) prediction lies outside "
            "(Lambda_1,Lambda_2). Use smaller epsilon for theorem-side validation."
        )
    return config.b * math.sqrt(k2), float(sigma)


def compute_shape_diagnostics(config: Config) -> ShapeDiagnostics:
    t = 2.0 * PI * np.arange(8192, dtype=float) / 8192
    X, Y, Xp, Yp, *_ = reference_obstacle_geometry(t, config)
    speed = np.hypot(Xp, Yp)
    dt = 2.0 * PI / len(t)
    area = 0.5 * dt * np.sum(X * Yp - Y * Xp)

    Xm, Ym, *_ = reference_obstacle_geometry(-t, config)
    x_scale = max(float(np.max(np.abs(X))), 1.0)
    y_scale = max(float(np.max(np.abs(Y))), 1.0)
    x_odd = float(np.max(np.abs(Xm + X)) / x_scale)
    y_even = float(np.max(np.abs(Ym - Y)) / y_scale)

    return ShapeDiagnostics(
        beta=float(config.shape_beta),
        area=float(abs(area)),
        max_abs_X=float(np.max(np.abs(X))),
        max_abs_Y=float(np.max(np.abs(Y))),
        min_reference_speed=float(np.min(speed)),
        x_odd_defect=x_odd,
        y_even_defect=y_even,
    )


def _solve_reference_mfs(
    config: Config,
    order: int,
    *,
    beta_override: float | None = None,
) -> tuple[MFSRow, np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the exterior Neumann problem of paper equation (2.1):

        Delta Psi = 0,
        dPsi/dn = n_2,
        grad Psi -> 0.

    The MFS representation uses
        G=-(1/2pi) log|x-s|
    with sources inside the obstacle and an exact zero-total-source constraint.

    For sum(c_j)=0 the MFS solution tends to zero at infinity, matching the
    normalization const=0 used after paper equation (2.2).
    """
    cfg = config if beta_override is None else replace(config, shape_beta=beta_override)
    N = int(order)
    if N < 8:
        raise ValueError("MFS order must be >= 8.")

    tc = 2.0 * PI * (np.arange(N, dtype=float) + 0.5) / N
    X, Y, Xp, Yp, *_ = reference_obstacle_geometry(tc, cfg)
    speed = np.hypot(Xp, Yp)
    if float(np.min(speed)) <= 1.0e-10:
        raise RuntimeError("Reference obstacle is not regular enough for MFS.")

    # CCW parametrization; paper's inward-looking normal:
    #   n=(-Y',X')/|r'|,   n_2=X'/|r'|.
    nx = -Yp / speed
    ny = Xp / speed
    g = ny.copy()

    ts = 2.0 * PI * np.arange(N, dtype=float) / N
    Xs, Ys, *_ = reference_obstacle_geometry(ts, cfg)
    rho = float(cfg.mfs_source_scale)
    sx = rho * Xs
    sy = rho * Ys

    dx = X[:, None] - sx[None, :]
    dy = Y[:, None] - sy[None, :]
    r2 = dx * dx + dy * dy
    if np.any(r2 <= np.finfo(float).tiny):
        raise RuntimeError("MFS source collided with boundary collocation point.")

    A = -(dx * nx[:, None] + dy * ny[:, None]) / (2.0 * PI * r2)

    # Enforce sum c=0 exactly by eliminating the final source strength.
    B = A[:, :-1] - A[:, -1][:, None]
    d, *_ = np.linalg.lstsq(B, g, rcond=cfg.mfs_rcond)
    c = np.concatenate((d, np.array([-float(np.sum(d))])))

    residual = np.linalg.norm(A @ c - g) / max(np.linalg.norm(g), 1.0e-30)
    total_source = float(abs(np.sum(c)))

    # Far-field dipoles:
    #   Psi ~ -nu*x/r^2 - mu*y/r^2.
    nu = -float(np.dot(c, sx)) / (2.0 * PI)
    mu = -float(np.dot(c, sy)) / (2.0 * PI)

    # Evaluate formula (2.7) directly on a fine periodic grid.
    Q = int(cfg.a1_quadrature_points)
    tq = 2.0 * PI * (np.arange(Q, dtype=float) + 0.5) / Q
    Xq, Yq, Xpq, Ypq, *_ = reference_obstacle_geometry(tq, cfg)
    dxq = Xq[:, None] - sx[None, :]
    dyq = Yq[:, None] - sy[None, :]
    rq = np.sqrt(dxq * dxq + dyq * dyq)
    if np.any(rq <= np.finfo(float).tiny):
        raise RuntimeError("MFS source collided with a1 quadrature point.")

    Psi = -(np.log(rq) @ c) / (2.0 * PI)
    dt = 2.0 * PI / Q
    area = 0.5 * dt * np.sum(Xq * Ypq - Yq * Xpq)

    # Formula (2.7).
    numerator_27 = dt * np.sum(
        (Yq * Xpq - 3.0 * Xq * Ypq) * (Yq - Psi)
    )
    denominator_27 = 2.0 * (2.0 * area + PI * mu)
    if abs(denominator_27) <= 1.0e-14:
        raise RuntimeError("Degenerate denominator in Theorem 2.3(iv) a1 formula.")
    a1_27 = float(numerator_27 / denominator_27)

    # Independent algebraic form (5.47), using the paper identity
    #   L0 Y = 1/2 (Y - Psi|_Gamma).
    # This MUST agree with (2.7) if the MFS normalization/sign conventions are
    # internally consistent.  It does not use the Example 5.4 -beta/12
    # approximation.
    L0Y = 0.5 * (Yq - Psi)
    numerator_547 = dt * np.sum(
        Yq * Xpq * L0Y - 3.0 * Xq * Ypq * L0Y
    )
    denominator_547 = 2.0 * area + PI * mu
    a1_547 = float(numerator_547 / denominator_547)
    a1_internal_relative_difference = abs(a1_27 - a1_547) / max(
        abs(a1_27), abs(a1_547), 1.0e-30
    )

    row = MFSRow(
        beta=float(cfg.shape_beta),
        order=N,
        source_scale=rho,
        area=float(area),
        mu=float(mu),
        nu=float(nu),
        a1_formula_2_7=a1_27,
        a1_formula_5_47=a1_547,
        a1_internal_relative_difference=float(a1_internal_relative_difference),
        a1_paper_example_leading=float(-cfg.shape_beta / 12.0),
        boundary_relative_residual=float(residual),
        zero_total_source_error=total_source,
    )
    return row, c, sx, sy


def run_mfs_convergence(
    config: Config,
) -> tuple[float, float, float, list[MFSRow], MFSRow]:
    # Circle calibration: mu=1, nu=0, a1=0.
    circle, *_ = _solve_reference_mfs(
        config,
        max(config.mfs_orders),
        beta_override=0.0,
    )
    if abs(circle.mu - 1.0) > config.mfs_circle_mu_tolerance:
        raise RuntimeError(
            f"MFS circle calibration failed: mu={circle.mu:.12g}, expected 1."
        )
    if abs(circle.a1_formula_2_7) > config.mfs_circle_a1_tolerance:
        raise RuntimeError(
            f"MFS circle a1 calibration failed: a1={circle.a1_formula_2_7:.12g}, expected 0."
        )
    if circle.a1_internal_relative_difference > 1.0e-8:
        raise RuntimeError(
            "MFS a1 formulas (2.7) and (5.47) are internally inconsistent on the circle."
        )

    rows = [_solve_reference_mfs(config, N)[0] for N in config.mfs_orders]
    final = rows[-1]

    if len(rows) >= 2:
        mu_change = abs(rows[-1].mu - rows[-2].mu) / max(abs(rows[-1].mu), 1e-30)
        a1_change = abs(rows[-1].a1_formula_2_7 - rows[-2].a1_formula_2_7) / max(
            abs(rows[-1].a1_formula_2_7), 1e-30
        )
        if mu_change > config.mfs_convergence_relative_tolerance:
            print(
                "WARNING: MFS mu last-step relative change "
                f"{mu_change:.3e} exceeds configured tolerance."
            )
        if a1_change > config.mfs_convergence_relative_tolerance:
            print(
                "WARNING: MFS a1 last-step relative change "
                f"{a1_change:.3e} exceeds configured tolerance."
            )

    if final.mu <= 0.0:
        raise RuntimeError(f"Theorem requires mu>0; MFS returned {final.mu}.")

    return (
        float(final.mu),
        float(final.nu),
        float(final.a1_formula_2_7),
        rows,
        circle,
    )


def geometry_admissibility(
    epsilon: float,
    a: float,
    config: Config,
) -> tuple[bool, str]:
    if epsilon <= 0.0:
        return False, "epsilon must be positive."

    t = 2.0 * PI * np.arange(4096, dtype=float) / 4096
    X, Y, *_ = obstacle_geometry(t, epsilon, a, config)
    if float(np.max(np.abs(Y))) >= config.b:
        return False, "Obstacle touches/intersects a waveguide wall."

    # Conservative radius bound for the local cylindrical Green representation.
    max_dx = 2.0 * float(np.max(np.abs(X)))
    shifted = np.asarray(Y, dtype=float) + config.b
    max_image_dy = 2.0 * float(np.max(np.abs(shifted)))
    max_r = math.hypot(max_dx, max_image_dy)
    green_limit = 0.99 * (4.0 * config.b)
    if max_r > green_limit:
        return False, (
            f"Conservative Green-series radius bound exceeded: "
            f"{max_r:.6g} > {green_limit:.6g}."
        )
    return True, "ok"


# ---------------------------------------------------------------------------
# Full-strip BEM
# ---------------------------------------------------------------------------


def boundary_nodes(M: int) -> np.ndarray:
    return (np.arange(M, dtype=float) + 0.5) * (2.0 * PI / M)


def make_green_functions(
    kb: complex,
    config: Config,
) -> tuple[Callable[..., complex], Callable[..., complex]]:
    b = float(config.b)
    d = 2.0 * b
    k = complex(kb / b)
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
            y + b,
            xi,
            eta + b,
            coefficients,
            k,
            d,
        )

    def green_regularized(x: float, y: float, xi: float, eta: float) -> complex:
        return lattice.greens_dirichlet_reg(
            x,
            y + b,
            xi,
            eta + b,
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
    a: float,
    config: Config,
    G: Callable[..., complex],
    G_regularized: Callable[..., complex],
) -> complex:
    x, y, *_ = obstacle_geometry(psi, epsilon, a, config)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = obstacle_geometry(
        theta, epsilon, a, config
    )

    x = float(x)
    y = float(y)
    xi = float(xi)
    eta = float(eta)
    xi_p = float(xi_p)
    eta_p = float(eta_p)
    xi_pp = float(xi_pp)
    eta_pp = float(eta_pp)

    w = math.hypot(xi_p, eta_p)
    h = config.finite_difference_step

    periodic_distance = abs(math.atan2(math.sin(psi - theta), math.cos(psi - theta)))
    if periodic_distance > 1.0e-12:
        G_xi, G_eta = source_derivatives(G, x, y, xi, eta, h)
        return xi_p * G_eta - eta_p * G_xi

    G_xi_reg, G_eta_reg = source_derivatives(
        G_regularized, x, y, xi, eta, h
    )
    geometric_term = (xi_pp * eta_p - eta_pp * xi_p) / (4.0 * PI * w**2)
    regularized_term = xi_p * G_eta_reg - eta_p * G_xi_reg
    return geometric_term + regularized_term


def assemble_full_matrix(
    kb: complex,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> np.ndarray:
    theta = boundary_nodes(M)
    G, G_regularized = make_green_functions(kb, config)
    K = np.empty((M, M), dtype=np.complex128)

    for i, psi in enumerate(theta):
        for j, source_theta in enumerate(theta):
            K[i, j] = weighted_normal_kernel(
                float(psi),
                float(source_theta),
                epsilon,
                a,
                config,
                G,
                G_regularized,
            )

    # 1/2 u = integral u dG/dn ds,
    # Nyström on the full [0,2pi) contour.
    return np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K


def full_singular_metrics(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[float, float, float]:
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    s = np.linalg.svd(A, compute_uv=False)
    smax = float(s[0])
    smin = float(s[-1])
    return smin, smax, smin / max(smax, np.finfo(float).tiny)


def full_singular_pair(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, np.ndarray]:
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    _, s, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= max(float(np.linalg.norm(vector)), 1.0e-30)

    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))

    smax = float(s[0])
    smin = float(s[-1])
    return smin, smax, smin / max(smax, np.finfo(float).tiny), vector


def relative_sv_from_sigma(
    sigma: float,
    epsilon: float,
    a: float,
    config: Config,
) -> float:
    kb = kb_from_sigma(float(sigma), config)
    if not (kb_1(config) < kb < kb_2(config)):
        return 1.0
    try:
        _, _, rel = full_singular_metrics(
            kb, epsilon, a, config.bem_M, config
        )
    except Exception:
        return 1.0
    return float(rel)


def boundary_x_even_residual(vector: np.ndarray) -> float:
    """
    For this parametrization t -> -t reflects x -> -x.
    On midpoint nodes this is j -> M-1-j.
    The Section 5.2 construction has theta(t) even.
    """
    reflected = vector[::-1]
    return float(
        np.linalg.norm(vector - reflected)
        / max(np.linalg.norm(vector), 1.0e-30)
    )


def near_zero_eigenpair(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[complex, np.ndarray, float]:
    """Return the eigenvalue of A closest to zero and its right eigenvector.

    A true BIC at real (k,a) makes the full outgoing BEM operator singular.
    For a simple root, one eigenvalue lambda_0(A) crosses zero in C.  Solving
    Re(lambda_0)=Im(lambda_0)=0 uses exactly the two real parameters available
    here: the O(epsilon^2) placement correction and sigma.
    """
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    eigvals, eigvecs = np.linalg.eig(A)
    idx = int(np.argmin(np.abs(eigvals)))
    lam = complex(eigvals[idx])
    vec = np.asarray(eigvecs[:, idx], dtype=np.complex128)
    vec /= max(float(np.linalg.norm(vec)), 1.0e-30)
    pivot = int(np.argmax(np.abs(vec)))
    if abs(vec[pivot]) > 0.0:
        vec *= np.exp(-1j * np.angle(vec[pivot]))
    smax = float(np.linalg.svd(A, compute_uv=False)[0])
    return lam, vec, smax


def _root_residual(
    x: np.ndarray,
    epsilon: float,
    a1: float,
    sigma_asym: float,
    config: Config,
) -> np.ndarray:
    a, sigma, kb = _joint_candidate_from_variables(
        epsilon, a1, sigma_asym, float(x[0]), float(x[1]), config
    )
    ok, _ = geometry_admissibility(epsilon, a, config)
    if not ok or not (kb_1(config) < kb < kb_2(config)):
        return np.array([1.0, 1.0], dtype=float)
    try:
        lam, _, smax = near_zero_eigenpair(
            kb, epsilon, a, config.bem_M, config
        )
    except Exception:
        return np.array([1.0, 1.0], dtype=float)
    scale = max(smax, 1.0)
    return np.array([lam.real / scale, lam.imag / scale], dtype=float)


def _variables_from_previous_mode(
    epsilon: float,
    a1: float,
    sigma_asym: float,
    previous_mode: TunedMode,
) -> np.ndarray:
    # Continue the scaled correction c=(a-eps*a1)/eps^2 and the multiplicative
    # sigma/asymptotic ratio, rather than raw a and sigma.
    c = float(previous_mode.scaled_a_correction)
    prev_ratio = previous_mode.sigma_numerical / max(
        previous_mode.sigma_asymptotic, 1.0e-30
    )
    log_ratio = math.log(max(prev_ratio, 1.0e-12))
    return np.array([c, log_ratio], dtype=float)


def _clip_seed(
    seed: np.ndarray,
    c_bound: float,
    log_bound: float,
) -> np.ndarray:
    margin = 1.0e-9
    return np.array([
        float(np.clip(seed[0], -c_bound + margin, c_bound - margin)),
        float(np.clip(seed[1], -log_bound + margin, log_bound - margin)),
    ])


def _solve_complex_eigenvalue_root(
    epsilon: float,
    a1: float,
    sigma_asym: float,
    c_bound: float,
    seeds: list[tuple[str, np.ndarray]],
    config: Config,
) -> tuple[Any, str, complex, float]:
    log_bound = math.log(config.sigma_search_factor)
    best: tuple[Any, str, complex, float] | None = None

    for source, raw_seed in seeds:
        seed = _clip_seed(np.asarray(raw_seed, dtype=float), c_bound, log_bound)
        try:
            result: Any = least_squares(
                lambda x: _root_residual(
                    np.asarray(x, dtype=float),
                    epsilon,
                    a1,
                    sigma_asym,
                    config,
                ),
                x0=seed,
                bounds=(
                    np.array([-c_bound, -log_bound], dtype=float),
                    np.array([ c_bound,  log_bound], dtype=float),
                ),
                xtol=config.root_xtol,
                ftol=config.root_ftol,
                gtol=config.root_gtol,
                max_nfev=config.root_max_nfev,
                x_scale='jac',
            )
        except Exception:
            continue

        a, sigma, kb = _joint_candidate_from_variables(
            epsilon, a1, sigma_asym, float(result.x[0]), float(result.x[1]), config
        )
        try:
            lam, _, smax = near_zero_eigenpair(
                kb, epsilon, a, config.bem_M, config
            )
            root_norm = abs(lam) / max(smax, 1.0)
        except Exception:
            continue

        candidate = (result, source, lam, float(root_norm))
        if best is None or candidate[3] < best[3]:
            best = candidate

        # Most cases should terminate on the continuation or coarse-SVD seed.
        # Extra multistarts are only a fallback, avoiding an unnecessary
        # multiplication of expensive BEM assemblies.
        if root_norm <= config.root_residual_tolerance:
            return candidate

    if best is None:
        raise RuntimeError('Complex-eigenvalue root solve failed for all seeds.')
    return best


# ---------------------------------------------------------------------------
# Joint placement/eigenvalue search
# ---------------------------------------------------------------------------


def _joint_candidate_from_variables(
    epsilon: float,
    a1: float,
    sigma_asym: float,
    correction: float,
    log_sigma_ratio: float,
    config: Config,
) -> tuple[float, float, float]:
    a = epsilon * a1 + correction * epsilon**2
    sigma = sigma_asym * math.exp(log_sigma_ratio)
    smin, smax = sigma_band(config)
    sigma = min(max(sigma, smin), smax)
    return float(a), float(sigma), float(kb_from_sigma(sigma, config))


def _joint_objective(
    x: np.ndarray,
    epsilon: float,
    a1: float,
    sigma_asym: float,
    config: Config,
) -> float:
    a, sigma, kb = _joint_candidate_from_variables(
        epsilon,
        a1,
        sigma_asym,
        float(x[0]),
        float(x[1]),
        config,
    )
    ok, _ = geometry_admissibility(epsilon, a, config)
    if not ok or not (kb_1(config) < kb < kb_2(config)):
        return 4.0
    _, _, rel = full_singular_metrics(kb, epsilon, a, config.bem_M, config)
    return float(math.log10(max(rel, np.finfo(float).tiny)))


def _polish_sigma(
    epsilon: float,
    a: float,
    sigma_seed: float,
    config: Config,
) -> tuple[float, float, float, float]:
    band_lo, band_hi = sigma_band(config)
    lo = max(band_lo, sigma_seed / config.polish_sigma_factor)
    hi = min(band_hi, sigma_seed * config.polish_sigma_factor)
    if not lo < hi:
        lo, hi = band_lo, band_hi

    result: Any = minimize_scalar(
        lambda s: math.log10(
            max(
                relative_sv_from_sigma(float(s), epsilon, a, config),
                np.finfo(float).tiny,
            )
        ),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    sigma = float(result.x) if result.success else float(sigma_seed)
    kb = kb_from_sigma(sigma, config)
    smin, smax, rel = full_singular_metrics(
        kb, epsilon, a, config.bem_M, config
    )
    return sigma, kb, smin, rel



def _spectral_drop_factor(
    epsilon: float,
    a: float,
    sigma: float,
    rel_center: float,
    config: Config,
) -> float:
    band_lo, band_hi = sigma_band(config)
    delta = config.drop_probe_relative_sigma
    left_sigma = max(band_lo, sigma * (1.0 - delta))
    right_sigma = min(band_hi, sigma * (1.0 + delta))
    rel_left = relative_sv_from_sigma(left_sigma, epsilon, a, config)
    rel_right = relative_sv_from_sigma(right_sigma, epsilon, a, config)
    return float(
        min(rel_left, rel_right)
        / max(rel_center, np.finfo(float).tiny)
    )


def tune_embedded_mode(
    epsilon: float,
    mu: float,
    a1: float,
    config: Config,
    previous_mode: TunedMode | None = None,
) -> TunedMode:
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, mu, config)
    a_leading = epsilon * a1

    _, _, leading_rel = full_singular_metrics(
        kb_asym, epsilon, a_leading, config.bem_M, config
    )

    log_bound = math.log(config.sigma_search_factor)
    chosen: tuple[Any, str, complex, float, float, bool] | None = None

    for correction_bound in config.placement_correction_factors:
        cap_in_c_units = config.placement_absolute_halfwidth_cap / max(
            epsilon**2, 1.0e-30
        )
        c_bound = min(float(correction_bound), cap_in_c_units)

        # Coarse SVD seed.  This is discovery only; it is never the final BIC.
        coarse: Any = minimize(
            lambda x: _joint_objective(
                np.asarray(x, dtype=float), epsilon, a1, sigma_asym, config
            ),
            x0=np.array([0.0, 0.0], dtype=float),
            method='Powell',
            bounds=[(-c_bound, c_bound), (-log_bound, log_bound)],
            options={
                'maxiter': config.coarse_seed_maxiter,
                'xtol': config.coarse_seed_xtol,
                'ftol': config.coarse_seed_ftol,
            },
        )

        seeds: list[tuple[str, np.ndarray]] = []
        continuation_used = False
        if config.use_epsilon_continuation and previous_mode is not None:
            seeds.append((
                'epsilon-continuation',
                _variables_from_previous_mode(
                    epsilon, a1, sigma_asym, previous_mode
                ),
            ))
            continuation_used = True

        seeds.append(('coarse-svd', np.asarray(coarse.x, dtype=float)))
        seeds.append(('theorem-leading', np.array([0.0, 0.0], dtype=float)))

        # Small deterministic multi-start cloud.
        for dc in config.root_multistart_c_offsets:
            for ds in config.root_multistart_log_sigma_offsets:
                if dc == 0.0 and ds == 0.0:
                    continue
                seeds.append((
                    f'multistart({dc:+.2g},{ds:+.2g})',
                    np.asarray(coarse.x, dtype=float) + np.array([dc, ds]),
                ))

        result, source, lam, root_norm = _solve_complex_eigenvalue_root(
            epsilon, a1, sigma_asym, c_bound, seeds, config
        )

        c = float(result.x[0])
        logr = float(result.x[1])
        c_interior = abs(c) <= config.placement_boundary_fraction * c_bound
        sigma_interior = abs(logr) <= config.placement_boundary_fraction * log_bound

        chosen = (
            result, source, lam, root_norm, c_bound,
            bool(continuation_used and source == 'epsilon-continuation'),
        )

        # Stop widening once a genuine interior complex root has been found.
        if c_interior and sigma_interior and root_norm <= config.root_residual_tolerance:
            break

    if chosen is None:
        raise RuntimeError('BIC tuning did not produce a complex-eigenvalue root candidate.')

    result, source, lam, root_norm, c_bound, continuation_seed_used = chosen
    correction = float(result.x[0])
    log_ratio = float(result.x[1])
    a, sigma, kb = _joint_candidate_from_variables(
        epsilon, a1, sigma_asym, correction, log_ratio, config
    )

    smin, smax, rel = full_singular_metrics(
        kb, epsilon, a, config.bem_M, config
    )
    scaled_correction = (a - a_leading) / max(epsilon**2, 1.0e-30)
    placement_interior = bool(
        abs(scaled_correction) <= config.placement_boundary_fraction * c_bound
        and abs(log_ratio) <= config.placement_boundary_fraction * log_bound
    )

    drop = _spectral_drop_factor(epsilon, a, sigma, rel, config)
    root_ok = bool(root_norm <= config.root_residual_tolerance)
    resolved = bool(
        placement_interior
        and root_ok
        and rel <= config.relative_near_singular_tolerance
        and drop >= config.minimum_drop_factor
        and kb_1(config) < kb < kb_2(config)
    )

    return TunedMode(
        epsilon=float(epsilon),
        a1=float(a1),
        a_leading=float(a_leading),
        a_numerical=float(a),
        scaled_a_correction=float(scaled_correction),
        sigma_asymptotic=float(sigma_asym),
        sigma_numerical=float(sigma),
        kb_numerical=float(kb),
        sigma_min=float(smin),
        sigma_max=float(smax),
        relative_singular_value=float(rel),
        drop_factor=float(drop),
        search_correction_bound=float(c_bound),
        search_minimum_is_interior=bool(placement_interior),
        resolved=resolved,
        leading_point_relative_singular_value=float(leading_rel),
        root_success=root_ok,
        root_residual_norm=float(root_norm),
        root_eigenvalue_real=float(lam.real),
        root_eigenvalue_imag=float(lam.imag),
        root_eigenvalue_abs=float(abs(lam)),
        root_nfev=int(getattr(result, 'nfev', -1)),
        root_seed_source=str(source),
        continuation_seed_used=bool(continuation_seed_used),
    )


# ---------------------------------------------------------------------------
# Full-field reconstruction and BIC diagnostics
# ---------------------------------------------------------------------------


def make_extended_field_green(
    kb: float,
    config: Config,
) -> Callable[..., complex]:
    b = config.b
    d = 2.0 * b
    period = 2.0 * d
    k = complex(kb / b)
    coeffs = lattice.lattice_sums(
        period,
        k,
        beta=0.0,
        M=config.lattice_terms,
        Lh=config.harmonic_order,
    )

    def wrap(value: float) -> float:
        return float((value + 0.5 * period) % period - 0.5 * period)

    def green(x: float, y: float, xi: float, eta: float) -> complex:
        yf = y + b
        ys = eta + b
        X = x - xi
        Y1 = wrap(yf - ys)
        Y2 = wrap(yf + ys)
        return (
            lattice.greens_periodic(X, Y1, coeffs, k, period)
            - lattice.greens_periodic(X, Y2, coeffs, k, period)
        )

    return green


def weighted_kernel_at_field_point(
    x: float,
    y: float,
    theta: float,
    epsilon: float,
    a: float,
    config: Config,
    G: Callable[..., complex],
) -> complex:
    xi, eta, xi_p, eta_p, *_ = obstacle_geometry(
        theta, epsilon, a, config
    )
    G_xi, G_eta = source_derivatives(
        G,
        float(x),
        float(y),
        float(xi),
        float(eta),
        config.finite_difference_step,
    )
    return float(xi_p) * G_eta - float(eta_p) * G_xi


def reconstruct_field_at_points(
    points: list[tuple[float, float]],
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    epsilon: float,
    a: float,
    config: Config,
    G: Callable[..., complex],
) -> np.ndarray:
    dtheta = 2.0 * PI / len(theta)
    values = np.empty(len(points), dtype=np.complex128)

    for i, (x, y) in enumerate(points):
        kernel = np.array(
            [
                weighted_kernel_at_field_point(
                    float(x),
                    float(y),
                    float(t),
                    epsilon,
                    a,
                    config,
                    G,
                )
                for t in theta
            ],
            dtype=np.complex128,
        )
        values[i] = dtheta * np.dot(kernel, boundary_vector)

    return values


def periodic_fourier_interpolate(
    theta_nodes: np.ndarray,
    values: np.ndarray,
    theta_targets: np.ndarray,
) -> np.ndarray:
    M = len(theta_nodes)
    theta0 = float(theta_nodes[0])
    coeffs = np.fft.fft(values) / M
    frequencies = np.fft.fftfreq(M, d=1.0 / M)
    phase = np.exp(
        1j * np.outer(np.asarray(theta_targets, dtype=float) - theta0, frequencies)
    )
    return phase @ coeffs


def offgrid_boundary_integral_residual(
    epsilon: float,
    a: float,
    kb: float,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    config: Config,
) -> float:
    """Off-grid BIE residual using oversampled DIAGNOSTIC quadrature.

    The spectral problem remains M=config.bem_M.  We only Fourier-interpolate
    the already-computed boundary vector to a denser source grid so this
    diagnostic is not limited by the 24-point trapezoid rule.
    """
    source_count = max(
        int(config.physical_boundary_source_oversample_factor) * len(theta),
        4 * len(theta),
    )
    source_theta = (
        (np.arange(source_count, dtype=float) + 0.5)
        * (2.0 * PI / source_count)
    )
    source_values = periodic_fourier_interpolate(
        theta, boundary_vector, source_theta
    )

    sample_count = max(
        int(config.physical_boundary_residual_samples),
        source_count,
    )
    targets = (
        (np.arange(sample_count, dtype=float) + 0.371)
        * (2.0 * PI / sample_count)
    ) % (2.0 * PI)

    u_targets = periodic_fourier_interpolate(theta, boundary_vector, targets)
    G, Greg = make_green_functions(complex(kb, 0.0), config)
    dtheta = 2.0 * PI / source_count

    residuals = np.empty(sample_count, dtype=np.complex128)
    rhs_values = np.empty(sample_count, dtype=np.complex128)

    for i, psi in enumerate(targets):
        kernel = np.array(
            [
                weighted_normal_kernel(
                    float(psi), float(src), epsilon, a, config, G, Greg
                )
                for src in source_theta
            ],
            dtype=np.complex128,
        )
        rhs = dtheta * np.dot(kernel, source_values)
        rhs_values[i] = rhs
        residuals[i] = 0.5 * u_targets[i] - rhs

    scale = max(
        float(np.max(np.abs(0.5 * u_targets))),
        float(np.max(np.abs(rhs_values))),
        float(np.max(np.abs(boundary_vector))),
        1.0e-30,
    )
    return float(np.max(np.abs(residuals)) / scale)


def first_open_channel_diagnostics(
    epsilon: float,
    a: float,
    kb: float,
    vector: np.ndarray,
    theta: np.ndarray,
    config: Config,
) -> dict[str, Any]:
    """Signed first-open-channel coefficient on several cross-sections.

    For Lambda_1<k^2<Lambda_2 the only open transverse channel is
        phi_1(y)=cos(pi y/(2b)).
    A true BIC must have zero complex coefficient, not merely a small positive
    norm ratio.  We therefore retain the complex coefficients and also test
    whether removing the expected propagating phase makes them cross-section
    independent.
    """
    G = make_extended_field_green(kb, config)
    b = config.b
    distances = np.asarray(config.physical_x_over_b, dtype=float) * b

    nodes, weights = leggauss(config.physical_quadrature_points)
    y_values = b * nodes
    y_weights = b * weights
    phi1 = np.cos(PI * y_values / (2.0 * b))
    norm_phi_sq = float(np.sum(y_weights * phi1**2))

    k = kb / b
    q1 = math.sqrt(max(k * k - lambda_1(config), 0.0))

    coeffs_left: list[complex] = []
    coeffs_right: list[complex] = []
    ratios_left: list[float] = []
    ratios_right: list[float] = []

    for distance in distances:
        for label, x in (("left", -float(distance)), ("right", float(distance))):
            field = reconstruct_field_at_points(
                [(x, float(y)) for y in y_values],
                vector, theta, epsilon, a, config, G,
            )
            coeff = complex(np.sum(y_weights * field * phi1) / norm_phi_sq)
            field_norm = math.sqrt(float(np.sum(y_weights * np.abs(field) ** 2)))
            ratio = abs(coeff) * math.sqrt(norm_phi_sq) / max(field_norm, 1.0e-30)
            if label == "left":
                coeffs_left.append(coeff)
                ratios_left.append(float(ratio))
            else:
                coeffs_right.append(coeff)
                ratios_right.append(float(ratio))

    def phase_consistency(coeffs: list[complex], side: str) -> tuple[float, complex]:
        if not coeffs:
            return math.nan, complex(math.nan, math.nan)
        arr = np.asarray(coeffs, dtype=np.complex128)
        # With the project Green-function convention the outgoing sign can be
        # implementation-dependent. Try both signs and keep the one with the
        # smaller relative cross-section variation. This is diagnostic only.
        candidates: list[tuple[float, np.ndarray]] = []
        for sign in (-1.0, 1.0):
            if side == "right":
                phase = np.exp(sign * 1j * q1 * distances)
            else:
                phase = np.exp(sign * 1j * q1 * distances)
            dephased = arr * phase
            mean = np.mean(dephased)
            variation = float(
                np.max(np.abs(dephased - mean))
                / max(float(np.max(np.abs(dephased))), 1.0e-30)
            )
            candidates.append((variation, dephased))
        variation, dephased = min(candidates, key=lambda item: item[0])
        return variation, complex(np.mean(dephased))

    left_consistency, left_signed = phase_consistency(coeffs_left, "left")
    right_consistency, right_signed = phase_consistency(coeffs_right, "right")

    return {
        "ratio_left": float(max(ratios_left) if ratios_left else math.nan),
        "ratio_right": float(max(ratios_right) if ratios_right else math.nan),
        "coefficient_left": left_signed,
        "coefficient_right": right_signed,
        "phase_consistency": float(max(left_consistency, right_consistency)),
        "raw_coefficients_left": coeffs_left,
        "raw_coefficients_right": coeffs_right,
    }


def physical_mode_diagnostics(
    epsilon: float,
    mode: TunedMode,
    config: Config,
) -> PhysicalDiagnostics:
    a = mode.a_numerical
    kb = mode.kb_numerical
    _, _, _, vector = full_singular_pair(
        kb, epsilon, a, config.bem_M, config
    )
    theta = boundary_nodes(config.bem_M)
    G = make_extended_field_green(kb, config)
    b = config.b

    boundary_even = boundary_x_even_residual(vector)

    nodes, weights = leggauss(config.physical_quadrature_points)
    y_values = b * nodes
    y_weights = b * weights
    distances = np.asarray(config.physical_x_over_b, dtype=float) * b

    left_norms: list[float] = []
    right_norms: list[float] = []
    cross_peak = 0.0

    for distance in distances:
        for label, x in (("left", -float(distance)), ("right", float(distance))):
            field = reconstruct_field_at_points(
                [(x, float(y)) for y in y_values],
                vector, theta, epsilon, a, config, G,
            )
            norm = math.sqrt(float(np.sum(y_weights * np.abs(field) ** 2)))
            if label == "left":
                left_norms.append(norm)
            else:
                right_norms.append(norm)
            cross_peak = max(cross_peak, float(np.max(np.abs(field))))

    left_arr = np.asarray(left_norms, dtype=float)
    right_arr = np.asarray(right_norms, dtype=float)
    monotone_left = bool(np.all(np.diff(left_arr) <= 1.0e-10 * max(left_arr[0], 1.0)))
    monotone_right = bool(np.all(np.diff(right_arr) <= 1.0e-10 * max(right_arr[0], 1.0)))

    def fit_decay(norms: np.ndarray) -> float:
        positive = np.maximum(norms, np.finfo(float).tiny)
        slope, _ = np.polyfit(distances, np.log(positive), 1)
        return float(-slope)

    decay_left = fit_decay(left_arr)
    decay_right = fit_decay(right_arr)
    expected = mode.sigma_numerical
    decay_error_left = abs(decay_left - expected) / max(expected, 1.0e-30)
    decay_error_right = abs(decay_right - expected) / max(expected, 1.0e-30)

    # Wall values.
    x0 = float(distances[0]) if len(distances) else 2.0 * b
    wall_points = [
        (-x0, -b), (-x0, b),
        (0.0, -b), (0.0, b),
        (x0, -b), (x0, b),
    ]
    wall_values = reconstruct_field_at_points(
        wall_points, vector, theta, epsilon, a, config, G
    )

    # Independent full-field x-reflection test.
    parity_points: list[tuple[float, float]] = []
    parity_pairs: list[tuple[int, int]] = []
    for x in (0.7 * x0, x0):
        for frac in (-0.65, -0.25, 0.20, 0.55):
            y = frac * b
            i = len(parity_points)
            parity_points.extend([(-x, y), (x, y)])
            parity_pairs.append((i, i + 1))
    parity_values = reconstruct_field_at_points(
        parity_points, vector, theta, epsilon, a, config, G
    )

    amplitude_scale = max(
        float(np.max(np.abs(vector))), cross_peak, 1.0e-30
    )
    wall_rel = float(np.max(np.abs(wall_values)) / amplitude_scale)
    field_even = 0.0
    for i, j in parity_pairs:
        field_even = max(
            field_even,
            float(abs(parity_values[i] - parity_values[j]) / amplitude_scale),
        )

    channel = first_open_channel_diagnostics(
        epsilon, a, kb, vector, theta, config
    )
    open_rel = max(channel["ratio_left"], channel["ratio_right"])
    left_coeff = complex(channel["coefficient_left"])
    right_coeff = complex(channel["coefficient_right"])
    phase_consistency = float(channel["phase_consistency"])

    bie_rel = offgrid_boundary_integral_residual(
        epsilon, a, kb, vector, theta, config
    )

    return PhysicalDiagnostics(
        epsilon=float(epsilon),
        a=float(a),
        kb=float(kb),
        full_M=int(config.bem_M),
        boundary_x_even_residual=float(boundary_even),
        field_x_even_residual=float(field_even),
        first_open_channel_relative_amplitude=float(open_rel),
        open_channel_left_real=float(left_coeff.real),
        open_channel_left_imag=float(left_coeff.imag),
        open_channel_right_real=float(right_coeff.real),
        open_channel_right_imag=float(right_coeff.imag),
        open_channel_phase_consistency=float(phase_consistency),
        wall_relative_residual=float(wall_rel),
        boundary_integral_relative_residual=float(bie_rel),
        decay_rate_left=float(decay_left),
        decay_rate_right=float(decay_right),
        expected_decay_rate=float(expected),
        decay_relative_error_left=float(decay_error_left),
        decay_relative_error_right=float(decay_error_right),
        monotone_decay_left=monotone_left,
        monotone_decay_right=monotone_right,
        boundary_x_even_verified=(
            boundary_even <= config.physical_x_even_relative_tolerance
        ),
        field_x_even_verified=(
            field_even <= config.physical_x_even_relative_tolerance
        ),
        open_channel_suppressed=(
            open_rel <= config.physical_open_channel_relative_tolerance
        ),
        open_channel_phase_consistent=(
            phase_consistency <= config.physical_channel_phase_consistency_tolerance
            or open_rel <= config.physical_open_channel_relative_tolerance
        ),
        walls_verified=(
            wall_rel <= config.physical_wall_relative_tolerance
        ),
        boundary_integral_verified=(
            bie_rel <= config.physical_boundary_residual_tolerance
        ),
        decay_verified=bool(
            monotone_left
            and monotone_right
            and decay_error_left <= config.physical_decay_relative_tolerance
            and decay_error_right <= config.physical_decay_relative_tolerance
        ),
    )


# ---------------------------------------------------------------------------
# Whole-real-band uniqueness screen at the tuned placement
# ---------------------------------------------------------------------------


def build_uniqueness_grid(
    epsilon: float,
    mu: float,
    config: Config,
) -> np.ndarray:
    k_lo = kb_1(config) + config.uniqueness_lower_k_margin
    k_hi = kb_2(config) - config.uniqueness_upper_k_margin

    linear = np.linspace(
        k_lo,
        k_hi,
        config.uniqueness_linear_points,
    )

    sigma_lo = sigma_from_kb(k_hi, config)
    sigma_hi = sigma_from_kb(k_lo, config)
    sigma_lo = max(sigma_lo, np.finfo(float).eps)
    sigmas = np.geomspace(
        sigma_lo,
        sigma_hi,
        config.uniqueness_log_points,
    )
    logarithmic = np.array(
        [kb_from_sigma(float(s), config) for s in sigmas],
        dtype=float,
    )

    _, sigma_asym = asymptotic_prediction(epsilon, mu, config)
    local_sigmas = sigma_asym * np.array(
        [0.20, 0.35, 0.55, 0.75, 1.0, 1.3, 1.8, 2.6, 4.0],
        dtype=float,
    )
    local_sigmas = local_sigmas[
        (local_sigmas > sigma_lo) & (local_sigmas < sigma_hi)
    ]
    local = np.array(
        [kb_from_sigma(float(s), config) for s in local_sigmas],
        dtype=float,
    )

    grid = np.unique(np.sort(np.concatenate((linear, logarithmic, local))))
    return grid[(grid > k_lo) & (grid < k_hi)]


def sampled_local_minimum_brackets(
    scan: np.ndarray,
    values: np.ndarray,
) -> list[tuple[float, float]]:
    brackets: list[tuple[float, float]] = []
    for i in range(1, len(scan) - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            brackets.append((float(scan[i - 1]), float(scan[i + 1])))
    return brackets


def refine_k_bracket(
    left: float,
    right: float,
    epsilon: float,
    a: float,
    config: Config,
) -> tuple[float, float, float, float, np.ndarray]:
    tiny = np.finfo(float).tiny

    def relsv(kb: float) -> float:
        _, _, rel = full_singular_metrics(
            float(kb), epsilon, a, config.bem_M, config
        )
        return rel

    left_rel = relsv(left)
    right_rel = relsv(right)
    result: Any = minimize_scalar(
        lambda k: math.log10(max(relsv(float(k)), tiny)),
        bounds=(left, right),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    kb = float(result.x)
    smin, _, rel, vector = full_singular_pair(
        kb, epsilon, a, config.bem_M, config
    )
    drop = min(left_rel, right_rel) / max(rel, tiny)
    return kb, sigma_from_kb(kb, config), smin, drop, vector


def whole_band_uniqueness_screen(
    epsilon: float,
    mu: float,
    mode: TunedMode,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, list[AdditionalCandidate]]:
    scan = build_uniqueness_grid(epsilon, mu, config)
    values = np.empty(len(scan), dtype=float)

    print(
        f"  whole-band real-axis screen: {len(scan)} points, "
        f"fixed M={config.bem_M}"
    )

    for i, kb in enumerate(scan):
        _, _, values[i] = full_singular_metrics(
            float(kb),
            epsilon,
            mode.a_numerical,
            config.bem_M,
            config,
        )
        if (i + 1) % 20 == 0 or i + 1 == len(scan):
            print(f"    {i+1:>3}/{len(scan)}")

    brackets = sampled_local_minimum_brackets(scan, values)
    brackets = [
        br for br in brackets
        if not (
            br[0] - config.uniqueness_exclusion_k_radius
            <= mode.kb_numerical
            <= br[1] + config.uniqueness_exclusion_k_radius
        )
    ][: config.uniqueness_max_candidates]

    additional: list[AdditionalCandidate] = []

    for left, right in brackets:
        kb, sigma, _, drop, vector = refine_k_bracket(
            left, right, epsilon, mode.a_numerical, config
        )
        _, _, rel = full_singular_metrics(
            kb,
            epsilon,
            mode.a_numerical,
            config.bem_M,
            config,
        )

        even_res = boundary_x_even_residual(vector)
        theta = boundary_nodes(config.bem_M)
        channel = first_open_channel_diagnostics(
            epsilon,
            mode.a_numerical,
            kb,
            vector,
            theta,
            config,
        )
        open_rel = max(channel["ratio_left"], channel["ratio_right"])

        # Do NOT require x-even parity here: uniqueness concerns any additional
        # embedded eigenvalue, regardless of its x-parity.
        resolved = bool(
            rel <= config.relative_near_singular_tolerance
            and drop >= config.minimum_drop_factor
            and open_rel <= config.physical_open_channel_relative_tolerance
        )

        additional.append(
            AdditionalCandidate(
                epsilon=float(epsilon),
                a=float(mode.a_numerical),
                kb=float(kb),
                sigma=float(sigma),
                relative_singular_value=float(rel),
                drop_factor=float(drop),
                boundary_even_residual=float(even_res),
                first_open_channel_relative_amplitude=float(open_rel),
                nonradiating_resolved_candidate=resolved,
            )
        )

    return scan, values, additional


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_epsilon(
    epsilon: float,
    shape: ShapeDiagnostics,
    mu: float,
    nu: float,
    a1: float,
    config: Config,
    previous_mode: TunedMode | None = None,
) -> tuple[
    ValidationResult,
    TunedMode,
    PhysicalDiagnostics | None,
    np.ndarray,
    np.ndarray,
    list[AdditionalCandidate],
]:
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, mu, config)

    mode = tune_embedded_mode(epsilon, mu, a1, config, previous_mode=previous_mode)
    ok_geometry, reason = geometry_admissibility(
        epsilon, mode.a_numerical, config
    )
    if not ok_geometry:
        raise RuntimeError(
            f"epsilon={epsilon:g}: tuned geometry inadmissible: {reason}"
        )

    physical: PhysicalDiagnostics | None = None
    physical_ok: bool | None = None
    if config.run_physical_diagnostics and mode.resolved:
        physical = physical_mode_diagnostics(epsilon, mode, config)
        physical_ok = bool(
            physical.boundary_x_even_verified
            and physical.field_x_even_verified
            and physical.open_channel_suppressed
            and physical.walls_verified
            and physical.boundary_integral_verified
            and physical.decay_verified
        )

    run_scan = bool(
        config.run_uniqueness_screen
        and epsilon <= config.uniqueness_epsilon_max + 1.0e-15
        and mode.resolved
        and physical_ok is True
    )

    if run_scan:
        scan, scan_values, additional = whole_band_uniqueness_screen(
            epsilon, mu, mode, config
        )
        additional_resolved = sum(
            c.nonradiating_resolved_candidate for c in additional
        )
        uniqueness_passed: bool | None = additional_resolved == 0
    else:
        scan = np.array([], dtype=float)
        scan_values = np.array([], dtype=float)
        additional = []
        additional_resolved = 0
        uniqueness_passed = None

    C = asymptotic_coefficient(mu, config)
    sigma_error = abs(mode.sigma_numerical - sigma_asym) / max(
        abs(sigma_asym), 1.0e-30
    )
    scaled_sigma_remainder = (
        abs(mode.sigma_numerical - C * epsilon**2)
        / max(epsilon**3 * abs(math.log(epsilon)), 1.0e-30)
    )

    result = ValidationResult(
        epsilon=float(epsilon),
        beta=float(config.shape_beta),
        area=float(shape.area),
        mu=float(mu),
        nu=float(nu),
        a1_formula_2_7=float(a1),
        a1_paper_example_leading=float(-config.shape_beta / 12.0),
        a_leading=float(epsilon * a1),
        a_numerical=float(mode.a_numerical),
        scaled_a_correction=float(mode.scaled_a_correction),
        placement_search_bound=float(mode.search_correction_bound),
        placement_minimum_is_interior=bool(mode.search_minimum_is_interior),
        lambda_1=float(lambda_1(config)),
        lambda_2=float(lambda_2(config)),
        kb_asymptotic=float(kb_asym),
        kb_numerical=float(mode.kb_numerical),
        sigma_asymptotic=float(sigma_asym),
        sigma_numerical=float(mode.sigma_numerical),
        sigma_over_epsilon_squared=float(
            mode.sigma_numerical / epsilon**2
        ),
        asymptotic_coefficient=float(C),
        scaled_sigma_remainder=float(scaled_sigma_remainder),
        relative_error_sigma=float(sigma_error),
        leading_point_relative_singular_value=float(
            mode.leading_point_relative_singular_value
        ),
        sigma_min_final=float(mode.sigma_min),
        relative_singular_value_final=float(mode.relative_singular_value),
        final_drop_factor=float(mode.drop_factor),
        root_success=bool(mode.root_success),
        root_residual_norm=float(mode.root_residual_norm),
        root_eigenvalue_abs=float(mode.root_eigenvalue_abs),
        root_seed_source=str(mode.root_seed_source),
        continuation_seed_used=bool(mode.continuation_seed_used),
        expected_bic_resolved=bool(mode.resolved),
        physical_bic_verified=physical_ok,
        uniqueness_screen_run=run_scan,
        additional_resolved_bics=int(additional_resolved),
        uniqueness_screen_passed=uniqueness_passed,
    )

    return result, mode, physical, scan, scan_values, additional


# ---------------------------------------------------------------------------
# Internal convergence at fixed M
# ---------------------------------------------------------------------------


def run_internal_convergence_study(
    summaries: list[ValidationResult],
    a1: float,
    config: Config,
) -> list[InternalConvergenceRow]:
    if not config.run_internal_convergence_study:
        return []

    target = min(
        summaries,
        key=lambda r: abs(r.epsilon - config.internal_convergence_epsilon),
    )
    eps = target.epsilon
    a = target.a_numerical
    sigma_seed = target.sigma_numerical

    def solve(cfg: Config) -> tuple[float, float, float]:
        sigma, kb, _, rel = _polish_sigma(eps, a, sigma_seed, cfg)
        return kb, sigma, rel

    reference_cfg = config
    kb_ref, sigma_ref, _ = solve(reference_cfg)

    rows: list[InternalConvergenceRow] = []

    def append(parameter: str, value: float, cfg: Config) -> None:
        kb, sigma, rel = solve(cfg)
        rows.append(
            InternalConvergenceRow(
                epsilon=float(eps),
                parameter=parameter,
                value=float(value),
                kb=float(kb),
                sigma_bem=float(sigma),
                relative_singular_value=float(rel),
                relative_kb_shift_from_reference=float(
                    abs(kb - kb_ref) / max(abs(kb_ref), 1e-30)
                ),
                relative_sigma_shift_from_reference=float(
                    abs(sigma - sigma_ref) / max(abs(sigma_ref), 1e-30)
                ),
            )
        )

    for h in config.finite_difference_steps_test:
        append(
            "finite_difference_step",
            h,
            replace(config, finite_difference_step=h),
        )

    for value in config.lattice_terms_test:
        append(
            "lattice_terms",
            float(value),
            replace(config, lattice_terms=int(value)),
        )

    for value in config.harmonic_orders_test:
        append(
            "harmonic_order",
            float(value),
            replace(config, harmonic_order=int(value)),
        )

    return rows


# ---------------------------------------------------------------------------
# Output / plots / report
# ---------------------------------------------------------------------------


def write_dataclass_csv(path: Path, rows: Sequence[Any]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0].keys()))
        writer.writeheader()
        writer.writerows(dictionaries)


def plot_shape(config: Config, output: Path) -> None:
    t = np.linspace(0.0, 2.0 * PI, 1200)
    X, Y, *_ = reference_obstacle_geometry(t, config)
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(X, Y)
    plt.axhline(0.0, linewidth=0.8)
    plt.axvline(0.0, linewidth=0.8)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Inflated y-axis-symmetric obstacle")
    plt.tight_layout()
    plt.savefig(output / "inflated_obstacle.png", dpi=config.plot_dpi)
    plt.close()


def plot_mfs_convergence(rows: list[MFSRow], output: Path, config: Config) -> None:
    N = np.array([r.order for r in rows], dtype=float)

    plt.figure(figsize=(7, 4.5))
    plt.plot(N, [r.mu for r in rows], "o-")
    plt.xlabel("MFS order")
    plt.ylabel(r"$\mu$")
    plt.tight_layout()
    plt.savefig(output / "mu_mfs_convergence.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(N, [r.a1_formula_2_7 for r in rows], "o-", label="formula (2.7)")
    plt.plot(N, [r.a1_formula_5_47 for r in rows], "s--", label="formula (5.47)")
    plt.axhline(
        -config.shape_beta / 12.0,
        linestyle="--",
        label=r"Example 5.4: $-\beta/12$",
    )
    plt.xlabel("MFS order")
    plt.ylabel(r"$a_1$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "a1_mfs_convergence.png", dpi=config.plot_dpi)
    plt.close()


def plot_summary(
    results: list[ValidationResult],
    output: Path,
    config: Config,
) -> None:
    eps = np.array([r.epsilon for r in results], dtype=float)
    kb_num = np.array([r.kb_numerical for r in results], dtype=float)
    kb_asym = np.array([r.kb_asymptotic for r in results], dtype=float)
    ratio = np.array([r.sigma_over_epsilon_squared for r in results], dtype=float)
    C = results[0].asymptotic_coefficient
    a_over_eps = np.array([r.a_numerical / r.epsilon for r in results])
    a1 = results[0].a1_formula_2_7
    correction = np.array([r.scaled_a_correction for r in results])

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, kb_num, "o-", label="BEM")
    plt.plot(eps, kb_asym, "s--", label="asymptotic")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$kb$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "kb_vs_epsilon.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, ratio, "o-", label=r"$\sigma_{\rm BEM}/\varepsilon^2$")
    plt.axhline(C, linestyle="--", label=r"$\pi^3\mu/b^3$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$\sigma/\varepsilon^2$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "sigma_over_epsilon_squared.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, a_over_eps, "o-", label=r"$a_{\rm num}/\varepsilon$")
    plt.axhline(a1, linestyle="--", label=r"$a_1$ from (2.7)")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$a/\varepsilon$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "placement_a_over_epsilon.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, correction, "o-")
    plt.axhline(0.0, linewidth=0.8)
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$(a-\varepsilon a_1)/\varepsilon^2$")
    plt.tight_layout()
    plt.savefig(output / "scaled_placement_correction.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(
        eps,
        [max(r.relative_singular_value_final, 1e-30) for r in results],
        "o-",
    )
    plt.axhline(
        config.relative_near_singular_tolerance,
        linestyle="--",
        label="spectral tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$s_{\min}/s_{\max}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "relative_singular_value.png", dpi=config.plot_dpi)
    plt.close()


def plot_uniqueness_scan(
    epsilon: float,
    scan: np.ndarray,
    values: np.ndarray,
    result: ValidationResult,
    output: Path,
    config: Config,
) -> None:
    if len(scan) == 0:
        return
    plt.figure(figsize=(8, 5))
    plt.semilogy(scan, values, "o-", markersize=3)
    plt.axvline(result.kb_numerical, linestyle="--", label="main BIC")
    plt.xlabel(r"$kb$")
    plt.ylabel(r"$s_{\min}/s_{\max}$")
    plt.title(rf"Whole-band BIC screen, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        output / f"whole_band_epsilon_{epsilon:.2f}.png",
        dpi=config.plot_dpi,
    )
    plt.close()


def plot_internal_convergence(
    rows: list[InternalConvergenceRow],
    output: Path,
    config: Config,
) -> None:
    for parameter in (
        "finite_difference_step",
        "lattice_terms",
        "harmonic_order",
    ):
        subset = [r for r in rows if r.parameter == parameter]
        if not subset:
            continue
        plt.figure(figsize=(7, 4.5))
        plt.plot(
            [r.value for r in subset],
            [r.relative_sigma_shift_from_reference for r in subset],
            "o-",
        )
        if parameter == "finite_difference_step":
            plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(parameter)
        plt.ylabel("relative sigma shift from reference")
        plt.tight_layout()
        plt.savefig(
            output / f"internal_convergence_{parameter}.png",
            dpi=config.plot_dpi,
        )
        plt.close()


def build_publication_claims(
    shape: ShapeDiagnostics,
    mu: float,
    nu: float,
    a1: float,
    summaries: list[ValidationResult],
    config: Config,
) -> list[PublicationClaimRow]:
    claims: list[PublicationClaimRow] = []

    symmetry_ok = bool(
        shape.x_odd_defect <= config.geometry_symmetry_tolerance
        and shape.y_even_defect <= config.geometry_symmetry_tolerance
        and abs(nu) <= config.nu_tolerance
    )
    claims.append(
        PublicationClaimRow(
            claim="Theorem 2.3(iv) symmetry assumptions: X odd, Y even, nu=0",
            status="supported" if symmetry_ok else "partial",
            evidence=(
                f"X-odd defect={shape.x_odd_defect:.3e}, "
                f"Y-even defect={shape.y_even_defect:.3e}, nu={nu:.3e}."
            ),
            limitation="Numerical symmetry/dipole checks do not replace the analytic symmetry argument.",
        )
    )

    small = [
        r for r in summaries
        if r.epsilon <= config.publication_asymptotic_epsilon_max + 1e-15
    ]

    placement_ok = bool(
        len(small) >= config.publication_min_asymptotic_points
        and all(r.expected_bic_resolved for r in small)
        and all(r.physical_bic_verified is True for r in small)
        and all(r.placement_minimum_is_interior for r in small)
        and np.all(np.isfinite([r.scaled_a_correction for r in small]))
    )
    correction_text = (
        f"scaled correction range "
        f"[{min(r.scaled_a_correction for r in small):.3g}, "
        f"{max(r.scaled_a_correction for r in small):.3g}]"
        if small else "no small-epsilon data"
    )
    claims.append(
        PublicationClaimRow(
            claim=r"Placement law a = epsilon a1 + O(epsilon^2)",
            status="supported" if placement_ok else "partial",
            evidence=(
                f"a1 from formula (2.7)={a1:.12g}; {correction_text}."
            ),
            limitation=(
                "Bounded scaled corrections on a finite epsilon sample are numerical "
                "consistency with O(epsilon^2), not a proof of the Big-O statement."
            ),
        )
    )

    physical_ok = [
        r for r in small
        if r.expected_bic_resolved and r.physical_bic_verified is True
    ]
    claims.append(
        PublicationClaimRow(
            claim="One non-radiating embedded trapped mode for sampled tuned geometries",
            status=(
                "supported"
                if len(small) >= config.publication_min_asymptotic_points
                and len(physical_ok) == len(small)
                else "partial"
            ),
            evidence=(
                f"{len(physical_ok)}/{len(small)} small-epsilon cases pass "
                "full-matrix singularity, x-even field, open-channel suppression, "
                "off-grid BIE and decay checks."
            ),
            limitation=(
                "Finite sampling supports the theorem numerically; it does not prove "
                "existence for every sufficiently small epsilon."
            ),
        )
    )

    if small:
        C = asymptotic_coefficient(mu, config)
        ratio_errors = np.array(
            [
                abs(r.sigma_over_epsilon_squared - C) / max(abs(C), 1e-30)
                for r in small
            ]
        )
        q = np.array([r.scaled_sigma_remainder for r in small])
        scaling_ok = bool(
            len(small) >= config.publication_min_asymptotic_points
            and np.all(np.isfinite(q))
            and ratio_errors[0] <= ratio_errors[-1]
        )
        claims.append(
            PublicationClaimRow(
                claim=r"sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)",
                status="supported" if scaling_ok else "partial",
                evidence=(
                    f"relative error of sigma/epsilon^2 from C changes from "
                    f"{ratio_errors[0]:.3%} at epsilon={small[0].epsilon:g} to "
                    f"{ratio_errors[-1]:.3%} at epsilon={small[-1].epsilon:g}; "
                    f"scaled-remainder range [{np.min(q):.3g},{np.max(q):.3g}]."
                ),
                limitation=(
                    "Finite asymptotic data support the stated scaling but do not prove "
                    "the remainder bound."
                ),
            )
        )

    scanned = [r for r in small if r.uniqueness_screen_run]
    uniqueness_ok = bool(
        scanned
        and all(r.physical_bic_verified is True for r in scanned)
        and all(r.uniqueness_screen_passed is True for r in scanned)
    )
    claims.append(
        PublicationClaimRow(
            claim="Uniqueness on the sampled real embedded band",
            status="supported" if uniqueness_ok else "partial",
            evidence=(
                f"{sum(r.uniqueness_screen_passed is True for r in scanned)}/"
                f"{len(scanned)} scanned small-epsilon cases found no additional "
                "resolved non-radiating real-axis candidate."
            ),
            limitation=(
                "A finite real-axis SVD screen is numerical evidence, not a theorem-wide "
                "proof and not a validated interval-arithmetic exclusion."
            ),
        )
    )

    claims.append(
        PublicationClaimRow(
            claim="Analyticity in epsilon and epsilon log epsilon",
            status="theoretical-only",
            evidence="This is established analytically in the paper, not by the finite computation.",
            limitation="No finite numerical experiment can prove analyticity.",
        )
    )

    return claims


def write_publication_report(
    output: Path,
    shape: ShapeDiagnostics,
    mfs_rows: list[MFSRow],
    circle: MFSRow,
    summaries: list[ValidationResult],
    claims: list[PublicationClaimRow],
    config: Config,
) -> None:
    final_mfs = mfs_rows[-1]
    lines = [
        "# Numerical validation of Theorem 2.3(iv)",
        "",
        "## Configuration",
        f"- b = {config.b}",
        f"- beta = {config.shape_beta}",
        f"- fixed BEM order M = {config.bem_M}",
        f"- epsilon values = {config.epsilon_values}",
        "",
        "## Independent reference-obstacle quantities",
        f"- area S = {shape.area:.12g}",
        f"- mu = {final_mfs.mu:.12g}",
        f"- nu = {final_mfs.nu:.3e}",
        f"- a1 from theorem formula (2.7) = {final_mfs.a1_formula_2_7:.12g}",
        f"- a1 from equivalent formula (5.47) = {final_mfs.a1_formula_5_47:.12g}",
        f"- relative (2.7)/(5.47) mismatch = {final_mfs.a1_internal_relative_difference:.3e}",
        f"- Example 5.4 leading value -beta/12 = {final_mfs.a1_paper_example_leading:.12g}",
        f"- circle calibration: mu={circle.mu:.12g}, a1={circle.a1_formula_2_7:.3e}",
        "",
        "## Methodological note",
        (
            "Unlike Theorem 2.3(iii), no parity reduction removes the continuous "
            "spectrum here. The reported BIC is therefore obtained from the full "
            "BEM matrix, with the O(epsilon^2) placement correction and sigma tuned "
            "numerically at fixed M."
        ),
        "",
        "## Claim audit",
    ]
    for claim in claims:
        lines.extend(
            [
                f"### {claim.claim}",
                f"- status: **{claim.status}**",
                f"- evidence: {claim.evidence}",
                f"- limitation: {claim.limitation}",
                "",
            ]
        )

    lines.extend(
        [
            "## Per-epsilon summary",
            "",
            "| eps | a_num | (a-eps a1)/eps^2 | kb | sigma | rel sv | open/BIC physical |",
            "|---:|---:|---:|---:|---:|---:|:---:|",
        ]
    )
    for r in summaries:
        lines.append(
            f"| {r.epsilon:.3f} | {r.a_numerical:.8g} | "
            f"{r.scaled_a_correction:.6g} | {r.kb_numerical:.10f} | "
            f"{r.sigma_numerical:.6g} | {r.relative_singular_value_final:.3e} | "
            f"{'PASS' if r.physical_bic_verified is True else 'FAIL/NA'} |"
        )

    (output / "publication_validation_report.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def print_result(
    result: ValidationResult,
    physical: PhysicalDiagnostics | None,
) -> None:
    print("\n--- Theorem 2.3(iv) validation result ---")
    print(f"epsilon = {result.epsilon:.6f}")
    print(f"a1 formula (2.7) = {result.a1_formula_2_7:.12e}")
    print(f"paper Example 5.4 -beta/12 = {result.a1_paper_example_leading:.12e}")
    print(f"leading placement eps*a1 = {result.a_leading:.12e}")
    print(f"numerical tuned placement = {result.a_numerical:.12e}")
    print(
        "(a-eps*a1)/eps^2 = "
        f"{result.scaled_a_correction:.8e}"
    )
    print(
        "placement search interior = "
        f"{'YES' if result.placement_minimum_is_interior else 'no'}"
    )
    print(f"kb asymptotic = {result.kb_asymptotic:.12f}")
    print(f"kb BEM = {result.kb_numerical:.12f}")
    print(f"sigma asym = {result.sigma_asymptotic:.8e}")
    print(f"sigma BEM = {result.sigma_numerical:.8e}")
    print(f"sigma/epsilon^2 = {result.sigma_over_epsilon_squared:.8e}")
    print(f"C=pi^3 mu/b^3 = {result.asymptotic_coefficient:.8e}")
    print(f"relative sigma error = {result.relative_error_sigma:.3%}")
    print(
        "relative singular value at leading (a,sigma) prediction = "
        f"{result.leading_point_relative_singular_value:.3e}"
    )
    print(
        "relative singular value after tuning = "
        f"{result.relative_singular_value_final:.3e}"
    )
    print(f"spectral drop = {result.final_drop_factor:.3e}")
    print(
        f"complex root: success={result.root_success}, "
        f"residual={result.root_residual_norm:.3e}, "
        f"|lambda0|={result.root_eigenvalue_abs:.3e}, "
        f"seed={result.root_seed_source}"
    )
    print(
        "embedded mode resolved = "
        f"{'PASS' if result.expected_bic_resolved else 'FAIL'}"
    )

    if physical is not None:
        print(f"boundary x-even residual = {physical.boundary_x_even_residual:.3e}")
        print(f"field x-even residual = {physical.field_x_even_residual:.3e}")
        print(
            "first open-channel amplitude = "
            f"{physical.first_open_channel_relative_amplitude:.3e}"
        )
        print(
            "signed open-channel coefficients L/R = "
            f"({physical.open_channel_left_real:+.3e}{physical.open_channel_left_imag:+.3e}j) / "
            f"({physical.open_channel_right_real:+.3e}{physical.open_channel_right_imag:+.3e}j); "
            f"phase consistency={physical.open_channel_phase_consistency:.3e}"
        )
        print(f"wall residual = {physical.wall_relative_residual:.3e}")
        print(f"off-grid BIE residual = {physical.boundary_integral_relative_residual:.3e}")
        print(
            f"decay left/right = "
            f"{physical.decay_rate_left:.6e}/{physical.decay_rate_right:.6e}"
        )
        physical_ok = bool(
            physical.boundary_x_even_verified
            and physical.field_x_even_verified
            and physical.open_channel_suppressed
            and physical.walls_verified
            and physical.boundary_integral_verified
            and physical.decay_verified
        )
        print(
            "physical embedded-mode checks = "
            f"{'PASS' if physical_ok else 'FAIL'}"
        )

    if result.uniqueness_screen_run:
        print(
            "whole-band uniqueness screen = "
            f"{'PASS' if result.uniqueness_screen_passed else 'FAIL'} "
            f"(additional resolved BICs={result.additional_resolved_bics})"
        )
    else:
        print("whole-band uniqueness screen = not run for this epsilon")


def main() -> None:
    config = CONFIG
    output = Path(config.output_directory)
    output.mkdir(parents=True, exist_ok=True)

    if abs(config.shape_beta) >= 1.0:
        raise RuntimeError(
            "Use |shape_beta|<1 in this script so the supplied reference contour "
            "remains safely regular/star-shaped for the MFS source construction."
        )
    if config.bem_M < 8:
        raise RuntimeError("bem_M is too small.")

    print(
        "=== Theorem 2.3(iv) numerical validation: "
        "full BEM + complex-root tuning + fixed-M SVD ==="
    )
    print(f"b = {config.b}")
    print(f"shape beta = {config.shape_beta}")
    print(f"fixed BEM order M = {config.bem_M}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"Lambda_2 = {lambda_2(config):.12f}")
    print(f"sqrt(Lambda_1)b = {kb_1(config):.12f}")
    print(f"sqrt(Lambda_2)b = {kb_2(config):.12f}")

    shape = compute_shape_diagnostics(config)
    write_dataclass_csv(output / "shape_diagnostics.csv", [shape])
    plot_shape(config, output)

    print("\n=== SHAPE / Y-AXIS SYMMETRY CHECK ===")
    print(f"reference area = {shape.area:.12g}")
    print(f"min |r'(t)| = {shape.min_reference_speed:.6e}")
    print(f"X-odd defect = {shape.x_odd_defect:.3e}")
    print(f"Y-even defect = {shape.y_even_defect:.3e}")

    if (
        shape.x_odd_defect > config.geometry_symmetry_tolerance
        or shape.y_even_defect > config.geometry_symmetry_tolerance
    ):
        raise RuntimeError("Reference obstacle failed y-axis symmetry check.")

    print("\n=== INDEPENDENT MFS: mu, nu, Psi, a1 ===")
    mu, nu, a1, mfs_rows, circle = run_mfs_convergence(config)
    write_dataclass_csv(output / "mfs_convergence.csv", mfs_rows)
    write_dataclass_csv(output / "mfs_circle_calibration.csv", [circle])
    plot_mfs_convergence(mfs_rows, output, config)

    print(
        f"circle calibration: mu={circle.mu:.12g}, "
        f"nu={circle.nu:.3e}, a1={circle.a1_formula_2_7:.3e}, "
        f"a1(5.47)={circle.a1_formula_5_47:.3e}"
    )
    for row in mfs_rows:
        print(
            f"  N={row.order:>4d}: mu={row.mu:.12g}, nu={row.nu:.3e}, "
            f"a1(2.7)={row.a1_formula_2_7:.12g}, "
            f"a1(5.47)={row.a1_formula_5_47:.12g}, "
            f"mismatch={row.a1_internal_relative_difference:.3e}, "
            f"boundary_res={row.boundary_relative_residual:.3e}"
        )

    print(f"selected mu = {mu:.12g}")
    print(f"selected nu = {nu:.3e} (must vanish by y-axis symmetry)")
    print(f"selected a1 from formula (2.7) = {a1:.12g}")
    print(
        f"selected a1 from formula (5.47) = {mfs_rows[-1].a1_formula_5_47:.12g}; "
        f"relative mismatch={mfs_rows[-1].a1_internal_relative_difference:.3e}"
    )
    print(
        "paper Example 5.4 leading diagnostic -beta/12 = "
        f"{-config.shape_beta/12.0:.12g}"
    )
    print(f"C = pi^3 mu/b^3 = {asymptotic_coefficient(mu, config):.12g}")

    if config.run_small_beta_a1_diagnostic:
        beta_rows: list[MFSRow] = []
        print("\n=== SMALL-BETA EXAMPLE 5.4 DIAGNOSTIC ===")
        for beta_value in config.small_beta_a1_values:
            row, *_ = _solve_reference_mfs(
                replace(config, shape_beta=float(beta_value)),
                config.small_beta_a1_order,
            )
            beta_rows.append(row)
            ratio = row.a1_formula_2_7 / beta_value
            print(
                f"  beta={beta_value:.4f}: a1(2.7)={row.a1_formula_2_7:.12g}, "
                f"a1/beta={ratio:.8g}, -1/12={-1/12:.8g}"
            )
        write_dataclass_csv(output / "a1_small_beta_diagnostic.csv", beta_rows)

    if abs(nu) > config.nu_tolerance:
        raise RuntimeError(
            f"Computed nu={nu:.3e} exceeds configured symmetry tolerance."
        )

    # Full complex-k smoke test inside the first open band.
    eps0 = config.epsilon_values[0]
    a0 = eps0 * a1
    kb_mid = 0.5 * (kb_1(config) + kb_2(config))
    ztest = complex(kb_mid, 1.0e-3)
    Atest = assemble_full_matrix(ztest, eps0, a0, config.bem_M, config)
    if not np.all(np.isfinite(Atest)):
        raise RuntimeError("Complex-k full BEM smoke test produced non-finite values.")
    print("\ncomplex-k full-matrix assembly: PASS")

    summaries: list[ValidationResult] = []
    modes: list[TunedMode] = []
    physical_rows: list[PhysicalDiagnostics] = []
    additional_rows: list[AdditionalCandidate] = []

    print("\n=== EMBEDDED-MODE / PLACEMENT SWEEP ===")
    print(f"epsilons = {config.epsilon_values}")
    previous_mode: TunedMode | None = None

    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.8f} ===")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, mu, config)
        print(f"predicted a leading = {epsilon*a1:.12e}")
        print(f"predicted kb = {kb_asym:.12f}")
        print(f"predicted sigma = {sigma_asym:.8e}")

        (
            result,
            mode,
            physical,
            scan,
            scan_values,
            additional,
        ) = validate_epsilon(
            epsilon,
            shape,
            mu,
            nu,
            a1,
            config,
            previous_mode=previous_mode,
        )

        summaries.append(result)
        modes.append(mode)
        if mode.resolved:
            previous_mode = mode
        if physical is not None:
            physical_rows.append(physical)
        additional_rows.extend(additional)

        print_result(result, physical)
        plot_uniqueness_scan(
            epsilon, scan, scan_values, result, output, config
        )

    write_dataclass_csv(output / "summary.csv", summaries)
    write_dataclass_csv(output / "tuned_modes.csv", modes)
    write_dataclass_csv(output / "physical_diagnostics.csv", physical_rows)
    write_dataclass_csv(output / "additional_candidates.csv", additional_rows)

    plot_summary(summaries, output, config)

    convergence = run_internal_convergence_study(
        summaries, a1, config
    )
    write_dataclass_csv(output / "internal_convergence.csv", convergence)
    plot_internal_convergence(convergence, output, config)

    small = [
        r for r in summaries
        if r.epsilon <= config.publication_asymptotic_epsilon_max + 1e-15
    ]
    write_dataclass_csv(
        output / "publication_asymptotic_window.csv",
        small,
    )

    claims = build_publication_claims(
        shape,
        mu,
        nu,
        a1,
        summaries,
        config,
    )
    write_dataclass_csv(output / "publication_claims.csv", claims)
    write_publication_report(
        output,
        shape,
        mfs_rows,
        circle,
        summaries,
        claims,
        config,
    )

    print("\n=== FINAL SUMMARY ===")
    resolved = sum(r.expected_bic_resolved for r in summaries)
    physical_ok = sum(r.physical_bic_verified is True for r in summaries)
    scanned = [r for r in summaries if r.uniqueness_screen_run]
    unique_scans = sum(r.uniqueness_screen_passed is True for r in scanned)

    print(f"resolved full-matrix embedded mode: {resolved}/{len(summaries)}")
    print(f"physical BIC diagnostics: {physical_ok}/{len(summaries)}")
    print(
        f"whole-band uniqueness scans passed: "
        f"{unique_scans}/{len(scanned)}"
    )
    print("Publication claim audit:")
    for claim in claims:
        print(f"  - {claim.status:>16}: {claim.claim}")
    print(f"\nFiles written to: {output.resolve()}")


if __name__ == "__main__":
    main()
