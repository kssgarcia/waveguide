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
         signed n=1 modal radiation coupling = 0,
     tracking the simple eigenvalue branch from each seed. The full COMPLEX
     matrix residual and both outgoing coefficients are certified afterwards;
     continuation in epsilon is used only from a physically certified mode;
  5. uses fixed-M SVD only as an independent local certification of that root;
  6. reconstructs the FULL field and verifies:
         - wall Dirichlet values,
         - off-grid boundary integral equation,
         - evenness under x -> -x (expected by the proof),
         - suppression of the first open transverse channel,
         - exponential decay with rate sigma;
  7. optionally screens a refined finite real grid strictly inside the band at
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
only as a separate paper-example diagnostic. Substitution of the printed
f0/f1 into (5.47) gives -beta/2; the last coefficient in that example conflicts
with its preceding formulas. No printed coefficient is imposed on the solver.

This is a numerical validation, not a proof of theorem-wide existence,
uniqueness, analyticity, or Big-O bounds.
"""

import csv
import json
from functools import lru_cache
import execution_harness as harness
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
ASSEMBLY_VERSION = "2026-10-13-Y-signed-channel-certified-v3"


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

    # Independent numerical evaluation of the theorem quantity a1.
    # Psi is obtained from the exterior-Laplace MFS and formula (2.7) is
    # integrated directly.  The -beta/12 expression is NEVER used to tune a.
    mfs_orders: tuple[int, ...] = (60, 80, 120, 160)
    mfs_source_scale: float = 0.65
    mfs_rcond: float = 1.0e-12
    a1_quadrature_points: int = 4096
    mfs_circle_mu_tolerance: float = 1.0e-6
    mfs_circle_a1_tolerance: float = 1.0e-8
    mfs_convergence_relative_tolerance: float = 5.0e-6

    # Robustness checks for a1 itself.  These are diagnostics only: the primary
    # a1 is still the converged direct evaluation of (2.7).
    run_a1_robustness_study: bool = True
    a1_robustness_source_scales: tuple[float, ...] = (0.55, 0.65, 0.75)
    a1_robustness_quadrature_points: tuple[int, ...] = (1024, 2048, 4096)
    a1_robustness_order: int = 160
    a1_robustness_relative_tolerance: float = 2.0e-4

    run_small_beta_a1_diagnostic: bool = True
    small_beta_a1_values: tuple[float, ...] = (0.02, 0.05, 0.10)
    small_beta_a1_order: int = 120

    # Independent a1 estimate from the numerically located BIC placement branch:
    # a_BIC/epsilon = a1 + O(epsilon).  This fit is reported separately and is
    # never fed back into the solver.
    placement_fit_epsilon_max: float = 0.10
    placement_fit_min_points: int = 3
    placement_fit_relative_tolerance: float = 0.10

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

    # Stage 2: solve the spectral crossing and signed modal coupling; certify
    # the full complex operator afterwards. Track the eigenvalue of A
    # closest to zero. Variables are c and log(sigma/sigma_asym), where
    # a=epsilon*a1+c*epsilon^2.
    root_max_nfev: int = 120
    root_jacobian_step: float = 1.0e-3
    root_jacobian_condition_limit: float = 1.0e10
    root_branch_overlap_minimum: float = 0.5
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

    # Spectral certification at fixed M.  We require both the usual condition-
    # like relative singular value and a smallest singular value normalized by
    # a typical row norm.  The second criterion prevents a huge smax from making
    # an otherwise non-small smin look acceptable.
    relative_near_singular_tolerance: float = 1.0e-4
    scaled_smin_tolerance: float = 1.0e-4
    singular_vector_residual_tolerance: float = 1.0e-8
    equilibrated_residual_tolerance: float = 1.0e-8
    equilibrated_smin_tolerance: float = 1.0e-8
    field_modal_terms: int = 128
    field_modal_relative_tolerance: float = 1.0e-12
    physical_modal_projection_tolerance: float = 1.0e-7
    mfs_offgrid_residual_tolerance: float = 1.0e-7
    mfs_identity_relative_tolerance: float = 1.0e-7
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
    uniqueness_refinement_levels: tuple[int, ...] = (1, 2)
    uniqueness_candidate_stability_tolerance: float = 1.0e-7

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
    output_directory: str = str(Path(__file__).resolve().parent / "results" / "theorem_2_3_iv" / "runs")
    execution_workers: int = 1
    execution_resume: bool = True


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
    offgrid_relative_residual: float = math.nan
    mu_boundary_identity: float = math.nan
    identity_relative_error: float = math.nan


@dataclass
class A1RobustnessRow:
    order: int
    source_scale: float
    quadrature_points: int
    mu: float
    nu: float
    a1_formula_2_7: float
    relative_difference_from_primary: float
    boundary_relative_residual: float


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
    scaled_singular_value: float
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
    normalized_residual: float = math.nan
    equilibrated_residual: float = math.nan
    equilibrated_sigma_min: float = math.nan
    root_jacobian_condition: float = math.inf
    modal_coupling_residual: float = math.nan
    root_optimizer_message: str = ""
    root_attempts_json: str = "[]"


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
    modal_projection_discrepancy: float = math.nan
    modal_projection_verified: bool = False


@dataclass
class AdditionalCandidate:
    epsilon: float
    a: float
    kb: float
    sigma: float
    relative_singular_value: float
    scaled_singular_value: float
    drop_factor: float
    boundary_even_residual: float
    open_channel_relative_amplitude: float
    nonradiating_resolved_candidate: bool
    status: str = "ambiguous"
    reason: str = ""
    optimizer_success: bool = False
    minimum_is_interior: bool = False
    scan_level: int = 0


@dataclass
class ScreenResult:
    scan: np.ndarray
    values: np.ndarray
    candidates: list[AdditionalCandidate]
    diagnostics: dict

    def __iter__(self):
        yield self.scan
        yield self.values
        yield self.candidates


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
    scaled_singular_value_final: float
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
    mode_exists_supported: bool = False
    uniqueness_supported: bool = False
    uniqueness_status: str = "ambiguous/not-established"
    screen_diagnostics_json: str = "{}"
    normalized_residual: float = math.nan
    equilibrated_residual: float = math.nan
    equilibrated_sigma_min: float = math.nan


@dataclass
class InternalConvergenceRow:
    epsilon: float
    parameter: str
    value: float
    kb: float
    sigma_bem: float
    relative_singular_value: float
    scaled_singular_value: float
    relative_kb_shift_from_reference: float
    relative_sigma_shift_from_reference: float
    a_numerical: float = math.nan
    a_shift_from_reference: float = math.nan
    physical_bic_verified: bool = False
    open_channel_relative_amplitude: float = math.nan
    normalized_residual: float = math.nan
    status: str = "ambiguous"


@dataclass
class PlacementFitDiagnostics:
    selection: str
    point_count: int
    epsilon_max: float
    a1_formula_2_7: float
    a1_fit_linear: float
    a2_fit_linear: float
    linear_rmse: float
    linear_relative_difference: float
    a1_fit_quadratic: float
    quadratic_relative_difference: float


@dataclass
class SmallBetaFitDiagnostics:
    point_count: int
    max_beta: float
    fitted_slope_through_origin: float
    paper_slope_minus_one_twelfth: float
    relative_difference_from_paper: float


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
    source_scale_override: float | None = None,
    quadrature_points_override: int | None = None,
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
    rho = float(cfg.mfs_source_scale if source_scale_override is None else source_scale_override)
    if not 0<rho<1:raise ValueError("MFS sources must lie inside the star-shaped obstacle")
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
    Q = int(cfg.a1_quadrature_points if quadrature_points_override is None else quadrature_points_override)
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

    # Algebraically equivalent cross-check (5.47), using the paper identity
    #   L0 Y = 1/2 (Y - Psi|_Gamma).
    # It uses the same computed Psi, so it is NOT an independent PDE solve.
    # Agreement checks normalization/sign conventions only.  The independent
    # numerical a1 check is the epsilon->0 BIC placement-branch intercept.
    L0Y = 0.5 * (Yq - Psi)
    numerator_547 = dt * np.sum(
        Yq * Xpq * L0Y - 3.0 * Xq * Ypq * L0Y
    )
    denominator_547 = 2.0 * area + PI * mu
    a1_547 = float(numerator_547 / denominator_547)
    a1_internal_relative_difference = abs(a1_27 - a1_547) / max(
        abs(a1_27), abs(a1_547), 1.0e-30
    )

    normal_x=-Ypq/np.hypot(Xpq,Ypq);normal_y=Xpq/np.hypot(Xpq,Ypq)
    r2q=dxq*dxq+dyq*dyq
    offgrid=-(dxq*normal_x[:,None]+dyq*normal_y[:,None])/(2*PI*r2q)@c
    offgrid_res=float(np.linalg.norm(offgrid-normal_y)/np.linalg.norm(normal_y))
    mu_identity=float((area+dt*np.sum(Psi*Xpq))/(2*PI))
    identity_error=abs(mu_identity-mu)/max(abs(mu),1e-30)
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
        offgrid_relative_residual=offgrid_res,mu_boundary_identity=mu_identity,
        identity_relative_error=float(identity_error),
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
            abs(rows[-1].a1_formula_2_7), config.mfs_circle_a1_tolerance
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

    for row in [circle,final]:
        if (not np.all(np.isfinite([row.mu,row.a1_formula_2_7,row.offgrid_relative_residual,row.identity_relative_error]))
            or row.offgrid_relative_residual>config.mfs_offgrid_residual_tolerance
            or row.identity_relative_error>config.mfs_identity_relative_tolerance):
            raise RuntimeError("MFS off-grid/identity certificate failed")
    if len(rows)<2:raise RuntimeError("MFS convergence needs at least two orders")
    if max(mu_change,a1_change)>config.mfs_convergence_relative_tolerance:
        raise RuntimeError("MFS mu/a1 order convergence not established")
    if final.mu <= 0.0:
        raise RuntimeError(f"Theorem requires mu>0; MFS returned {final.mu}.")

    return (
        float(final.mu),
        float(final.nu),
        float(final.a1_formula_2_7),
        rows,
        circle,
    )


def run_a1_robustness_study(
    primary_a1: float,
    config: Config,
) -> list[A1RobustnessRow]:
    """Check formula (2.7) against MFS source placement and quadrature resolution.

    This is deliberately independent of the waveguide/BEM eigenvalue calculation.
    It only probes the exterior Laplace problem that defines Psi, mu and therefore a1.
    """
    if not config.run_a1_robustness_study:
        return []

    rows: list[A1RobustnessRow] = []
    for rho in config.a1_robustness_source_scales:
        for q in config.a1_robustness_quadrature_points:
            row, *_ = _solve_reference_mfs(
                config,
                config.a1_robustness_order,
                source_scale_override=float(rho),
                quadrature_points_override=int(q),
            )
            rel = abs(row.a1_formula_2_7 - primary_a1) / max(abs(primary_a1), 1.0e-30)
            rows.append(
                A1RobustnessRow(
                    order=int(config.a1_robustness_order),
                    source_scale=float(rho),
                    quadrature_points=int(q),
                    mu=float(row.mu),
                    nu=float(row.nu),
                    a1_formula_2_7=float(row.a1_formula_2_7),
                    relative_difference_from_primary=float(rel),
                    boundary_relative_residual=float(row.boundary_relative_residual),
                )
            )

    worst = max((r.relative_difference_from_primary for r in rows), default=0.0)
    if worst > config.a1_robustness_relative_tolerance:
        print(
            "WARNING: a1 robustness spread exceeds tolerance: "
            f"max relative difference={worst:.3e}."
        )
    if (not rows or not np.all(np.isfinite([r.relative_difference_from_primary for r in rows]))
        or worst>config.a1_robustness_relative_tolerance):
        raise RuntimeError("a1 source-scale/quadrature robustness not established")
    return rows


def fit_small_beta_a1(rows: list[MFSRow]) -> SmallBetaFitDiagnostics | None:
    if not rows:
        return None
    beta = np.asarray([r.beta for r in rows], dtype=float)
    values = np.asarray([r.a1_formula_2_7 for r in rows], dtype=float)
    denom = float(np.dot(beta, beta))
    if denom <= 0.0:
        return None
    slope = float(np.dot(beta, values) / denom)
    paper_slope = -1.0 / 12.0
    rel = abs(slope - paper_slope) / abs(paper_slope)
    return SmallBetaFitDiagnostics(
        point_count=len(rows),
        max_beta=float(np.max(beta)),
        fitted_slope_through_origin=slope,
        paper_slope_minus_one_twelfth=paper_slope,
        relative_difference_from_paper=float(rel),
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

    # The periodic Green representation is evaluated after wrapping each vertical
    # image separation into the nearest representative of the exact 4b period.
    # Therefore |Delta y_wrapped| <= 2b.  Keep a conservative radius check only
    # for the remaining horizontal extent.
    max_dx = 2.0 * float(np.max(np.abs(X)))
    max_wrapped_dy = 2.0 * float(config.b)
    max_r = math.hypot(max_dx, max_wrapped_dy)
    green_limit = 0.99 * (4.0 * config.b)
    if max_r > green_limit:
        return False, (
            f"Wrapped Green-series radius bound exceeded: "
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
    """Dirichlet strip Green function with exact periodic reduction in y.

    The lattice cylindrical expansion is local in radius even though the periodic
    Green function itself is exactly periodic.  Evaluating an unwrapped image
    separation near one full period can therefore magnify harmonic truncation
    error.  We always reduce both vertical separations modulo the exact 4b period
    before calling greens_periodic/greens_periodic_reg.
    """
    b = float(config.b)
    d = 2.0 * b
    period = 2.0 * d
    k = complex(kb / b)
    coefficients = lattice.lattice_sums(
        period,
        k,
        beta=0.0,
        M=config.lattice_terms,
        Lh=config.harmonic_order,
    )

    def wrap_periodic_y(value: float) -> float:
        return float((value + 0.5 * period) % period - 0.5 * period)

    def green(x: float, y: float, xi: float, eta: float) -> complex:
        yf = float(y) + b
        ys = float(eta) + b
        X = float(x) - float(xi)
        Y1 = wrap_periodic_y(float(y) - float(eta))
        Y2 = wrap_periodic_y(float(y) + float(eta) + 2*b)
        term1 = lattice.greens_periodic(X, Y1, coefficients, k, period)
        term2 = lattice.greens_periodic(X, Y2, coefficients, k, period)
        return term1 - term2

    def green_regularized(x: float, y: float, xi: float, eta: float) -> complex:
        yf = float(y) + b
        ys = float(eta) + b
        X = float(x) - float(xi)
        Y1 = wrap_periodic_y(float(y) - float(eta))
        Y2 = wrap_periodic_y(float(y) + float(eta) + 2*b)
        term1 = lattice.greens_periodic_reg(X, Y1, coefficients, k, period)
        term2 = lattice.greens_periodic(X, Y2, coefficients, k, period)
        return term1 - term2

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


@lru_cache(maxsize=4096)
def _kernel_geometry(t: float, epsilon: float, a: float, config: Config) -> tuple[float, ...]:
    return tuple(float(v) for v in obstacle_geometry(t, epsilon, a, config))


def weighted_normal_kernel(
    psi: float,
    theta: float,
    epsilon: float,
    a: float,
    config: Config,
    G: Callable[..., complex],
    G_regularized: Callable[..., complex],
) -> complex:
    x, y, *_ = _kernel_geometry(psi, epsilon, a, config)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = _kernel_geometry(
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


def assemble_full_matrix_scalar(
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


def assemble_full_matrix(
    kb: complex,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> np.ndarray:
    period = 4*config.b
    k = complex(kb)/config.b
    coeffs = lattice.lattice_sums(period,k,beta=0.,M=config.lattice_terms,Lh=config.harmonic_order)
    X,Y,Dx,Dy,Iy,xp,yp,curvature,plans = _assembly_geometry(epsilon,a,M,config)
    shifted = (_periodic_array(X,Y,coeffs,k,False,plans[0])
               - _periodic_array(X,Iy,coeffs,k,False,plans[1]))
    h=config.finite_difference_step
    gx=(shifted[0]-shifted[1])/(2*h);gy=(shifted[2]-shifted[3])/(2*h)
    weighted=xp[None,:]*gy-yp[None,:]*gx
    reg=(_periodic_array(Dx,Dy,coeffs,k,True,plans[2])
         - _periodic_array(Dx,plans[4],coeffs,k,False,plans[3]))
    gx=(reg[0]-reg[1])/(2*h);gy=(reg[2]-reg[3])/(2*h)
    weighted[np.diag_indices(M)]=curvature+xp*gy-yp*gx
    return np.eye(M,dtype=complex)-(4*PI/M)*weighted


def _periodic_array(X,Y,coeffs,k,regularized,plan):
    radius,angular=plan
    argument=k*radius
    result=lattice.J0(argument) if regularized else lattice.H0(argument)
    for order,coefficient in enumerate(coeffs):
        result=result+(1 if order==0 else 2)*coefficient*lattice.jv(order,argument)*angular[order]
    return -.25j*result


@lru_cache(maxsize=16)
def _assembly_geometry(epsilon,a,M,config):
    t=boundary_nodes(M);x,y,xp,yp,xpp,ypp=obstacle_geometry(t,epsilon,a,config)
    h=config.finite_difference_step;period=4*config.b
    xi=np.stack([x+h,x-h,x,x]);eta=np.stack([y,y,y+h,y-h])
    X=x[None,:,None]-xi[:,None,:];Y=y[None,:,None]-eta[:,None,:]
    Iy=y[None,:,None]+eta[:,None,:]+2*config.b
    Dx=x[None,:]-xi;Dy=y[None,:]-eta;Di=y[None,:]+eta+2*config.b
    def prepare(xx,yy):
        yy=(yy+period/2)%period-period/2
        radius=np.hypot(xx,yy)
        if np.any(radius>.99*period):raise ValueError("Outside lattice Green disk")
        angle=PI/2-np.arctan2(yy,xx)
        return radius,np.stack([np.cos(n*angle) for n in range(2*config.harmonic_order+1)])
    plans=(prepare(X,Y),prepare(X,Iy),prepare(Dx,Dy),prepare(Dx,Di),Di)
    curvature=(xpp*yp-ypp*xp)/(4*PI*(xp*xp+yp*yp))
    return X,Y,Dx,Dy,Iy,xp,yp,curvature,plans



def _matrix_reference_scale(A: np.ndarray) -> float:
    """Typical row scale used to prevent one huge singular value from hiding smin."""
    row_norms = np.linalg.norm(A, axis=1)
    finite = row_norms[np.isfinite(row_norms) & (row_norms > 0.0)]
    if finite.size!=len(row_norms) or not np.all(np.isfinite(A)):
        raise ValueError("Nonfinite or zero row scale")
    return max(float(np.median(finite)), np.finfo(float).tiny)


def full_singular_metrics(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, float]:
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    s = np.linalg.svd(A, compute_uv=False)
    smax = float(s[0])
    smin = float(s[-1])
    rel = smin / max(smax, np.finfo(float).tiny)
    scaled = smin / _matrix_reference_scale(A)
    return smin, smax, float(rel), float(scaled)


def full_singular_pair(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, float, np.ndarray]:
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    _, s, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= max(float(np.linalg.norm(vector)), 1.0e-30)

    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))

    smax = float(s[0])
    smin = float(s[-1])
    rel = smin / max(smax, np.finfo(float).tiny)
    scaled = smin / _matrix_reference_scale(A)
    return smin, smax, float(rel), float(scaled), vector


def relative_sv_from_sigma(
    sigma: float,
    epsilon: float,
    a: float,
    config: Config,
) -> float:
    kb = kb_from_sigma(float(sigma), config)
    if not (kb_1(config) < kb < kb_2(config)):
        return math.nan
    try:
        _, _, rel, scaled = full_singular_metrics(
            kb, epsilon, a, config.bem_M, config
        )
    except Exception:
        return math.nan
    return float(rel)


def singularity_measure_from_sigma(
    sigma: float,
    epsilon: float,
    a: float,
    config: Config,
) -> float:
    """Combined scale-robust singularity measure used for local discovery."""
    kb = kb_from_sigma(float(sigma), config)
    if not (kb_1(config) < kb < kb_2(config)):
        return math.nan
    try:
        _, _, rel, scaled = full_singular_metrics(
            kb, epsilon, a, config.bem_M, config
        )
    except Exception:
        return math.nan
    return float(max(rel, scaled))


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


def modal_open_channel_coefficients(kb, epsilon, a, vector, config):
    """Exact n=1 outgoing coefficient of the represented double-layer density.

    G_n=-exp(-q|x-xi|) phi_n(y)phi_n(eta)/(2b*q), q=-i*p1.
    Both sides use exp(+i*p1*distance); the source derivative gives the signs
    below. No y-quadrature, harmonic expansion or finite differences is used.
    """
    t=boundary_nodes(len(vector));xi,eta,xp,yp,*_=obstacle_geometry(t,epsilon,a,config)
    q1=math.sqrt((kb/config.b)**2-lambda_1(config));p=PI/(2*config.b)
    phi=np.cos(p*eta);dphi=-p*np.sin(p*eta);dt=2*PI/len(t)
    right=dt*np.sum(vector*np.exp(-1j*q1*xi)*(xp*dphi+1j*q1*yp*phi))
    left=dt*np.sum(vector*np.exp(+1j*q1*xi)*(xp*dphi-1j*q1*yp*phi))
    factor=-1j/(2*config.b*q1)
    return complex(factor*left),complex(factor*right)


def certify_full_mode(kb,epsilon,a,M,interior,drop,config):
    A=assemble_full_matrix(complex(kb),epsilon,a,M,config)
    if not np.all(np.isfinite(A)):raise ValueError("Nonfinite mode matrix")
    _,sv,vh=np.linalg.svd(A);v=vh.conj().T[:,-1];rows=np.linalg.norm(A,axis=1)
    if sv[0]<=0 or np.any(rows<=np.finfo(float).tiny):raise ValueError("Degenerate mode matrix")
    B=A/rows[:,None];bs=np.linalg.svd(B,compute_uv=False)
    residual=float(np.linalg.norm(A@v)/(sv[0]*np.linalg.norm(v)))
    eq=float(np.linalg.norm(B@v)/(bs[0]*np.linalg.norm(v)))
    supported=bool(np.all(np.isfinite([residual,eq,bs[-1],drop])) and interior
                   and drop>=config.minimum_drop_factor
                   and sv[-1]/sv[0]<=config.relative_near_singular_tolerance
                   and sv[-1]/_matrix_reference_scale(A)<=config.scaled_smin_tolerance
                   and residual<=config.singular_vector_residual_tolerance
                   and eq<=config.equilibrated_residual_tolerance
                   and bs[-1]<=config.equilibrated_smin_tolerance)
    return dict(supported=supported,normalized_residual=residual,
                equilibrated_residual=eq,equilibrated_sigma_min=float(bs[-1]))


def near_zero_eigenpair(
    kb: float,
    epsilon: float,
    a: float,
    M: int,
    config: Config,
    anchor: np.ndarray | None = None,
) -> tuple[complex, np.ndarray, float]:
    """Return the eigenvalue of A closest to zero and its right eigenvector.

    A true BIC at real (k,a) makes the full outgoing BEM operator singular.
    Each root solve tracks its seeded eigenvector overlap. The root uses a
    signed radiation constraint; the full complex eigenvalue is still tested
    independently before spectral certification.
    """
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, a, M, config)
    eigvals, eigvecs = np.linalg.eig(A)
    if anchor is None:
        idx = int(np.argmin(np.abs(eigvals)))
    else:
        overlaps=np.abs(anchor.conj()@eigvecs)/np.linalg.norm(eigvecs,axis=0)
        idx=int(np.argmax(overlaps))
        if overlaps[idx]<config.root_branch_overlap_minimum:
            raise ValueError("Lost tracked eigenvalue branch")
    lam = complex(eigvals[idx])
    vec = np.asarray(eigvecs[:, idx], dtype=np.complex128)
    vec /= max(float(np.linalg.norm(vec)), 1.0e-30)
    if anchor is not None:
        vec *= np.exp(-1j*np.angle(np.vdot(anchor,vec)))
    else:
        pivot = int(np.argmax(np.abs(vec)))
        if abs(vec[pivot]) > 0.0:
            vec *= np.exp(-1j * np.angle(vec[pivot]))
    scale = _matrix_reference_scale(A)
    return lam, vec, scale


def _root_residual(x,epsilon,a1,sigma_asym,config,anchor=None):
    """Spectral crossing plus signed radiation coupling in the x-even branch.

    For a real x-even density the modal source integral is real (left=right),
    and its zero cancels the outgoing n=1 wave. This avoids solving its squared
    leakage through Im(lambda). The final COMPLEX matrix residual and both
    COMPLEX outgoing coefficients remain independent acceptance requirements.
    """
    a,sigma,kb=_joint_candidate_from_variables(epsilon,a1,sigma_asym,float(x[0]),float(x[1]),config)
    ok,reason=geometry_admissibility(epsilon,a,config)
    if not ok:raise ValueError(reason)
    lam,vec,scale=near_zero_eigenpair(kb,epsilon,a,config.bem_M,config,anchor)
    q1=math.sqrt((kb/config.b)**2-lambda_1(config))
    left,right=modal_open_channel_coefficients(kb,epsilon,a,vec.real.astype(complex),config)
    # Undo -i/(2b*q1), then normalize the O(epsilon^2) source coupling.
    coupling=(1j*2*config.b*q1*(left+right)/2).real/max(epsilon**2,1e-30)
    residual=np.array([lam.real/scale,coupling],float)
    if not np.all(np.isfinite(residual)):raise ValueError("Nonfinite root residual")
    return residual


def _root_jacobian(x,epsilon,a1,sigma_asym,config,anchor=None):
    h=config.root_jacobian_step
    columns=[]
    for j in range(2):
        step=np.eye(2)[j]*h
        plus=x+step;minus=x-step
        if j==1:
            low,high=_log_sigma_bounds(sigma_asym,config)
            plus[1]=min(plus[1],high);minus[1]=max(minus[1],low)
        columns.append((_root_residual(plus,epsilon,a1,sigma_asym,config,anchor)
                        -_root_residual(minus,epsilon,a1,sigma_asym,config,anchor))/(plus[j]-minus[j]))
    return np.column_stack(columns)


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


def _log_sigma_bounds(sigma_asym,config):
    low,high=sigma_band(config);limit=math.log(config.sigma_search_factor)
    bounds=(max(-limit,math.log(low/sigma_asym)+1e-8),
            min(limit,math.log(high/sigma_asym)-1e-8))
    if not bounds[0]<bounds[1]:raise ValueError("Empty physical sigma search window")
    return bounds


def _clip_seed(
    seed: np.ndarray,
    c_bound: float,
    log_bounds: tuple[float,float],
) -> np.ndarray:
    margin = 1.0e-9
    return np.array([
        float(np.clip(seed[0], -c_bound + margin, c_bound - margin)),
        float(np.clip(seed[1], log_bounds[0] + margin, log_bounds[1] - margin)),
    ])


def _solve_complex_eigenvalue_root(epsilon,a1,sigma_asym,c_bound,seeds,config):
    """Track each seeded branch and retain complete failure/endpoint diagnostics."""
    log_bounds=_log_sigma_bounds(sigma_asym,config);attempts=[];results=[]
    for source,raw_seed in seeds:
        if source.startswith('multistart') and any(r[4] for r in results):break
        seed=_clip_seed(np.asarray(raw_seed,float),c_bound,log_bounds)
        try:
            a,_,kb=_joint_candidate_from_variables(epsilon,a1,sigma_asym,*seed,config)
            _,anchor,_=near_zero_eigenpair(kb,epsilon,a,config.bem_M,config)
            result=least_squares(lambda x:_root_residual(x,epsilon,a1,sigma_asym,config,anchor),
                seed,jac=lambda x:_root_jacobian(x,epsilon,a1,sigma_asym,config,anchor),
                bounds=([-c_bound,log_bounds[0]],[c_bound,log_bounds[1]]),
                xtol=config.root_xtol,ftol=config.root_ftol,gtol=config.root_gtol,
                max_nfev=config.root_max_nfev,x_scale='jac')
            a,_,kb=_joint_candidate_from_variables(epsilon,a1,sigma_asym,*result.x,config)
            lam,vec,scale=near_zero_eigenpair(kb,epsilon,a,config.bem_M,config,anchor)
            norm=abs(lam)/scale
            jac=_root_jacobian(result.x,epsilon,a1,sigma_asym,config,anchor)
            condition=float(np.linalg.cond(jac));coupling=float(abs(_root_residual(result.x,epsilon,a1,sigma_asym,config,anchor)[1]))
            result.jacobian_condition=condition;result.coupling_residual=coupling
            accepted=bool(result.success and np.isfinite(norm) and norm<=config.root_residual_tolerance
                          and condition<=config.root_jacobian_condition_limit
                          and coupling<=config.root_residual_tolerance)
            attempts.append(dict(source=source,optimizer_success=bool(result.success),message=str(result.message),
                root_norm=float(norm),coupling_residual=coupling,jacobian_condition=condition,
                correction=float(result.x[0]),log_sigma_ratio=float(result.x[1]),accepted=accepted))
            results.append((result,source,lam,float(norm),accepted))
        except Exception as exc:
            attempts.append(dict(source=source,optimizer_success=False,accepted=False,error=f'{type(exc).__name__}: {exc}'))
        # Evaluate ALL priority seeds. Only skip fallback once one is certified
        # by the spectral and signed-coupling endpoints, never by lambda alone.
        if source.startswith('multistart') and any(r[4] for r in results):
            break
    if not results:raise RuntimeError('All root attempts failed: '+json.dumps(attempts))
    candidates=[r for r in results if r[4]] or results
    chosen=min(candidates,key=lambda r:(not r[4],r[3]+getattr(r[0],'coupling_residual',math.inf)))
    chosen[0].attempts=attempts
    return chosen[:4]


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
    if not smin < sigma < smax:
        raise ValueError("Spectral variable outside the open band; no clipping")
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
    _, _, rel, scaled = full_singular_metrics(
        kb, epsilon, a, config.bem_M, config
    )
    measure = max(rel, scaled)
    if not math.isfinite(measure):
        raise ValueError("Invalid coarse spectral objective")
    return float(math.log10(max(measure, np.finfo(float).tiny)))


def _polish_sigma(
    epsilon: float,
    a: float,
    sigma_seed: float,
    config: Config,
) -> tuple[float, float, float, float, float]:
    band_lo, band_hi = sigma_band(config)
    lo = max(band_lo, sigma_seed / config.polish_sigma_factor)
    hi = min(band_hi, sigma_seed * config.polish_sigma_factor)
    if not lo < hi:
        lo, hi = band_lo, band_hi

    result: Any = minimize_scalar(
        lambda s: math.log10(
            max(
                singularity_measure_from_sigma(float(s), epsilon, a, config),
                np.finfo(float).tiny,
            )
        ),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    if not result.success:raise RuntimeError("Scalar sigma polishing failed")
    sigma = float(result.x)
    kb = kb_from_sigma(sigma, config)
    smin, smax, rel, scaled = full_singular_metrics(
        kb, epsilon, a, config.bem_M, config
    )
    return sigma, kb, smin, rel, scaled



def _spectral_drop_factor(
    epsilon: float,
    a: float,
    sigma: float,
    measure_center: float,
    config: Config,
) -> float:
    band_lo, band_hi = sigma_band(config)
    delta = config.drop_probe_relative_sigma
    left_sigma = max(band_lo, sigma * (1.0 - delta))
    right_sigma = min(band_hi, sigma * (1.0 + delta))
    measure_left = singularity_measure_from_sigma(left_sigma, epsilon, a, config)
    measure_right = singularity_measure_from_sigma(right_sigma, epsilon, a, config)
    return float(
        min(measure_left, measure_right)
        / max(measure_center, np.finfo(float).tiny)
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

    _, _, leading_rel, leading_scaled = full_singular_metrics(
        kb_asym, epsilon, a_leading, config.bem_M, config
    )

    log_low,log_high = _log_sigma_bounds(sigma_asym,config)
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
            bounds=[(-c_bound, c_bound), (log_low, log_high)],
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

        seeds.append(('theorem-leading', np.array([0.0, 0.0], dtype=float)))
        seeds.append(('coarse-svd', np.asarray(coarse.x, dtype=float)))

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
        sigma_interior = config.placement_boundary_fraction*log_low <= logr <= config.placement_boundary_fraction*log_high

        chosen = (
            result, source, lam, root_norm, c_bound,
            bool(continuation_used and source == 'epsilon-continuation'),
        )

        # Stop widening once a genuine interior complex root has been found.
        if (c_interior and sigma_interior and result.success
                and result.jacobian_condition <= config.root_jacobian_condition_limit
                and result.coupling_residual <= config.root_residual_tolerance
                and root_norm <= config.root_residual_tolerance):
            break

    if chosen is None:
        raise RuntimeError('BIC tuning did not produce a complex-eigenvalue root candidate.')

    result, source, lam, root_norm, c_bound, continuation_seed_used = chosen
    correction = float(result.x[0])
    log_ratio = float(result.x[1])
    a, sigma, kb = _joint_candidate_from_variables(
        epsilon, a1, sigma_asym, correction, log_ratio, config
    )

    smin, smax, rel, scaled = full_singular_metrics(
        kb, epsilon, a, config.bem_M, config
    )
    scaled_correction = (a - a_leading) / max(epsilon**2, 1.0e-30)
    placement_interior = bool(
        abs(scaled_correction) <= config.placement_boundary_fraction * c_bound
        and config.placement_boundary_fraction*log_low <= log_ratio <= config.placement_boundary_fraction*log_high
    )

    measure = max(rel, scaled)
    drop = _spectral_drop_factor(epsilon, a, sigma, measure, config)
    root_ok = bool(result.success and math.isfinite(root_norm)
                   and root_norm <= config.root_residual_tolerance
                   and result.coupling_residual <= config.root_residual_tolerance
                   and result.jacobian_condition <= config.root_jacobian_condition_limit)
    certificate=certify_full_mode(kb,epsilon,a,config.bem_M,placement_interior,drop,config)
    resolved = bool(
        certificate["supported"]
        and placement_interior
        and root_ok
        and rel <= config.relative_near_singular_tolerance
        and scaled <= config.scaled_smin_tolerance
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
        scaled_singular_value=float(scaled),
        drop_factor=float(drop),
        search_correction_bound=float(c_bound),
        search_minimum_is_interior=bool(placement_interior),
        resolved=resolved,
        leading_point_relative_singular_value=float(leading_rel),
        root_success=bool(result.success),
        root_residual_norm=float(root_norm),
        root_eigenvalue_real=float(lam.real),
        root_eigenvalue_imag=float(lam.imag),
        root_eigenvalue_abs=float(abs(lam)),
        root_nfev=int(getattr(result, 'nfev', -1)),
        root_seed_source=str(source),
        continuation_seed_used=bool(continuation_seed_used),
        normalized_residual=certificate['normalized_residual'],
        equilibrated_residual=certificate['equilibrated_residual'],
        equilibrated_sigma_min=certificate['equilibrated_sigma_min'],
        root_jacobian_condition=result.jacobian_condition,
        modal_coupling_residual=result.coupling_residual,
        root_optimizer_message=str(result.message),root_attempts_json=json.dumps(result.attempts),
    )


# ---------------------------------------------------------------------------
# Full-field reconstruction and BIC diagnostics
# ---------------------------------------------------------------------------


def make_extended_field_green(kb: float, config: Config) -> Callable[..., complex]:
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
        if abs(X)>=b:
            n=np.arange(1,config.field_modal_terms+1,dtype=float);p=n*PI/d
            q=np.sqrt(p*p-k*k+0j)
            q=np.where(p<k.real,-1j*np.sqrt(k*k-p*p+0j),q)
            terms=np.exp(-q*abs(X))*np.sin(p*yf)*np.sin(p*ys)/q
            scale=max(float(np.sum(np.abs(terms))),np.finfo(float).tiny)
            if (not np.all(np.isfinite(terms)) or np.sum(np.abs(terms[len(terms)//2:]))>config.field_modal_relative_tolerance*scale):
                raise ValueError("Modal field Green did not converge")
            return complex(-np.sum(terms)/d)
        Y1 = wrap(y - eta)
        Y2 = wrap(y + eta + 2*b)
        return (
            lattice.greens_periodic(X, Y1, coeffs, k, period)
            - lattice.greens_periodic(X, Y2, coeffs, k, period)
        )

    green.kb=kb
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


def _modal_field(points,vector,theta,epsilon,a,kb,config):
    n=np.arange(1,config.field_modal_terms+1,dtype=float);p=n*PI/(2*config.b)
    k=kb/config.b;q=np.sqrt(p*p-k*k+0j)
    q=np.where(p<k,-1j*np.sqrt(k*k-p*p+0j),q)
    xi,eta,xp,yp,*_=obstacle_geometry(theta,epsilon,a,config)
    points=np.asarray(points,float);X=points[:,0,None]-xi[None,:]
    if np.any(abs(X)<config.b):raise ValueError("Modal reconstruction requires separated sections")
    exp=np.exp(-abs(X)[:,:,None]*q[None,None,:])
    source=(xp[:,None]*p[None,:]*np.cos((eta[:,None]+config.b)*p[None,:])/q[None,:]
            -np.sign(X)[:,:,None]*yp[None,:,None]*np.sin((eta[:,None]+config.b)*p[None,:]))
    terms=-(2*PI/len(theta))/(2*config.b)*np.sum(exp*source*vector[None,:,None],axis=1)
    terms*=np.sin((points[:,1,None]+config.b)*p[None,:])
    if not np.all(np.isfinite(terms)):raise ValueError("Nonfinite modal field")
    tail=np.sum(abs(terms[:,len(n)//2:]),axis=1)
    scale=np.maximum(np.sum(abs(terms),axis=1),np.finfo(float).tiny)
    if np.any(tail>config.field_modal_relative_tolerance*scale):raise ValueError("Modal field tail did not converge")
    return np.sum(terms,axis=1)


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

    far=[i for i,(x,y) in enumerate(points) if np.min(abs(x-obstacle_geometry(theta,epsilon,a,config)[0]))>=config.b]
    if far:
        values[far]=_modal_field([points[i] for i in far],boundary_vector,theta,epsilon,a,float(G.kb),config)
    for i, (x, y) in enumerate(points):
        if i in far:continue
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


def _offgrid_kernel_array(x,y,xi,eta,xp,yp,kb,config):
    """Batched version of the unchanged centered source-difference kernel."""
    h=config.finite_difference_step;period=4*config.b;k=complex(kb)/config.b
    coefficients=lattice.lattice_sums(period,k,beta=0.,M=config.lattice_terms,Lh=config.harmonic_order)
    xis=np.stack([xi+h,xi-h,xi,xi]);etas=np.stack([eta,eta,eta+h,eta-h])
    X=x[None,:,None]-xis[:,None,:]
    Y=y[None,:,None]-etas[:,None,:];IY=y[None,:,None]+etas[:,None,:]+2*config.b
    def prepare(Y):
        Y=(Y+period/2)%period-period/2;radius=np.hypot(X,Y)
        if np.any(radius<=np.finfo(float).tiny):raise ValueError("Off-grid kernel source collision")
        if np.any(radius>.99*period):raise ValueError("Outside lattice expansion disk")
        angle=PI/2-np.arctan2(Y,X)
        return radius,np.stack([np.cos(n*angle) for n in range(2*config.harmonic_order+1)])
    values=(_periodic_array(X,Y,coefficients,k,False,prepare(Y))
            -_periodic_array(X,IY,coefficients,k,False,prepare(IY)))
    gx=(values[0]-values[1])/(2*h);gy=(values[2]-values[3])/(2*h)
    kernel=xp[None,:]*gy-yp[None,:]*gx
    if not np.all(np.isfinite(kernel)):raise ValueError("Nonfinite off-grid kernel")
    return kernel


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
    x,y,*_=obstacle_geometry(targets,epsilon,a,config)
    xi,eta,xp,yp,*_=obstacle_geometry(source_theta,epsilon,a,config)
    kernel=_offgrid_kernel_array(x,y,xi,eta,xp,yp,kb,config)
    rhs_values=(2*PI/source_count)*kernel@source_values
    residuals=0.5*u_targets-rhs_values

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

    norms_left=[];norms_right=[]
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
                norms_left.append(field_norm)
                coeffs_left.append(coeff)
                ratios_left.append(float(ratio))
            else:
                norms_right.append(field_norm)
                coeffs_right.append(coeff)
                ratios_right.append(float(ratio))

    def phase_consistency(coeffs,side):
        if not coeffs:return math.nan,complex(math.nan,math.nan)
        arr=np.asarray(coeffs,complex)
        dephased=arr*np.exp(-1j*q1*distances)
        variation=float(np.max(abs(dephased-np.mean(dephased)))/max(float(np.max(abs(dephased))),1e-30))
        return variation,complex(np.mean(dephased))

    left_consistency, left_signed = phase_consistency(coeffs_left, "left")
    right_consistency, right_signed = phase_consistency(coeffs_right, "right")

    modal_left,modal_right=modal_open_channel_coefficients(kb,epsilon,a,vector,config)
    discrepancies=[]
    for coeffs,norms,modal in [(coeffs_left,norms_left,modal_left),(coeffs_right,norms_right,modal_right)]:
        for distance,coeff,norm in zip(distances,coeffs,norms):
            discrepancies.append(abs(coeff*np.exp(-1j*q1*distance)-modal)*math.sqrt(norm_phi_sq)/max(norm,1e-30))
    return {
        "modal_projection_discrepancy":float(max(discrepancies,default=math.inf)),
        "modal_coefficient_left":modal_left,"modal_coefficient_right":modal_right,
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
    _, _, _, _, vector = full_singular_pair(
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
        modal_projection_discrepancy=channel["modal_projection_discrepancy"],
        modal_projection_verified=channel["modal_projection_discrepancy"]<=config.physical_modal_projection_tolerance,
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


def refine_k_bracket(left,right,epsilon,a,config):
    # Work in sigma: bounded optimization in k loses small-cutoff separations
    # through its sqrt(machine-epsilon)*abs(k) stopping floor.
    sigma_lo=sigma_from_kb(right,config);sigma_hi=sigma_from_kb(left,config)
    def measure(sig):
        kb=kb_from_sigma(float(sig),config)
        _,_,rel,scaled=full_singular_metrics(kb,epsilon,a,config.bem_M,config)
        value=max(rel,scaled)
        if not np.isfinite(value):raise ValueError("Nonfinite refinement objective")
        return value
    endpoints=[measure(sigma_lo),measure(sigma_hi)]
    result=minimize_scalar(lambda sig:math.log10(max(measure(sig),np.finfo(float).tiny)),
            bounds=(sigma_lo,sigma_hi),method='bounded',options={'xatol':config.minimizer_xatol})
    if not result.success:raise RuntimeError('Refinement failed: '+str(result.message))
    sigma=float(result.x)
    if not sigma_lo<sigma<sigma_hi:raise RuntimeError("Minimum is not interior")
    kb=kb_from_sigma(sigma,config)
    smin,_,rel,scaled,vector=full_singular_pair(kb,epsilon,a,config.bem_M,config)
    drop=min(endpoints)/max(rel,scaled,np.finfo(float).tiny)
    if not np.all(np.isfinite([sigma,kb,smin,rel,scaled,drop])):raise ValueError("Incomplete candidate metrics")
    return kb,sigma,smin,rel,scaled,drop,vector


def candidate_matches_primary(kb,sigma,vector,epsilon,mode,config):
    # scipy's bounded refinement has a sqrt(eps)*sigma stopping floor even
    # when xatol is smaller. Propagate it to kb, rather than using a fixed
    # exclusion radius. An independently recomputed vector must also agree.
    resolution=max(10*config.minimizer_xatol,
                   8*np.sqrt(np.finfo(float).eps)*config.b**2*sigma**2/max(abs(kb),1e-30),
                   5*mode.normalized_residual*config.b**2*sigma**2/max(abs(kb),1e-30))
    if abs(kb-mode.kb_numerical)>resolution:return False
    *_,primary=full_singular_pair(mode.kb_numerical,epsilon,mode.a_numerical,config.bem_M,config)
    overlap=abs(np.vdot(primary,vector))/(np.linalg.norm(primary)*np.linalg.norm(vector))
    return bool(np.isfinite(overlap) and overlap>=1-1e-7)


def whole_band_uniqueness_screen(epsilon,mu,mode,config):
    """Converged sampled screen, never a theorem-wide uniqueness certificate."""
    all_candidates=[];diagnostics={'levels':[],'complete':False,'sampled_screen_passed':False,
        'uniqueness_supported':False,'limitation':'Finite grid and omitted cutoff intervals; no certified global exclusion.'}
    previous_roots=None;stable=False
    scan=np.array([],float);values=np.array([],float)
    for level in config.uniqueness_refinement_levels:
        cfg=replace(config,uniqueness_linear_points=config.uniqueness_linear_points*level,
                    uniqueness_log_points=config.uniqueness_log_points*level)
        scan=build_uniqueness_grid(epsilon,mu,cfg);values=np.full(len(scan),math.nan)
        errors=[];rows=[];roots=[]
        print(f'  sampled real-band screen: level={level}, points={len(scan)}')
        for i,kb in enumerate(scan):
            try:
                _,_,rel,scaled=full_singular_metrics(float(kb),epsilon,mode.a_numerical,config.bem_M,config)
                values[i]=max(rel,scaled)
                if not np.isfinite(values[i]):errors.append(f'nonfinite@{kb:.16g}')
            except Exception as exc:errors.append(f'{type(exc).__name__}@{kb:.16g}: {exc}')
        brackets=sampled_local_minimum_brackets(scan,values)
        for index,(left,right) in enumerate(brackets):
            base=dict(epsilon=float(epsilon),a=float(mode.a_numerical),kb=(left+right)/2,
                sigma=math.nan,relative_singular_value=math.nan,scaled_singular_value=math.nan,
                drop_factor=math.nan,boundary_even_residual=math.nan,open_channel_relative_amplitude=math.nan,
                nonradiating_resolved_candidate=False,scan_level=level)
            if index>=config.uniqueness_max_candidates:
                rows.append(AdditionalCandidate(**base,status='pending/candidate-limit',reason='Not evaluated; blocks screen approval'))
                continue
            try:
                kb,sigma,_,rel,scaled,drop,v=refine_k_bracket(left,right,epsilon,mode.a_numerical,config)
                base.update(kb=kb,sigma=sigma,relative_singular_value=rel,scaled_singular_value=scaled,
                            drop_factor=drop,boundary_even_residual=boundary_x_even_residual(v))
                cert=certify_full_mode(kb,epsilon,mode.a_numerical,config.bem_M,True,drop,config)
                # Local deduplication AFTER refinement, not a broad pre-exclusion.
                known=candidate_matches_primary(kb,sigma,v,epsilon,mode,config)
                if known:
                    status='known-root';reason='Matches independently certified primary root within local residual uncertainty'
                elif max(rel,scaled)>max(config.relative_near_singular_tolerance,config.scaled_smin_tolerance):
                    status='excluded-by-sampled-residual';reason='Refined sampled valley is not near singular'
                elif not cert['supported']:
                    status='ambiguous';reason='Near singular valley lacks complete spectral certificate'
                else:
                    candidate_mode=replace(mode,kb_numerical=kb,sigma_numerical=sigma)
                    phys=physical_mode_diagnostics(epsilon,candidate_mode,config)
                    base['open_channel_relative_amplitude']=phys.first_open_channel_relative_amplitude
                    # Global uniqueness must allow both x-parities.
                    physical=bool(phys.open_channel_suppressed and phys.open_channel_phase_consistent
                        and phys.walls_verified and phys.boundary_integral_verified and phys.decay_verified
                        and phys.modal_projection_verified)
                    if physical:
                        status='additional-mode-supported';reason='Full spectral and physical checks passed'
                        base['nonradiating_resolved_candidate']=True;roots.append(kb)
                    elif not phys.open_channel_suppressed and phys.modal_projection_verified:
                        status='radiating-candidate';reason='Independently checked open channel is nonzero at configured resolution'
                    else:status='ambiguous';reason='Physical diagnostics incomplete or inconclusive'
                rows.append(AdditionalCandidate(**base,status=status,reason=reason,optimizer_success=True,minimum_is_interior=True))
            except Exception as exc:
                rows.append(AdditionalCandidate(**base,status='ambiguous/refinement-failed',reason=f'{type(exc).__name__}: {exc}'))
        unresolved=[r for r in rows if r.status.startswith(('ambiguous','pending'))]
        complete=bool(len(scan)>=3 and not errors and not unresolved and np.all(np.isfinite(values)))
        stable=bool(previous_roots is not None and len(previous_roots)==len(roots)
                    and (not roots or np.max(abs(np.sort(previous_roots)-np.sort(roots)))<=config.uniqueness_candidate_stability_tolerance))
        diagnostics['levels'].append(dict(level=level,points=len(scan),brackets=len(brackets),
            evaluated=min(len(brackets),config.uniqueness_max_candidates),unresolved=len(unresolved),
            known_roots=sum(r.status=='known-root' for r in rows),additional_modes=roots,
            errors=errors,complete=complete,candidate_sets_stable=stable))
        all_candidates.extend(rows);previous_roots=roots
    diagnostics['complete']=bool(diagnostics['levels'] and all(r['complete'] for r in diagnostics['levels']))
    diagnostics['sampled_screen_passed']=bool(diagnostics['complete'] and stable and not previous_roots)
    diagnostics['candidate_sets_stable']=stable
    diagnostics['lower_omitted_interval']=[kb_1(config),kb_1(config)+config.uniqueness_lower_k_margin]
    diagnostics['upper_omitted_interval']=[kb_2(config)-config.uniqueness_upper_k_margin,kb_2(config)]
    return ScreenResult(scan,values,all_candidates,diagnostics)


def physical_is_supported(physical):
    if physical is None:return False
    checks=('boundary_x_even_verified','field_x_even_verified','open_channel_suppressed',
            'open_channel_phase_consistent','walls_verified','boundary_integral_verified','decay_verified',
            'modal_projection_verified')
    metrics=('boundary_x_even_residual','field_x_even_residual','first_open_channel_relative_amplitude',
             'wall_relative_residual','boundary_integral_relative_residual','decay_relative_error_left',
             'decay_relative_error_right','modal_projection_discrepancy')
    return bool(all(getattr(physical,k,False) is True for k in checks)
                and np.all(np.isfinite([getattr(physical,k,math.nan) for k in metrics])))


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
        physical_ok = physical_is_supported(physical)

    run_scan = bool(
        config.run_uniqueness_screen
        and epsilon <= config.uniqueness_epsilon_max + 1.0e-15
        and mode.resolved
        and physical_ok is True
    )

    if run_scan:
        screen = whole_band_uniqueness_screen(
            epsilon, mu, mode, config
        )
        scan,scan_values,additional=screen
        screen_diagnostics=screen.diagnostics
        additional_resolved = sum(
            c.nonradiating_resolved_candidate
            for c in additional if c.scan_level==config.uniqueness_refinement_levels[-1]
        )
        uniqueness_passed: bool | None = True if screen_diagnostics["sampled_screen_passed"] else None
    else:
        scan = np.array([], dtype=float)
        scan_values = np.array([], dtype=float)
        additional = []
        additional_resolved = 0
        uniqueness_passed = None
        screen_diagnostics={"complete":False,"reason":"Screen not run"}

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
        scaled_singular_value_final=float(mode.scaled_singular_value),
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
        mode_exists_supported=bool(mode.resolved and physical_ok is True),
        uniqueness_supported=False,uniqueness_status="ambiguous/not-established",
        screen_diagnostics_json=json.dumps(screen_diagnostics),
        normalized_residual=mode.normalized_residual,equilibrated_residual=mode.equilibrated_residual,
        equilibrated_sigma_min=mode.equilibrated_sigma_min,
    )

    directory=Path(config.output_directory)/'certificates'/f'epsilon_{epsilon:.8f}'
    harness.atomic_json(directory/'mode.json',asdict(mode))
    harness.atomic_json(directory/'certificate.json',asdict(result))
    harness.atomic_json(directory/'screen.json',screen_diagnostics)
    if physical is not None:harness.atomic_json(directory/'physical.json',asdict(physical))
    return result, mode, physical, scan, scan_values, additional


def build_placement_a1_fits(
    summaries: list[ValidationResult],
    a1_formula: float,
    config: Config,
) -> list[PlacementFitDiagnostics]:
    """Estimate a1 independently from the numerically located placement branch.

    The theorem implies a_BIC/epsilon = a1 + O(epsilon).  Therefore the
    epsilon->0 intercept of a_num/epsilon supplies a numerical estimate of a1
    that is independent of the MFS evaluation of formula (2.7).  We report a
    provisional fit using all resolved full-matrix roots and, when enough points
    exist, a stricter fit using only physically verified BICs.
    """
    fits: list[PlacementFitDiagnostics] = []

    for selection, physical_only in (("resolved-roots", False), ("physical-bics", True)):
        rows = [
            r for r in summaries
            if r.epsilon <= config.placement_fit_epsilon_max + 1.0e-15
            and r.expected_bic_resolved
            and np.isfinite(r.a_numerical)
            and r.epsilon > 0.0
            and (not physical_only or r.physical_bic_verified is True)
        ]
        if len(rows) < config.placement_fit_min_points:
            continue

        rows = sorted(rows, key=lambda r: r.epsilon)
        x = np.asarray([r.epsilon for r in rows], dtype=float)
        y = np.asarray([r.a_numerical / r.epsilon for r in rows], dtype=float)

        linear = np.polyfit(x, y, 1)
        slope = float(linear[0])
        intercept = float(linear[1])
        pred = slope * x + intercept
        rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
        rel = abs(intercept - a1_formula) / max(abs(a1_formula), 1.0e-30)

        if len(rows) >= max(4, config.placement_fit_min_points + 1):
            quadratic = np.polyfit(x, y, 2)
            q_intercept = float(quadratic[2])
            q_rel = abs(q_intercept - a1_formula) / max(abs(a1_formula), 1.0e-30)
        else:
            q_intercept = math.nan
            q_rel = math.nan

        fits.append(
            PlacementFitDiagnostics(
                selection=selection,
                point_count=len(rows),
                epsilon_max=float(max(x)),
                a1_formula_2_7=float(a1_formula),
                a1_fit_linear=intercept,
                a2_fit_linear=slope,
                linear_rmse=rmse,
                linear_relative_difference=float(rel),
                a1_fit_quadratic=q_intercept,
                quadratic_relative_difference=float(q_rel),
            )
        )

    return fits


# ---------------------------------------------------------------------------
# Internal convergence at fixed M
# ---------------------------------------------------------------------------


def run_internal_convergence_study(summaries,a1,config):
    """Retune BOTH placement and frequency, then repeat physical certificates."""
    if not config.run_internal_convergence_study:return []
    target=min(summaries,key=lambda r:abs(r.epsilon-config.internal_convergence_epsilon))
    eps=target.epsilon;mu=target.mu;rows=[]
    reference=tune_embedded_mode(eps,mu,a1,config)
    variations=[('finite_difference_step',h,replace(config,finite_difference_step=h)) for h in config.finite_difference_steps_test]
    variations += [('lattice_terms',float(n),replace(config,lattice_terms=n)) for n in config.lattice_terms_test]
    variations += [('harmonic_order',float(n),replace(config,harmonic_order=n)) for n in config.harmonic_orders_test]
    for name,value,cfg in variations:
        row=InternalConvergenceRow(eps,name,float(value),math.nan,math.nan,math.nan,math.nan,math.nan,math.nan)
        try:
            mode=tune_embedded_mode(eps,mu,a1,cfg,previous_mode=reference)
            physical=physical_mode_diagnostics(eps,mode,cfg) if mode.resolved else None
            verified=physical_is_supported(physical)
            row.kb=mode.kb_numerical;row.sigma_bem=mode.sigma_numerical
            row.relative_singular_value=mode.relative_singular_value;row.scaled_singular_value=mode.scaled_singular_value
            row.relative_kb_shift_from_reference=abs(row.kb-reference.kb_numerical)/abs(reference.kb_numerical)
            row.relative_sigma_shift_from_reference=abs(row.sigma_bem-reference.sigma_numerical)/abs(reference.sigma_numerical)
            row.a_numerical=mode.a_numerical;row.a_shift_from_reference=abs(mode.a_numerical-reference.a_numerical)
            row.physical_bic_verified=verified;row.normalized_residual=mode.normalized_residual
            row.open_channel_relative_amplitude=physical.first_open_channel_relative_amplitude if physical else math.nan
            row.status='certified' if verified and mode.resolved else 'ambiguous/not-established'
        except Exception as exc:row.status=f'ambiguous/{type(exc).__name__}: {exc}'
        rows.append(row)
        print(f'  retuned convergence {name}={value:g}: {row.status}',flush=True)
    return rows



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
    plt.plot(N, [abs(r.mu-rows[-1].mu) for r in rows], "o-")
    plt.yscale("symlog",linthresh=1e-13)
    plt.xlabel("MFS order")
    plt.ylabel("|mu_N - mu_last| (zero retained)")
    plt.tight_layout()
    plt.savefig(output / "mu_mfs_convergence.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(N, [r.a1_formula_2_7 for r in rows], "o-", label="formula (2.7)")
    plt.plot(N, [r.a1_formula_5_47 for r in rows], "s--", label="formula (5.47)")
    plt.axhline(
        -config.shape_beta / 12.0,
        linestyle="--",
        label=r"Printed example: $-\beta/12$ (inconsistent; not a target)",
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
        [
            max(r.relative_singular_value_final, r.scaled_singular_value_final, 1e-30)
            for r in results
        ],
        "o-",
    )
    plt.axhline(
        max(config.relative_near_singular_tolerance, config.scaled_smin_tolerance),
        linestyle="--",
        label="auxiliary tolerance only (not the full certificate)",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel("auxiliary combined singularity measure")
    plt.title("Auxiliary diagnostic; normalized/equilibrated certificates required")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "relative_singular_value.png", dpi=config.plot_dpi)
    plt.close()

    plt.figure(figsize=(8,5))
    for field,label,tolerance in [("normalized_residual","normalized",config.singular_vector_residual_tolerance),
                                   ("equilibrated_residual","row-equilibrated",config.equilibrated_residual_tolerance),
                                   ("equilibrated_sigma_min","scaled smin",config.equilibrated_smin_tolerance)]:
        plt.semilogy(eps,[getattr(r,field)/tolerance for r in results],"o-",label=label)
    plt.axhline(1.,color="black",linestyle="--",label="acceptance limit")
    plt.xlabel("epsilon");plt.ylabel("metric / its configured tolerance")
    plt.title("Spectral certificates (physical checks are additional)")
    plt.legend();plt.tight_layout();plt.savefig(output/"effective_spectral_certificates.png",dpi=config.plot_dpi);plt.close()


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
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    refined=max(result.relative_singular_value_final,result.scaled_singular_value_final)
    for ax in axes:
        ax.semilogy(scan,values,"o-",markersize=3,label="sampled grid")
        ax.semilogy([result.kb_numerical],[refined],"r*",markersize=12,label="refined certified root")
        ax.axvline(result.kb_numerical,linestyle="--",color="grey")
        ax.set_xlabel("kb");ax.set_ylabel("auxiliary singularity measure");ax.legend()
    axes[0].axvspan(kb_1(config),scan[0],color="grey",alpha=.2)
    axes[0].axvspan(scan[-1],kb_2(config),color="grey",alpha=.2,label="omitted margins")
    delta=max(5*abs(kb_2(config)-result.kb_numerical),1e-6)
    axes[1].set_xlim(max(kb_1(config),result.kb_numerical-delta),min(kb_2(config),result.kb_numerical+delta))
    axes[1].set_title("Near-cutoff zoom; line joins samples only")
    fig.suptitle(f"Sampled-band screen, epsilon={epsilon:.2f}; global uniqueness not established")
    fig.tight_layout();fig.savefig(output/f"whole_band_epsilon_{epsilon:.2f}.png",dpi=config.plot_dpi);plt.close(fig)


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
        plt.yscale("symlog",linthresh=1e-12)
        passed=[r for r in subset if r.physical_bic_verified and r.status=="certified"]
        failed=[r for r in subset if r not in passed]
        plt.scatter([r.value for r in passed],[r.relative_sigma_shift_from_reference for r in passed],marker="o",color="green",zorder=4,label="certified")
        plt.scatter([r.value for r in failed],[r.relative_sigma_shift_from_reference for r in failed],marker="x",s=70,color="red",zorder=5,label="ambiguous / failed certificate")
        plt.legend()
        plt.xlabel(parameter)
        plt.ylabel("relative sigma shift (zero retained)")
        plt.title("Frequency stability alone does not certify a mode")
        plt.tight_layout()
        plt.savefig(
            output / f"internal_convergence_{parameter}.png",
            dpi=config.plot_dpi,
        )
        plt.close()


def plot_placement_a1_fit(
    summaries: list[ValidationResult],
    a1: float,
    output: Path,
    config: Config,
) -> None:
    rows = [
        r for r in summaries
        if r.epsilon <= config.placement_fit_epsilon_max + 1.0e-15
        and r.expected_bic_resolved
        and np.isfinite(r.a_numerical)
        and r.epsilon > 0.0
    ]
    if not rows:
        return
    rows = sorted(rows, key=lambda r: r.epsilon)
    x = np.asarray([r.epsilon for r in rows], dtype=float)
    y = np.asarray([r.a_numerical / r.epsilon for r in rows], dtype=float)

    plt.figure(figsize=(7.5, 4.8))
    plt.plot(x, y, "o", label="resolved full-matrix roots")
    physical = [r for r in rows if r.physical_bic_verified is True]
    if physical:
        xp = np.asarray([r.epsilon for r in physical], dtype=float)
        yp = np.asarray([r.a_numerical / r.epsilon for r in physical], dtype=float)
        plt.plot(xp, yp, "s", label="physical BIC PASS")
    if len(rows) >= config.placement_fit_min_points:
        slope, intercept = np.polyfit(x, y, 1)
        xx = np.linspace(0.0, max(x), 200)
        plt.plot(xx, slope * xx + intercept, "--", label="linear epsilon->0 fit")
    plt.axhline(a1, linestyle=":", label="a1 from formula (2.7)")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$a_{\rm num}/\varepsilon$")
    plt.title("Independent placement-branch estimate of a1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "placement_a1_fit.png", dpi=config.plot_dpi)
    plt.close()


def build_publication_claims(
    shape: ShapeDiagnostics,
    mu: float,
    nu: float,
    a1: float,
    a1_robustness: list[A1RobustnessRow],
    placement_fits: list[PlacementFitDiagnostics],
    small_beta_fit: SmallBetaFitDiagnostics | None,
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

    a1_worst = max(
        (r.relative_difference_from_primary for r in a1_robustness),
        default=0.0,
    )
    a1_robust_ok = bool(
        a1_robustness
        and np.all(np.isfinite([r.relative_difference_from_primary for r in a1_robustness]))
        and a1_worst <= config.a1_robustness_relative_tolerance
    )
    claims.append(
        PublicationClaimRow(
            claim="Direct numerical evaluation of a1 from theorem formula (2.7)",
            status="supported" if a1_robust_ok else "partial",
            evidence=(
                f"a1={a1:.12g} from the exterior-Laplace MFS; "
                f"max source-scale/quadrature relative variation={a1_worst:.3e}."
            ),
            limitation=(
                "The (5.47) value is an algebraic cross-check using the same Psi, not a fully "
                "independent PDE discretization. The BIC-branch intercept below is the independent check."
            ),
        )
    )

    small = [
        r for r in summaries
        if r.epsilon <= config.publication_asymptotic_epsilon_max + 1e-15
    ]
    physical_small = [
        r for r in small
        if r.expected_bic_resolved and r.physical_bic_verified is True
    ]

    physical_fit = next(
        (f for f in placement_fits if f.selection == "physical-bics"),
        None,
    )
    resolved_fit = next(
        (f for f in placement_fits if f.selection == "resolved-roots"),
        None,
    )
    fit_for_evidence = physical_fit if physical_fit is not None else resolved_fit
    fit_ok = bool(
        physical_fit is not None
        and physical_fit.linear_relative_difference <= config.placement_fit_relative_tolerance
    )

    correction_text = (
        f"scaled correction range "
        f"[{min(r.scaled_a_correction for r in physical_small):.3g}, "
        f"{max(r.scaled_a_correction for r in physical_small):.3g}]"
        if physical_small else "no physically verified small-epsilon placement data"
    )
    if fit_for_evidence is not None:
        fit_text = (
            f"{fit_for_evidence.selection} extrapolation gives "
            f"a1={fit_for_evidence.a1_fit_linear:.12g} "
            f"({fit_for_evidence.linear_relative_difference:.3%} from formula (2.7))"
        )
    else:
        fit_text = "insufficient verified points for an independent epsilon->0 intercept"

    placement_ok = bool(
        len(small) >= config.publication_min_asymptotic_points
        and len(physical_small) == len(small)
        and all(r.placement_minimum_is_interior for r in physical_small)
        and np.all(np.isfinite([r.scaled_a_correction for r in physical_small]))
        and a1_robust_ok
        and fit_ok
    )
    claims.append(
        PublicationClaimRow(
            claim=r"Finite-sample consistency with placement law a = epsilon a1 + O(epsilon^2)",
            status="supported" if placement_ok else "partial",
            evidence=(
                f"a1 from formula (2.7)={a1:.12g}; {correction_text}; {fit_text}."
            ),
            limitation=(
                "The full-BEM placement intercept is independent of the MFS a1 calculation, but "
                "finite epsilon sampling still cannot prove the O(epsilon^2) remainder."
            ),
        )
    )

    claims.append(
        PublicationClaimRow(
            claim="One non-radiating embedded trapped mode for sampled tuned geometries",
            status=(
                "supported"
                if len(small) >= config.publication_min_asymptotic_points
                and len(physical_small) == len(small)
                else "partial"
            ),
            evidence=(
                f"{len(physical_small)}/{len(small)} small-epsilon cases pass full-matrix "
                "singularity (relative and row-scaled smin), x-even field, signed open-channel "
                "suppression/phase consistency, off-grid BIE and decay checks."
            ),
            limitation=(
                "Finite sampling supports the theorem numerically; it does not prove existence "
                "for every sufficiently small epsilon."
            ),
        )
    )

    if physical_small:
        C = asymptotic_coefficient(mu, config)
        ratio_errors = np.array(
            [
                abs(r.sigma_over_epsilon_squared - C) / max(abs(C), 1e-30)
                for r in physical_small
            ]
        )
        q = np.array([r.scaled_sigma_remainder for r in physical_small])
        scaling_ok = bool(
            len(physical_small) >= config.publication_min_asymptotic_points
            and np.all(np.isfinite(q))
            and ratio_errors[0] <= ratio_errors[-1]
        )
        claims.append(
            PublicationClaimRow(
                claim=r"Finite-sample consistency with sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)",
                status="supported" if scaling_ok else "partial",
                evidence=(
                    f"Using physically verified BICs only, relative error of sigma/epsilon^2 "
                    f"from C changes from {ratio_errors[0]:.3%} at epsilon={physical_small[0].epsilon:g} "
                    f"to {ratio_errors[-1]:.3%} at epsilon={physical_small[-1].epsilon:g}; "
                    f"scaled-remainder range [{np.min(q):.3g},{np.max(q):.3g}]."
                ),
                limitation=(
                    "Finite asymptotic data support the stated scaling but do not prove the remainder bound."
                ),
            )
        )
    else:
        claims.append(
            PublicationClaimRow(
                claim=r"Finite-sample consistency with sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)",
                status="not-established",
                evidence="No physically verified small-epsilon BICs available for the asymptotic comparison.",
                limitation="Spectral near-roots alone are not used to certify a trapped-mode asymptotic law.",
            )
        )

    scanned = [r for r in small if r.uniqueness_screen_run]
    scan_ok = bool(
        len(scanned)==len(small) and len(small)>0
        and all(r.physical_bic_verified is True for r in scanned)
        and all(r.uniqueness_screen_passed is True for r in scanned)
    )
    claims.append(
        PublicationClaimRow(
            claim="No additional resolved non-radiating candidate on the sampled real band",
            status="supported" if scan_ok else "partial",
            evidence=(
                f"{sum(r.uniqueness_screen_passed is True for r in scanned)}/"
                f"{len(scanned)} scanned small-epsilon cases found no additional candidate "
                "passing the strengthened singularity and open-channel filters."
            ),
            limitation=(
                "This is a finite real-axis screen and must not be interpreted as a proof of theorem-wide uniqueness."
            ),
        )
    )

    claims.append(PublicationClaimRow(claim="Uniqueness of the embedded mode on the full open band",
        status="ambiguous/not-established",evidence="A finite real-axis screen is not a certified spectral count.",
        limitation="Unscreened cutoff intervals and finite-grid exclusion are not covered by a global certificate."))
    if small_beta_fit is not None:
        example_ok = bool(small_beta_fit.relative_difference_from_paper <= 0.10)
        claims.append(
            PublicationClaimRow(
                claim="Example 2.4 / 5.4 small-beta coefficient a1 ~ -beta/12",
                status="supported" if example_ok else "unresolved-inconsistency",
                evidence=(
                    f"Direct (2.7) computations give fitted a1/beta slope "
                    f"{small_beta_fit.fitted_slope_through_origin:.8g} versus paper -1/12="
                    f"{small_beta_fit.paper_slope_minus_one_twelfth:.8g}."
                ),
                limitation=(
                    "Substituting the printed f0, f1 into (5.47) gives -beta/2, not -beta/12; "
                    "the MFS agrees with that substitution. The displayed example is internally inconsistent; "
                    "(2.7) remains the numerical definition, without forcing a printed reference value."
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
    a1_robustness: list[A1RobustnessRow],
    placement_fits: list[PlacementFitDiagnostics],
    small_beta_fit: SmallBetaFitDiagnostics | None,
    summaries: list[ValidationResult],
    claims: list[PublicationClaimRow],
    config: Config,
) -> None:
    final_mfs = mfs_rows[-1]
    robustness_worst = max(
        (r.relative_difference_from_primary for r in a1_robustness),
        default=math.nan,
    )
    lines = [
        "# Numerical validation of Theorem 2.3(iv)",
        "",
        "## Configuration",
        f"- assembly version = `{ASSEMBLY_VERSION}`",
        "- Result provenance: corrected/v3. Historical embedded_validation and v2 artifacts are legacy/pre-fix.",
        "- Existence and full-band uniqueness are separate. The finite screen never certifies full-band uniqueness.",
        "- Formula (2.7) is retained; the printed Example 5.4 coefficient disagrees with its own (5.47) substitution.",
        f"- b = {config.b}",
        f"- beta = {config.shape_beta}",
        f"- fixed BEM order M = {config.bem_M}",
        f"- epsilon values = {config.epsilon_values}",
        "- Green assembly uses exact vertical periodic reduction before evaluating the local lattice expansion.",
        "",
        "## Independent reference-obstacle quantities",
        f"- area S = {shape.area:.12g}",
        f"- mu = {final_mfs.mu:.12g}",
        f"- nu = {final_mfs.nu:.3e}",
        f"- a1 from theorem formula (2.7) = {final_mfs.a1_formula_2_7:.12g}",
        f"- a1 from algebraically equivalent formula (5.47) = {final_mfs.a1_formula_5_47:.12g}",
        f"- relative (2.7)/(5.47) mismatch = {final_mfs.a1_internal_relative_difference:.3e}",
        f"- a1 MFS source-scale/quadrature worst relative variation = {robustness_worst:.3e}",
        f"- Example 5.4 leading value -beta/12 = {final_mfs.a1_paper_example_leading:.12g}",
        f"- circle calibration: mu={circle.mu:.12g}, a1={circle.a1_formula_2_7:.3e}",
        "",
        "The value from (5.47) is an algebraic cross-check using the same computed Psi. "
        "The independent numerical check of a1 is the epsilon->0 placement-branch intercept below.",
        "",
        "## Independent a1 estimate from the BIC placement branch",
    ]

    if placement_fits:
        lines.extend([
            "| selection | points | eps max | a1 intercept | rel. diff. vs (2.7) | linear RMSE | quadratic intercept |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ])
        for fit in placement_fits:
            q = f"{fit.a1_fit_quadratic:.12g}" if np.isfinite(fit.a1_fit_quadratic) else "--"
            lines.append(
                f"| {fit.selection} | {fit.point_count} | {fit.epsilon_max:.3g} | "
                f"{fit.a1_fit_linear:.12g} | {fit.linear_relative_difference:.3%} | "
                f"{fit.linear_rmse:.3e} | {q} |"
            )
    else:
        lines.append("Insufficient resolved small-epsilon placements for an independent a1 intercept fit.")

    lines.extend(["", "## Small-beta Example 5.4 diagnostic"])
    if small_beta_fit is not None:
        lines.append(
            f"- fitted slope a1/beta = {small_beta_fit.fitted_slope_through_origin:.12g}; "
            f"paper slope -1/12 = {small_beta_fit.paper_slope_minus_one_twelfth:.12g}; "
            f"relative difference = {small_beta_fit.relative_difference_from_paper:.3%}."
        )
        lines.append(
            "- This discrepancy is reported as unresolved; the script does not replace formula (2.7) by -beta/12."
        )
    else:
        lines.append("Small-beta diagnostic disabled or unavailable.")

    lines.extend([
        "",
        "## Methodological note",
        (
            "Unlike Theorem 2.3(iii), no parity reduction removes the continuous spectrum here. "
            "The full BEM matrix is used. Formula (2.7) supplies only the leading placement "
            "coefficient: a=epsilon*a1+c*epsilon^2, with c determined numerically. The final "
            "acceptance requires both relative and row-scaled singularity criteria plus physical "
            "non-radiation diagnostics."
        ),
        "",
        "## Claim audit",
    ])
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
            "| eps | a_num | (a-eps a1)/eps^2 | kb | sigma | rel sv | row-scaled smin | physical |",
            "|---:|---:|---:|---:|---:|---:|---:|:---:|",
        ]
    )
    for r in summaries:
        lines.append(
            f"| {r.epsilon:.3f} | {r.a_numerical:.8g} | "
            f"{r.scaled_a_correction:.6g} | {r.kb_numerical:.10f} | "
            f"{r.sigma_numerical:.6g} | {r.relative_singular_value_final:.3e} | "
            f"{r.scaled_singular_value_final:.3e} | "
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
    print(
        "row-scaled smallest singular value = "
        f"{result.scaled_singular_value_final:.3e}"
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
            and physical.open_channel_phase_consistent
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
            f"{'PASS (sampled screen only)' if result.uniqueness_screen_passed else 'AMBIGUOUS'} "
            f"(additional resolved BICs={result.additional_resolved_bics})"
        )
    else:
        print("sampled-band screen = not run for this epsilon")
    print(f"mode_exists_supported={result.mode_exists_supported}; uniqueness_supported={result.uniqueness_supported} ({result.uniqueness_status})")


def _run_validation(config: Config) -> None:
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
    print(f"assembly version = {ASSEMBLY_VERSION}")
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

    print("\n=== a1 ROBUSTNESS: SOURCE SCALE / QUADRATURE ===")
    a1_robustness = run_a1_robustness_study(a1, config)
    write_dataclass_csv(output / "a1_robustness.csv", a1_robustness)
    if a1_robustness:
        worst_a1 = max(r.relative_difference_from_primary for r in a1_robustness)
        print(f"max relative a1 variation = {worst_a1:.3e}")
        for row in a1_robustness:
            print(
                f"  rho={row.source_scale:.2f}, Q={row.quadrature_points:>5d}: "
                f"a1={row.a1_formula_2_7:.12g}, rel={row.relative_difference_from_primary:.3e}, "
                f"MFS_res={row.boundary_relative_residual:.3e}"
            )
    else:
        print("a1 robustness study disabled")

    beta_rows: list[MFSRow] = []
    small_beta_fit: SmallBetaFitDiagnostics | None = None
    if config.run_small_beta_a1_diagnostic:
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
        small_beta_fit = fit_small_beta_a1(beta_rows)
        if small_beta_fit is not None:
            write_dataclass_csv(output / "a1_small_beta_fit.csv", [small_beta_fit])
            print(
                "  fitted small-beta slope a1/beta = "
                f"{small_beta_fit.fitted_slope_through_origin:.10g}; "
                f"paper -1/12={small_beta_fit.paper_slope_minus_one_twelfth:.10g}; "
                f"relative difference={small_beta_fit.relative_difference_from_paper:.3%}"
            )

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
        ) = harness.case(
            sys.modules[__name__], config, f"epsilon_{epsilon:.8f}",
            (epsilon, shape, mu, nu, a1), {"previous_mode": previous_mode},
        )

        summaries.append(result)
        modes.append(mode)
        if result.mode_exists_supported:
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

    placement_fits = build_placement_a1_fits(summaries, a1, config)
    write_dataclass_csv(output / "placement_a1_fit.csv", placement_fits)
    plot_placement_a1_fit(summaries, a1, output, config)
    print("\n=== INDEPENDENT a1 FROM NUMERICAL BIC PLACEMENT BRANCH ===")
    if placement_fits:
        for fit in placement_fits:
            print(
                f"{fit.selection}: n={fit.point_count}, "
                f"a1_fit={fit.a1_fit_linear:.12g}, "
                f"formula(2.7)={fit.a1_formula_2_7:.12g}, "
                f"relative difference={fit.linear_relative_difference:.3%}, "
                f"RMSE={fit.linear_rmse:.3e}"
            )
    else:
        print("insufficient resolved points for a placement-branch a1 fit")

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
        a1_robustness,
        placement_fits,
        small_beta_fit,
        summaries,
        config,
    )
    write_dataclass_csv(output / "publication_claims.csv", claims)
    write_publication_report(
        output,
        shape,
        mfs_rows,
        circle,
        a1_robustness,
        placement_fits,
        small_beta_fit,
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


def main() -> None:
    harness.run_entrypoint(sys.modules[__name__], CONFIG, _run_validation)


if __name__ == "__main__":
    main()
