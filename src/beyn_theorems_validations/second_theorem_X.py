from __future__ import annotations

"""
Numerical validation of Theorem 2.3(iii) from

    P. Zhevandrov et al.,
    "Discrete and embedded trapped modes in a plane quantum waveguide
     with a small obstacle: exact solutions" (March 4, 2025).

Target statement (Theorem 2.3(iii))
-----------------------------------
For an obstacle symmetric with respect to the x-axis,

    a = 0,
    X(t) even,
    Y(t) odd,

there is, for sufficiently small epsilon, a unique eigenvalue embedded in
[Lambda_1, Lambda_2),

    k^2 = Lambda_2 - sigma^2,

with

    sigma = epsilon^2 * pi^3 * mu / b^3 + O(epsilon^3 log epsilon).

The key numerical idea is to exploit the exact odd-in-y symmetry mentioned in
Remark 2.5 of the paper.  In the odd subspace the first transverse (even) open
channel disappears, so [Lambda_1, Lambda_2) becomes a discrete spectral window.
We therefore:

  1. assemble the same full-strip BEM operator as in the Theorem 2.1 code;
  2. project it exactly onto the odd-in-y parity subspace;
  3. use Beyn on that reduced operator for global discovery/counting;
  4. refine each candidate in sigma with fixed-M SVD;
  5. reconstruct the full odd field and verify wall values, odd symmetry,
     suppression of the first open channel, the off-grid BIE, and exponential
     decay with rate sigma;
  6. compute the shape dipole strength mu independently from the inflated
     obstacle via a Laplace exterior-Neumann Method of Fundamental Solutions
     (MFS), including a circle calibration and a convergence table;
  7. compare sigma/epsilon^2 with pi^3 mu / b^3 and report the scaled remainder

         |sigma_num - C epsilon^2| / (epsilon^3 |log epsilon|).

This is a numerical validation of the numerically testable content.  It does
NOT constitute a proof of analyticity or theorem-wide uniqueness.

Project dependency
------------------
The script uses the same project-level ``lattice_sums.py`` module as the
existing Theorem 2.1 validation code.
"""

import csv
from functools import lru_cache
import execution_harness as harness
import math
import os
import sys
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, ClassVar, cast

# Avoid BLAS oversubscription on macOS / multi-process project environments.
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
from scipy.optimize import minimize_scalar  # noqa: E402

# Same import convention as the existing waveguide scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice  # noqa: E402

PI = np.pi
ASSEMBLY_VERSION = "theorem-2.3-iii-periodic-green-v3"


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Config:
    # Waveguide / obstacle.
    b: float = 1.0
    a: float = 0.0  # Theorem 2.3(iii) requires a=0 exactly.
    shape_beta: float = 0.50

    # Physical small-obstacle sweep.  The first five points are the principal
    # small-epsilon reporting window; larger values are retained as finite-size
    # exploration.
    epsilon_values: tuple[float, ...] = (
        0.02, 0.04, 0.06, 0.08, 0.10,
        0.12, 0.14, 0.16, 0.18, 0.20,
    )
    publication_asymptotic_epsilon_max: float = 0.10
    publication_min_asymptotic_points: int = 4

    # Green function / BEM.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Independent computation of the dipole strength mu for the INFLATED
    # obstacle (epsilon removed).  The MFS sources are scaled toward the origin.
    mu_mfs_orders: tuple[int, ...] = (60, 80, 120, 160)
    mu_mfs_source_scale: float = 0.65
    mu_mfs_rcond: float = 1.0e-12
    mu_circle_calibration_tolerance: float = 1.0e-6
    mu_convergence_relative_tolerance: float = 5.0e-6

    # Beyn on the odd-in-y reduced operator.  M here is the FULL obstacle
    # boundary order.  The odd operator has dimension M/2.
    bem_M: int = 24
    beyn_quadrature_levels: tuple[int, ...] = (256, 512, 1024)
    beyn_adaptive_quadrature_levels: tuple[int, ...] = (2048, 4096)
    beyn_probe_dim: int = 8
    beyn_random_seed: int = 1729

    beyn_rank_relative_tolerance: float = 1.0e-8
    beyn_rank_absolute_tolerance: float = 1.0e-12
    beyn_rank_gap_threshold: float = 1.0e3
    beyn_empty_s0_tolerance: float = 1.0e-9

    # The contour encloses almost all of [sqrt(Lambda_1)b, sqrt(Lambda_2)b).
    # The left margin avoids evaluating exactly at Lambda_1 in the unreduced
    # full Green function.  The odd projection removes the open even channel.
    beyn_lower_threshold_margin: float = 1.0e-4
    beyn_cutoff_margin: float = 1.0e-5
    beyn_min_cutoff_margin: float = 1.0e-10
    beyn_cutoff_margin_fraction_of_predicted_gap: float = 0.25
    beyn_check_tighter_cutoff_margin: bool = True
    beyn_tighter_cutoff_margin_factor: float = 0.20
    beyn_imag_half_height: float = 3.0e-2

    beyn_real_axis_tolerance: float = 5.0e-5
    beyn_cluster_tolerance: float = 5.0e-5
    beyn_near_contour_real_tolerance: float = 5.0e-5
    beyn_candidate_convergence_tolerance: float = 2.0e-4
    beyn_s0_change_warning_tolerance: float = 0.15

    # Local SVD certification uses the single full boundary order bem_M.
    local_sigma_seed_bracket_factor: float = 4.0
    local_sigma_scan_points: int = 120
    minimizer_sigma_xatol: float = 1.0e-12
    minimum_drop_factor: float = 100.0
    near_singular_tolerance: float = 1.0e-6
    singular_vector_residual_tolerance: float = 1.0e-8
    equilibrated_residual_tolerance: float = 1.0e-8
    parity_operator_tolerance: float = 1.0e-7
    full_vector_residual_tolerance: float = 1.0e-8
    beyn_moment_convergence_tolerance: float = 1.0e-3
    beyn_linear_solve_tolerance: float = 1.0e-10
    mu_boundary_residual_tolerance: float = 1.0e-7
    mu_identity_relative_tolerance: float = 1.0e-6
    field_modal_relative_tolerance: float = 1.0e-12
    field_modal_terms: int = 128
    relative_near_singular_tolerance: float = 1.0e-4

    # Physical diagnostics for the embedded mode.
    run_physical_diagnostics: bool = True
    physical_quadrature_points: int = 96
    physical_x_over_b: tuple[float, ...] = (1.5, 2.0, 2.5, 3.0)
    physical_boundary_residual_samples: int = 128
    physical_wall_relative_tolerance: float = 1.0e-6
    physical_centerline_relative_tolerance: float = 1.0e-6
    physical_parity_relative_tolerance: float = 1.0e-6
    physical_open_channel_relative_tolerance: float = 1.0e-5
    physical_boundary_residual_tolerance: float = 5.0e-3
    physical_decay_relative_tolerance: float = 0.20

    # One representative discretization convergence study.
    run_internal_convergence_study: bool = True
    internal_convergence_epsilon: float = 0.10
    internal_sigma_bracket_factor: float = 1.20
    finite_difference_steps_test: tuple[float, ...] = (
        1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4,
    )
    lattice_terms_test: tuple[int, ...] = (25, 50, 100, 200)
    harmonic_orders_test: tuple[int, ...] = (5, 10, 15, 20, 30, 40)

    plot_dpi: int = 220
    output_directory: str = str(Path(__file__).resolve().parent / "results" / "theorem_2_3_iii" / "runs")
    execution_workers: int = 4
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
    x_even_defect: float
    y_odd_defect: float


@dataclass
class MuMFSRow:
    beta: float
    order: int
    source_scale: float
    mu: float
    nu: float
    boundary_relative_residual: float
    zero_total_source_error: float
    offgrid_relative_residual: float = math.nan
    identity_relative_error: float = math.nan


@dataclass
class BeynDiagnostics:
    epsilon: float
    full_M: int
    reduced_dimension: int
    quadrature_points: int
    estimated_rank: int
    threshold_rank: int
    gap_rank: int
    selected_gap_ratio: float
    leading_s0_singular_value: float
    max_linear_solve_relative_residual: float
    raw_eigenvalues: int
    accepted_candidates: int
    near_contour_seeds: int
    s0_relative_change: float = math.nan
    s1_relative_change: float = math.nan
    moments_converged: bool = False
    rank_stable: bool = False
    candidate_stable: bool = False
    _s0: ClassVar[Any] = None
    _s1: ClassVar[Any] = None


@dataclass
class BeynEigenvalueRow:
    epsilon: float
    quadrature_points: int
    raw_index: int
    real_part: float
    imag_part: float
    inside_contour: bool
    inside_target_band: bool
    near_real: bool
    real_span_overrun: float
    near_contour_seed: bool
    accepted_for_refinement: bool


@dataclass
class BeynConvergenceRow:
    epsilon: float
    quadrature_points: int
    estimated_rank: int
    s0_sv1: float
    s0_sv2: float
    s0_sv3: float
    s0_relative_change_from_previous: float
    accepted_candidates: int
    primary_candidate_real: float
    primary_candidate_imag: float
    max_candidate_shift_from_previous: float
    rank_stable_from_previous: bool
    s1_relative_change_from_previous: float = math.nan
    moments_converged: bool = False


@dataclass
class DiscoverySeed:
    real: float
    imag: float
    source: str


@dataclass
class ModeRefinementRow:
    epsilon: float
    candidate_index: int
    seed_source: str
    full_M: int
    reduced_dimension: int
    kb: float
    sigma_bem: float
    sigma_min: float
    sigma_max: float
    relative_singular_value: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_previous: float


@dataclass
class ModeResult:
    epsilon: float
    candidate_index: int
    seed_source: str
    kb_numerical: float
    sigma_numerical: float
    sigma_min_final: float
    sigma_max_final: float
    relative_singular_value_final: float
    final_drop_factor: float
    final_relative_mesh_change: float  # NaN for compatibility: no M sweep.
    resolved: bool
    normalized_residual: float = math.nan
    equilibrated_residual: float = math.nan
    equilibrated_sigma_min: float = math.nan
    full_vector_residual: float = math.nan
    parity_leakage_residual: float = math.nan


@dataclass
class PhysicalDiagnostics:
    epsilon: float
    kb: float
    full_M: int
    wall_relative_residual: float
    centerline_relative_residual: float
    odd_parity_relative_residual: float
    first_open_channel_relative_amplitude: float
    boundary_integral_relative_residual: float
    decay_rate_left: float
    decay_rate_right: float
    expected_decay_rate: float
    decay_relative_error_left: float
    decay_relative_error_right: float
    monotone_decay_left: bool
    monotone_decay_right: bool
    walls_verified: bool
    centerline_verified: bool
    parity_verified: bool
    open_channel_suppressed: bool
    boundary_integral_verified: bool
    decay_verified: bool


@dataclass
class ValidationResult:
    epsilon: float
    beta: float
    mu: float
    nu: float
    asymptotic_coefficient: float
    lambda_1: float
    lambda_2: float
    kb_lower_threshold: float
    kb_upper_threshold: float
    effective_cutoff_margin: float
    kb_asymptotic: float
    sigma_asymptotic: float
    beyn_final_quadrature_points: int
    beyn_estimated_rank: int
    beyn_rank_stable: bool
    beyn_final_s0_relative_change: float
    tighter_margin_rank_consistent: bool | None
    local_refinement_seed_count: int
    resolved_mode_count: int
    kb_numerical: float
    k2_numerical: float
    sigma_numerical: float
    sigma_over_epsilon_squared: float
    scaled_asymptotic_remainder: float
    relative_error_sigma: float
    sigma_min_final: float
    relative_singular_value_final: float
    final_drop_factor: float
    final_relative_mesh_change: float  # NaN for compatibility: no M sweep.
    operator_parity_commutator_residual: float
    unique_embedded_mode_supported: bool
    wall_relative_residual: float
    centerline_relative_residual: float
    odd_parity_relative_residual: float
    first_open_channel_relative_amplitude: float
    boundary_integral_relative_residual: float
    decay_rate_left: float
    decay_rate_right: float
    physical_embedded_mode_verified: bool | None
    mode_exists_supported: bool = False
    odd_uniqueness_supported: bool = False
    full_problem_uniqueness_supported: bool = False
    diagnostic_complete: bool = False
    beyn_threshold_rank: int = -1
    beyn_gap_rank: int = -1
    beyn_final_s1_relative_change: float = math.nan
    beyn_moments_converged: bool = False
    normalized_residual: float = math.nan
    equilibrated_residual: float = math.nan
    equilibrated_sigma_min: float = math.nan
    full_vector_residual: float = math.nan
    parity_leakage_residual: float = math.nan
    assembly_version: str = ASSEMBLY_VERSION


@dataclass
class InternalConvergenceRow:
    epsilon: float
    parameter: str
    value: float
    full_M: int
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
# Geometry and theorem quantities
# ---------------------------------------------------------------------------


def reference_obstacle_geometry(
    t: np.ndarray | float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Inflated obstacle Gamma (epsilon removed).

    X(t) = cos t - beta/2 cos(2t)   is even.
    Y(t) = sin t - beta/2 sin(2t)   is odd.
    """
    t = np.asarray(t)
    beta = float(config.shape_beta)

    X = np.cos(t) - 0.5 * beta * np.cos(2.0 * t)
    Y = np.sin(t) - 0.5 * beta * np.sin(2.0 * t)

    Xp = -np.sin(t) + beta * np.sin(2.0 * t)
    Yp = np.cos(t) - beta * np.cos(2.0 * t)

    Xpp = -np.cos(t) + 2.0 * beta * np.cos(2.0 * t)
    Ypp = -np.sin(t) + 2.0 * beta * np.sin(2.0 * t)
    return X, Y, Xp, Yp, Xpp, Ypp


def obstacle_geometry(
    t: np.ndarray | float,
    epsilon: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Physical x-axis-symmetric obstacle required by Theorem 2.3(iii).

    X(t) = eps [cos t - beta/2 cos(2t)]   (even)
    Y(t) = eps [sin t - beta/2 sin(2t)]   (odd)

    The center is a=0 exactly.
    """
    X, Y, Xp, Yp, Xpp, Ypp = reference_obstacle_geometry(t, config)
    eps = float(epsilon)
    return eps * X, eps * Y, eps * Xp, eps * Yp, eps * Xpp, eps * Ypp


def boundary_nodes(M: int) -> np.ndarray:
    if M % 2 != 0:
        raise ValueError("All full boundary orders M must be even for exact y-parity pairing.")
    return (np.arange(M, dtype=float) + 0.5) * (2.0 * PI / M)


def odd_lift_matrix(M: int) -> np.ndarray:
    """Orthonormal basis Q for vectors satisfying u(2pi-t)=-u(t)."""
    if M % 2 != 0:
        raise ValueError("M must be even for odd-in-y reduction.")
    Q = np.zeros((M, M // 2), dtype=np.complex128)
    for j in range(M // 2):
        partner = M - 1 - j
        Q[j, j] = 1.0 / math.sqrt(2.0)
        Q[partner, j] = -1.0 / math.sqrt(2.0)
    return Q


def reflection_matrix(M: int) -> np.ndarray:
    R = np.zeros((M, M), dtype=np.complex128)
    for j in range(M):
        R[M - 1 - j, j] = 1.0
    return R


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


def asymptotic_prediction(
    epsilon: float,
    mu: float,
    config: Config,
) -> tuple[float, float]:
    sigma = asymptotic_coefficient(mu, config) * float(epsilon) ** 2
    k2 = lambda_2(config) - sigma**2
    if k2 <= lambda_1(config):
        raise ValueError(
            "The leading asymptotic value has left [Lambda_1,Lambda_2). "
            "Use smaller epsilon for the theorem-side validation."
        )
    return config.b * math.sqrt(k2), sigma


def sigma_from_kb(kb: float, config: Config) -> float:
    cutoff = kb_2(config)
    if kb >= cutoff:
        return 0.0
    gap_sq_dimensionless = max((cutoff - kb) * (cutoff + kb), 0.0)
    return math.sqrt(gap_sq_dimensionless) / config.b


def kb_from_sigma(sigma: float, config: Config) -> float:
    sigma = float(sigma)
    cutoff = kb_2(config)
    sigma_b = sigma * config.b
    radicand = cutoff**2 - sigma_b**2
    if radicand <= 0.0:
        return 0.0
    root = math.sqrt(radicand)
    gap = sigma_b**2 / max(cutoff + root, np.finfo(float).tiny)
    return cutoff - gap


def sigma_target_band_limits(config: Config) -> tuple[float, float]:
    """Sigma interval corresponding to k in (sqrt(Lambda_1),sqrt(Lambda_2))."""
    k_left = kb_1(config) + config.beyn_lower_threshold_margin
    k_right = kb_2(config) - config.beyn_min_cutoff_margin
    sigma_min = sigma_from_kb(k_right, config)
    sigma_max = sigma_from_kb(k_left, config)
    return float(sigma_min), float(sigma_max)


def effective_cutoff_margin(epsilon: float, mu: float, config: Config) -> float:
    kb_asym, _ = asymptotic_prediction(epsilon, mu, config)
    gap = max(kb_2(config) - kb_asym, config.beyn_min_cutoff_margin)
    return float(
        min(
            config.beyn_cutoff_margin,
            max(
                config.beyn_min_cutoff_margin,
                config.beyn_cutoff_margin_fraction_of_predicted_gap * gap,
            ),
        )
    )


def geometry_admissibility(
    epsilon: float,
    config: Config,
) -> tuple[bool, str]:
    if abs(config.a) > 100.0 * np.finfo(float).eps:
        return False, "Theorem 2.3(iii) requires a=0 exactly."
    if epsilon <= 0.0:
        return False, "epsilon must be positive."

    t = np.linspace(0.0, 2.0 * PI, 4096, endpoint=False)
    X, Y, *_ = obstacle_geometry(t, epsilon, config)
    if float(np.max(np.abs(Y))) >= config.b:
        return False, "Obstacle intersects/touches a waveguide wall."

    # Conservative bound for the local cylindrical Green-series representation.
    max_dx = 2.0 * float(np.max(np.abs(X)))
    max_image_dy = 2.0 * (config.b + float(np.max(np.abs(Y))))
    max_r = math.hypot(max_dx, max_image_dy)
    green_limit = 0.99 * (4.0 * config.b)
    if max_r > green_limit:
        return False, (
            f"Conservative Green-series radius bound exceeded: {max_r:.6g} > "
            f"{green_limit:.6g}."
        )
    return True, "ok"


def compute_shape_diagnostics(config: Config) -> ShapeDiagnostics:
    t = np.linspace(0.0, 2.0 * PI, 8192, endpoint=False)
    X, Y, Xp, Yp, *_ = reference_obstacle_geometry(t, config)
    speed = np.hypot(Xp, Yp)
    area = 0.5 * np.trapezoid(X * Yp - Y * Xp, t)

    # Symmetry check at the same sample values.
    Xm, Ym, *_ = reference_obstacle_geometry(-t, config)
    x_scale = max(float(np.max(np.abs(X))), 1.0)
    y_scale = max(float(np.max(np.abs(Y))), 1.0)
    x_even_defect = float(np.max(np.abs(Xm - X)) / x_scale)
    y_odd_defect = float(np.max(np.abs(Ym + Y)) / y_scale)

    return ShapeDiagnostics(
        beta=float(config.shape_beta),
        area=float(abs(area)),
        max_abs_X=float(np.max(np.abs(X))),
        max_abs_Y=float(np.max(np.abs(Y))),
        min_reference_speed=float(np.min(speed)),
        x_even_defect=x_even_defect,
        y_odd_defect=y_odd_defect,
    )


# ---------------------------------------------------------------------------
# Independent mu computation via exterior Laplace MFS
# ---------------------------------------------------------------------------


def compute_mu_mfs(
    config: Config,
    order: int,
    *,
    beta_override: float | None = None,
) -> MuMFSRow:
    """Compute the dipole strengths (mu,nu) in paper equation (2.2).

    We solve the inflated-obstacle exterior Neumann problem

        Delta Psi = 0,
        dPsi/dn = n_2,
        grad Psi -> 0,

    with an MFS expansion in logarithmic fundamental solutions whose sources are
    placed inside the obstacle.  The zero-total-source constraint removes the
    logarithmic far-field term.  The remaining dipole moment gives

        Psi ~ -mu*y/r^2 - nu*x/r^2.

    The sign convention is calibrated automatically on the unit circle, where
    the exact value is mu=1, nu=0.
    """
    cfg = config if beta_override is None else replace(config, shape_beta=beta_override)
    N = int(order)
    if N < 8:
        raise ValueError("MFS order must be at least 8.")

    tc = 2.0 * PI * (np.arange(N, dtype=float) + 0.5) / N
    X, Y, Xp, Yp, *_ = reference_obstacle_geometry(tc, cfg)
    speed = np.hypot(Xp, Yp)
    if float(np.min(speed)) <= 1.0e-10:
        raise ValueError("Reference obstacle is not regular enough for MFS.")

    # Inward-looking normal of the paper for a CCW parametrization:
    # n = (-Y', X')/|r'|, hence n_2 = X'/|r'|.
    nx = -Yp / speed
    ny = Xp / speed
    g = ny.copy()

    ts = 2.0 * PI * np.arange(N, dtype=float) / N
    Xs, Ys, *_ = reference_obstacle_geometry(ts, cfg)
    rho = float(cfg.mu_mfs_source_scale)
    sx = rho * Xs
    sy = rho * Ys

    dx = X[:, None] - sx[None, :]
    dy = Y[:, None] - sy[None, :]
    r2 = dx * dx + dy * dy
    if np.any(r2 <= np.finfo(float).tiny):
        raise ValueError("MFS source collided with a collocation point.")

    # G=-(1/2pi) log|x-s|, grad_x G=-(x-s)/(2pi|x-s|^2).
    A = -(dx * nx[:, None] + dy * ny[:, None]) / (2.0 * PI * r2)

    # Enforce sum(c_j)=0 exactly by eliminating the final coefficient.
    B = A[:, :-1] - A[:, -1][:, None]
    d, *_ = np.linalg.lstsq(B, g, rcond=cfg.mu_mfs_rcond)
    c = np.concatenate((d, np.array([-float(np.sum(d))])))

    residual = np.linalg.norm(A @ c - g) / max(np.linalg.norm(g), 1.0e-30)
    total_source = float(abs(np.sum(c)))

    # With sum c=0,
    #   Psi(x) ~ (x sum c*sx + y sum c*sy)/(2pi r^2),
    # so comparison with -nu*x/r^2-mu*y/r^2 gives:
    nu = -float(np.dot(c, sx)) / (2.0 * PI)
    mu = -float(np.dot(c, sy)) / (2.0 * PI)

    # Independent target points, not collocation points, and paper (2.3).
    tt=2*PI*(np.arange(4*N,dtype=float)+.173)/(4*N)
    xx,yy,xxp,yyp,*_=reference_obstacle_geometry(tt,cfg);ss=np.hypot(xxp,yyp)
    nxx=-yyp/ss;nyy=xxp/ss
    ddx=xx[:,None]-sx;ddy=yy[:,None]-sy;rr=ddx*ddx+ddy*ddy
    off=-(ddx*nxx[:,None]+ddy*nyy[:,None])/(2*PI*rr)@c
    offgrid=float(np.linalg.norm(off-nyy)/np.linalg.norm(nyy))
    psi=-(np.log(rr)/(4*PI))@c
    dt=2*PI/len(tt);area=.5*dt*np.sum(xx*yyp-yy*xxp)
    identity=(area+dt*np.sum(nyy*psi*ss))/(2*PI)
    identity_error=float(abs(identity-mu)/max(abs(mu),1e-30))
    return MuMFSRow(
        beta=float(cfg.shape_beta),
        order=N,
        source_scale=rho,
        mu=mu,
        nu=nu,
        boundary_relative_residual=float(residual),
        zero_total_source_error=total_source,
        offgrid_relative_residual=offgrid,identity_relative_error=identity_error,
    )


def run_mu_convergence(config: Config) -> tuple[float, float, list[MuMFSRow], MuMFSRow]:
    # Circle calibration is an independent sign/normalization check.
    circle = compute_mu_mfs(
        config,
        max(config.mu_mfs_orders),
        beta_override=0.0,
    )
    if abs(circle.mu - 1.0) > config.mu_circle_calibration_tolerance:
        raise RuntimeError(
            f"MFS circle calibration failed: mu={circle.mu:.12g}, expected 1."
        )

    rows = [compute_mu_mfs(config, N) for N in config.mu_mfs_orders]
    finest = rows[-1]
    if len(rows) >= 2:
        rel = abs(rows[-1].mu - rows[-2].mu) / max(abs(rows[-1].mu), 1.0e-30)
        if rel > config.mu_convergence_relative_tolerance:
            raise RuntimeError(f"mu MFS not converged: {rel:.3e}")

    for row in rows+[circle]:
        if (not np.all(np.isfinite([row.mu,row.nu,row.offgrid_relative_residual,row.identity_relative_error]))
                or row.offgrid_relative_residual>config.mu_boundary_residual_tolerance
                or row.identity_relative_error>config.mu_identity_relative_tolerance):
            raise RuntimeError("MFS off-grid/identity diagnostics failed")
    if finest.mu <= 0.0:
        raise RuntimeError(f"Computed mu must be positive; got {finest.mu}.")
    return float(finest.mu), float(finest.nu), rows, circle


# ---------------------------------------------------------------------------
# Full-strip Green function and BEM operator
# ---------------------------------------------------------------------------


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

    period = 2*d
    def wrap(value):
        return (value + period/2) % period - period/2
    def green(x, y, xi, eta):
        return (lattice.greens_periodic(x-xi, wrap(y-eta), coefficients, k, period)
                - lattice.greens_periodic(x-xi, wrap(y+eta+2*b), coefficients, k, period))
    def green_regularized(x, y, xi, eta):
        return (lattice.greens_periodic_reg(x-xi, wrap(y-eta), coefficients, k, period)
                - lattice.greens_periodic(x-xi, wrap(y+eta+2*b), coefficients, k, period))

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
def _kernel_geometry(t: float, epsilon: float, config: Config) -> tuple[float, ...]:
    return tuple(float(v) for v in obstacle_geometry(t, epsilon, config))


def weighted_normal_kernel(
    psi: float,
    theta: float,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
    G_regularized: Callable[..., complex],
) -> complex:
    x, y, *_ = _kernel_geometry(psi, epsilon, config)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = _kernel_geometry(theta, epsilon, config)

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

    # Shifted midpoint nodes are distinct unless psi==theta from the same grid.
    periodic_distance = abs(math.atan2(math.sin(psi - theta), math.cos(psi - theta)))
    if periodic_distance > 1.0e-12:
        G_xi, G_eta = source_derivatives(G, x, y, xi, eta, h)
        return xi_p * G_eta - eta_p * G_xi

    G_xi_reg, G_eta_reg = source_derivatives(G_regularized, x, y, xi, eta, h)
    geometric_term = (xi_pp * eta_p - eta_pp * xi_p) / (4.0 * PI * w**2)
    regularized_term = xi_p * G_eta_reg - eta_p * G_xi_reg
    return geometric_term + regularized_term


def assemble_full_matrix(
    kb: complex,
    epsilon: float,
    M: int,
    config: Config,
) -> np.ndarray:
    period = 4*config.b
    k = complex(kb)/config.b
    coeffs = lattice.lattice_sums(period,k,beta=0.,M=config.lattice_terms,Lh=config.harmonic_order)
    X,Y,Dx,Dy,Iy,xp,yp,curvature,plans = _assembly_geometry(epsilon,M,config)
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
def _assembly_geometry(epsilon,M,config):
    t=boundary_nodes(M);x,y,xp,yp,xpp,ypp=obstacle_geometry(t,epsilon,config)
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


def assemble_odd_matrix(
    kb: complex,
    epsilon: float,
    M: int,
    config: Config,
) -> np.ndarray:
    A = assemble_full_matrix(kb, epsilon, M, config)
    Q = odd_lift_matrix(M)
    return Q.conj().T @ A @ Q


def parity_commutator_residual(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> float:
    A = assemble_full_matrix(complex(kb, 0.0), epsilon, M, config)
    R = reflection_matrix(M)
    numerator = np.linalg.norm(A @ R - R @ A)
    denominator = max(np.linalg.norm(A), 1.0e-30)
    return float(numerator / denominator)


def odd_singular_metrics(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float]:
    A = assemble_odd_matrix(complex(kb, 0.0), epsilon, M, config)
    s = np.linalg.svd(A, compute_uv=False)
    smax = float(s[0])
    smin = float(s[-1])
    return smin, smax, smin / max(smax, np.finfo(float).tiny)


def odd_singular_pair(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, np.ndarray, np.ndarray]:
    A = assemble_odd_matrix(complex(kb, 0.0), epsilon, M, config)
    _, s, Vh = np.linalg.svd(A, full_matrices=False)
    vred = Vh.conj().T[:, -1]
    vred /= max(float(np.linalg.norm(vred)), 1.0e-30)
    Q = odd_lift_matrix(M)
    vfull = Q @ vred
    pivot = int(np.argmax(np.abs(vfull)))
    if abs(vfull[pivot]) > 0.0:
        phase = np.exp(-1j * np.angle(vfull[pivot]))
        vred *= phase
        vfull *= phase
    smax = float(s[0])
    smin = float(s[-1])
    return smin, smax, smin / max(smax, np.finfo(float).tiny), vred, vfull


# ---------------------------------------------------------------------------
# Beyn global discovery on the odd sector
# ---------------------------------------------------------------------------


def contour_geometry(
    epsilon: float,
    mu: float,
    config: Config,
    *,
    cutoff_margin_override: float | None = None,
) -> tuple[float, float, float, float, float, float]:
    left = kb_1(config) + config.beyn_lower_threshold_margin
    margin = (
        effective_cutoff_margin(epsilon, mu, config)
        if cutoff_margin_override is None
        else float(cutoff_margin_override)
    )
    right = kb_2(config) - margin
    if not (kb_1(config) < left < right < kb_2(config)):
        raise ValueError("Invalid Theorem 2.3(iii) Beyn contour.")
    center = 0.5 * (left + right)
    rx = 0.5 * (right - left)
    ry = config.beyn_imag_half_height
    return left, right, center, rx, ry, margin


def ellipse_point(
    theta: float,
    epsilon: float,
    mu: float,
    config: Config,
    *,
    cutoff_margin_override: float | None = None,
) -> tuple[complex, complex]:
    _, _, center, rx, ry, _ = contour_geometry(
        epsilon, mu, config, cutoff_margin_override=cutoff_margin_override
    )
    z = center + rx * math.cos(theta) + 1j * ry * math.sin(theta)
    dz = -rx * math.sin(theta) + 1j * ry * math.cos(theta)
    return complex(z), complex(dz)


def point_inside_ellipse(
    z: complex,
    epsilon: float,
    mu: float,
    config: Config,
    *,
    cutoff_margin_override: float | None = None,
    tolerance: float = 1.0e-9,
) -> bool:
    _, _, center, rx, ry, _ = contour_geometry(
        epsilon, mu, config, cutoff_margin_override=cutoff_margin_override
    )
    q = ((z.real - center) / rx) ** 2 + (z.imag / ry) ** 2
    return bool(q <= 1.0 + tolerance)


def probing_matrix(reduced_dim: int, config: Config) -> np.ndarray:
    if config.beyn_probe_dim >= reduced_dim:
        raise ValueError(
            "beyn_probe_dim must be smaller than the odd reduced matrix dimension."
        )
    rng = np.random.default_rng(config.beyn_random_seed)
    V = rng.standard_normal((reduced_dim, config.beyn_probe_dim)) + 1j * rng.standard_normal(
        (reduced_dim, config.beyn_probe_dim)
    )
    Q, _ = np.linalg.qr(V)
    return Q[:, : config.beyn_probe_dim]


def estimate_beyn_rank(
    singular_values: np.ndarray,
    config: Config,
) -> tuple[int, int, int, float]:
    s = np.asarray(singular_values, dtype=float)
    if len(s) == 0 or not np.isfinite(s[0]) or s[0] <= config.beyn_empty_s0_tolerance:
        return 0, 0, 0, math.nan

    threshold = max(
        config.beyn_rank_absolute_tolerance,
        config.beyn_rank_relative_tolerance * float(s[0]),
    )
    threshold_rank = int(np.sum(s > threshold))
    if threshold_rank == 0:
        return 0, 0, 0, math.nan

    threshold_rank = min(threshold_rank, len(s))
    n_gaps = min(threshold_rank, len(s) - 1)
    if n_gaps <= 0:
        return threshold_rank, threshold_rank, threshold_rank, math.inf

    floor = max(config.beyn_rank_absolute_tolerance, np.finfo(float).tiny)
    ratios = s[:n_gaps] / np.maximum(s[1 : n_gaps + 1], floor)
    idx = int(np.argmax(ratios))
    gap_rank = idx + 1
    selected_gap = float(ratios[idx])
    rank = gap_rank if selected_gap >= config.beyn_rank_gap_threshold else threshold_rank
    return int(rank), int(threshold_rank), int(gap_rank), selected_gap


def cluster_real_candidates(values: list[complex], config: Config) -> list[complex]:
    if not values:
        return []
    values = sorted(values, key=lambda z: z.real)
    clusters: list[list[complex]] = [[values[0]]]
    for z in values[1:]:
        if abs(z.real - clusters[-1][-1].real) <= config.beyn_cluster_tolerance:
            clusters[-1].append(z)
        else:
            clusters.append([z])
    return [sum(c) / len(c) for c in clusters]


def contour_quadrature(points):
    edges=np.unique(np.r_[0.,np.geomspace(1e-12,PI/2,16),np.linspace(.01,PI/2,17)])
    order=max(2,math.ceil(points/(4*(len(edges)-1))))
    x,w=leggauss(order);nodes=[];weights=[]
    for quadrant in range(4):
        for lo,hi in zip(edges[:-1],edges[1:]):
            t=lo+(hi-lo)*(x+1)/2
            nodes.extend(quadrant*PI/2+(t if quadrant%2==0 else PI/2-t))
            weights.extend(w*(hi-lo)/2)
    return np.asarray(nodes),np.asarray(weights)


def beyn_discover(
    epsilon: float,
    quadrature_points: int,
    mu: float,
    config: Config,
    *,
    cutoff_margin_override: float | None = None,
) -> tuple[np.ndarray, np.ndarray, BeynDiagnostics, list[BeynEigenvalueRow], np.ndarray]:
    M = int(config.bem_M)
    reduced_dim = M // 2
    Nq = int(quadrature_points)
    V = probing_matrix(reduced_dim, config)

    S0 = np.zeros((reduced_dim, config.beyn_probe_dim), dtype=np.complex128)
    S1 = np.zeros_like(S0)
    max_residual = 0.0

    angles,weights=contour_quadrature(Nq)
    for theta,angular_weight in zip(angles,weights):
        z, dz = ellipse_point(
            theta,
            epsilon,
            mu,
            config,
            cutoff_margin_override=cutoff_margin_override,
        )
        A = assemble_odd_matrix(z, epsilon, M, config)
        if not np.all(np.isfinite(A)):raise ValueError("Nonfinite contour matrix")
        X = np.linalg.solve(A, V)
        rel_res = np.linalg.norm(A @ X - V) / max(np.linalg.norm(V), 1.0e-30)
        if not np.isfinite(rel_res):raise ValueError("Nonfinite contour solve")
        max_residual = max(max_residual, float(rel_res))
        weight = dz * angular_weight / (2j * PI)
        S0 += weight * X
        S1 += weight * z * X

    U, s, Vh = np.linalg.svd(S0, full_matrices=False)
    rank, threshold_rank, gap_rank, selected_gap = estimate_beyn_rank(s, config)
    if rank == 0:
        eigenvalues = np.array([], dtype=np.complex128)
    else:
        Ur = U[:, :rank]
        Wr = Vh[:rank, :].conj().T
        sr = s[:rank]
        B = Ur.conj().T @ S1 @ Wr @ np.diag(1.0 / sr)
        eigenvalues = np.linalg.eigvals(B)

    left, right, *_ = contour_geometry(
        epsilon, mu, config, cutoff_margin_override=cutoff_margin_override
    )
    rows: list[BeynEigenvalueRow] = []
    accepted = 0
    near_count = 0
    for idx, eig0 in enumerate(eigenvalues):
        eig = complex(eig0)
        inside = point_inside_ellipse(
            eig,
            epsilon,
            mu,
            config,
            cutoff_margin_override=cutoff_margin_override,
        )
        inside_band = bool(kb_1(config) < eig.real < kb_2(config))
        inside_span = bool(left < eig.real < right)
        near_real = bool(abs(eig.imag) <= config.beyn_real_axis_tolerance)
        overrun = max(left - eig.real, eig.real - right, 0.0)
        accepted_here = bool(inside and inside_band and inside_span and near_real)
        near_seed = bool(
            not accepted_here
            and inside_band
            and near_real
            and overrun <= config.beyn_near_contour_real_tolerance
        )
        accepted += int(accepted_here)
        near_count += int(near_seed)
        rows.append(
            BeynEigenvalueRow(
                epsilon=float(epsilon),
                quadrature_points=Nq,
                raw_index=idx,
                real_part=float(eig.real),
                imag_part=float(eig.imag),
                inside_contour=inside,
                inside_target_band=inside_band,
                near_real=near_real,
                real_span_overrun=float(overrun),
                near_contour_seed=near_seed,
                accepted_for_refinement=accepted_here,
            )
        )

    diag = BeynDiagnostics(
        epsilon=float(epsilon),
        full_M=M,
        reduced_dimension=reduced_dim,
        quadrature_points=Nq,
        estimated_rank=rank,
        threshold_rank=threshold_rank,
        gap_rank=gap_rank,
        selected_gap_ratio=selected_gap,
        leading_s0_singular_value=float(s[0]) if len(s) else math.nan,
        max_linear_solve_relative_residual=float(max_residual),
        raw_eigenvalues=len(eigenvalues),
        accepted_candidates=accepted,
        near_contour_seeds=near_count,
    )
    diag._s0=S0;diag._s1=S1
    return eigenvalues, s, diag, rows, S0


def candidate_values(
    eigenvalues: np.ndarray,
    rows: list[BeynEigenvalueRow],
    config: Config,
) -> tuple[list[complex], list[complex]]:
    strict = [
        complex(eigenvalues[r.raw_index]) for r in rows if r.accepted_for_refinement
    ]
    near = [
        complex(eigenvalues[r.raw_index]) for r in rows if r.near_contour_seed
    ]
    return cluster_real_candidates(strict, config), cluster_real_candidates(near, config)


def candidate_set_shift(previous: list[complex], current: list[complex]) -> float:
    if len(previous) != len(current) or not current:
        return math.nan
    p = sorted(previous, key=lambda z: z.real)
    c = sorted(current, key=lambda z: z.real)
    return float(max(abs(a.real - b.real) for a, b in zip(p, c, strict=True)))


def aitken_delta_squared(values: list[float]) -> float:
    if len(values) < 3:
        return math.nan
    x0, x1, x2 = values[-3:]
    den = x2 - 2.0 * x1 + x0
    scale = max(abs(x0), abs(x1), abs(x2), 1.0)
    if abs(den) <= 100.0 * np.finfo(float).eps * scale:
        return math.nan
    return float(x0 - (x1 - x0) ** 2 / den)


def run_beyn_convergence_study(
    epsilon: float,
    mu: float,
    config: Config,
) -> tuple[
    np.ndarray,
    np.ndarray,
    BeynDiagnostics,
    list[BeynEigenvalueRow],
    list[DiscoverySeed],
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
    float,
]:
    levels = list(config.beyn_quadrature_levels)
    diagnostics: list[BeynDiagnostics] = []
    all_rows: list[BeynEigenvalueRow] = []
    convergence: list[BeynConvergenceRow] = []
    previous_S0: np.ndarray | None = None
    previous_S1: np.ndarray | None = None
    previous_candidates: list[complex] = []
    raw_real_history: list[float] = []

    final_eigs = np.array([], dtype=np.complex128)
    final_s = np.array([], dtype=float)
    final_diag: BeynDiagnostics | None = None
    final_rows: list[BeynEigenvalueRow] = []
    final_strict: list[complex] = []
    final_near: list[complex] = []

    def run_level(Nq: int) -> None:
        nonlocal previous_S0, previous_S1, previous_candidates
        nonlocal final_eigs, final_s, final_diag, final_rows, final_strict, final_near

        eigs, s, diag, rows, S0 = beyn_discover(epsilon, Nq, mu, config)
        strict, near = candidate_values(eigs, rows, config)
        primary_pool = strict if strict else near
        primary = primary_pool[0] if primary_pool else None
        if primary is None and len(eigs):
            # Diagnostic raw value nearest the upper threshold.
            raw = min((complex(z) for z in eigs), key=lambda z: abs(z.real - kb_2(config)))
        else:
            raw = primary
        if raw is not None and np.isfinite(raw.real):
            raw_real_history.append(float(raw.real))

        s0_change = math.nan
        if previous_S0 is not None:
            s0_change = float(
                np.linalg.norm(S0 - previous_S0) / max(np.linalg.norm(S0), 1.0e-30)
            )
        shift = candidate_set_shift(previous_candidates, strict)
        rank_stable = bool(
            diagnostics and diagnostics[-1].estimated_rank == diag.estimated_rank
        )
        s1_change=math.nan if previous_S1 is None else float(np.linalg.norm(diag._s1-previous_S1)/max(np.linalg.norm(diag._s1),1e-30))
        diag.s0_relative_change=s0_change;diag.s1_relative_change=s1_change
        diag.moments_converged=bool(np.isfinite(s0_change) and np.isfinite(s1_change)
                                    and max(s0_change,s1_change)<=config.beyn_moment_convergence_tolerance)
        diag.rank_stable=rank_stable
        diag.candidate_stable=bool(strict and np.isfinite(shift) and shift<=config.beyn_candidate_convergence_tolerance)
        convergence.append(
            BeynConvergenceRow(
                epsilon=float(epsilon),
                quadrature_points=int(Nq),
                estimated_rank=int(diag.estimated_rank),
                s0_sv1=float(s[0]) if len(s) > 0 else math.nan,
                s0_sv2=float(s[1]) if len(s) > 1 else math.nan,
                s0_sv3=float(s[2]) if len(s) > 2 else math.nan,
                s0_relative_change_from_previous=s0_change,
                accepted_candidates=len(strict),
                primary_candidate_real=(float(primary.real) if primary is not None else math.nan),
                primary_candidate_imag=(float(primary.imag) if primary is not None else math.nan),
                max_candidate_shift_from_previous=shift,
                rank_stable_from_previous=rank_stable,
                s1_relative_change_from_previous=s1_change,
                moments_converged=diag.moments_converged,
            )
        )
        diagnostics.append(diag)
        all_rows.extend(rows)
        previous_S0 = S0
        previous_S1 = diag._s1
        previous_candidates = strict
        final_eigs, final_s, final_diag, final_rows = eigs, s, diag, rows
        final_strict, final_near = strict, near

    for Nq in levels:
        run_level(int(Nq))

    # Escalate only when the discovery is not yet clean enough.
    for Nq in config.beyn_adaptive_quadrature_levels:
        rank_stable = bool(
            len(convergence) >= 2
            and convergence[-1].estimated_rank == convergence[-2].estimated_rank
        )
        if (final_strict and rank_stable and final_diag.moments_converged
                and final_diag.candidate_stable and final_diag.threshold_rank==final_diag.gap_rank==final_diag.estimated_rank
                and final_diag.max_linear_solve_relative_residual<=config.beyn_linear_solve_tolerance):
            break
        run_level(int(Nq))

    assert final_diag is not None

    # Preserve strict candidates from every level; a candidate disappearing
    # under refinement is not silently erased from the diagnostic.
    accepted_history=[complex(row.real_part,row.imag_part) for row in all_rows if row.accepted_for_refinement]
    all_strict=cluster_real_candidates(accepted_history,config)
    seeds: list[DiscoverySeed] = [
        DiscoverySeed(float(z.real),float(z.imag),"strict-beyn") for z in all_strict
    ]
    if not seeds:
        seeds.extend(
            DiscoverySeed(float(z.real), float(z.imag), "near-contour") for z in final_near
        )

    aitken = math.nan
    if not seeds and final_diag.estimated_rank == 1 and len(raw_real_history) >= 3:
        aitken = aitken_delta_squared(raw_real_history)
        if np.isfinite(aitken) and kb_1(config) < aitken < kb_2(config):
            seeds.append(DiscoverySeed(float(aitken), 0.0, "aitken"))

    return (
        final_eigs,
        final_s,
        final_diag,
        final_rows,
        seeds,
        diagnostics,
        all_rows,
        convergence,
        float(aitken),
    )


# ---------------------------------------------------------------------------
# Local sigma-SVD refinement
# ---------------------------------------------------------------------------


def absolute_sv_from_sigma(
    sigma: float,
    epsilon: float,
    M: int,
    config: Config,
) -> float:
    kb = kb_from_sigma(sigma, config)
    return odd_singular_metrics(kb, epsilon, M, config)[0]


def refine_candidate_for_M(
    epsilon: float,
    M: int,
    sigma_left: float,
    sigma_right: float,
    config: Config,
) -> tuple[float, float, float, float, float, float, bool]:
    if not (0.0 < sigma_left < sigma_right):
        raise ValueError("Invalid sigma bracket.")

    # Coarse local scan prevents a bounded minimizer from locking to an edge.
    grid = np.linspace(sigma_left, sigma_right, 41)
    vals = np.array(
        [absolute_sv_from_sigma(float(s), epsilon, M, config) for s in grid],
        dtype=float,
    )
    local = [
        i for i in range(1, len(grid) - 1)
        if vals[i] < vals[i - 1] and vals[i] < vals[i + 1]
    ]
    if local:
        idx = min(local, key=lambda i: vals[i])
    else:
        idx = int(np.argmin(vals))
    lo = float(grid[max(0, idx - 1)])
    hi = float(grid[min(len(grid) - 1, idx + 1)])
    if not lo < hi:
        lo, hi = sigma_left, sigma_right

    result = cast(Any, minimize_scalar(
        lambda s: absolute_sv_from_sigma(float(s), epsilon, M, config),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": config.minimizer_sigma_xatol, "maxiter": 200},
    ))
    if not result.success:raise RuntimeError(f"Sigma minimization failed: {result.message}")
    if not np.all(np.isfinite(vals)):raise ValueError("Nonfinite local scan")
    sigma_star = float(result.x)
    # scipy bounded uses sqrt(eps)*abs(x) as part of its stopping rule.
    # Center a second minimization to reach the requested absolute sigma
    # resolution without changing the singular-value objective or tolerances.
    width=hi-lo
    polish=cast(Any,minimize_scalar(
        lambda t:absolute_sv_from_sigma(sigma_star+width*float(t),epsilon,M,config),
        bounds=((lo-sigma_star)/width,(hi-sigma_star)/width),method="bounded",
        options={"xatol":config.minimizer_sigma_xatol/width,"maxiter":200}))
    if not polish.success:raise RuntimeError(f"Centered sigma minimization failed: {polish.message}")
    polished_sigma=sigma_star+width*float(polish.x)
    if absolute_sv_from_sigma(polished_sigma,epsilon,M,config)<absolute_sv_from_sigma(sigma_star,epsilon,M,config):
        sigma_star=polished_sigma
    kb_star = kb_from_sigma(sigma_star, config)
    smin, smax, rel = odd_singular_metrics(kb_star, epsilon, M, config)

    left_val = absolute_sv_from_sigma(sigma_left, epsilon, M, config)
    right_val = absolute_sv_from_sigma(sigma_right, epsilon, M, config)
    drop = min(left_val, right_val) / max(smin, np.finfo(float).tiny)

    edge_margin = 1.0e-3 * (sigma_right - sigma_left)
    interior = bool(
        sigma_star > sigma_left + edge_margin
        and sigma_star < sigma_right - edge_margin
    )
    if not np.all(np.isfinite([kb_star,sigma_star,smin,smax,rel,drop,left_val,right_val])):
        raise ValueError("Nonfinite refinement diagnostics")
    return kb_star, sigma_star, smin, smax, rel, float(drop), interior


def sigma_bracket_from_seed(seed: DiscoverySeed, config: Config) -> tuple[float, float]:
    sigma_seed = sigma_from_kb(seed.real, config)
    global_left, global_right = sigma_target_band_limits(config)
    factor = max(config.local_sigma_seed_bracket_factor, 1.01)
    left = max(global_left, sigma_seed / factor)
    right = min(global_right, sigma_seed * factor)
    if not left < sigma_seed < right:
        width = max(0.1 * sigma_seed, 1.0e-6)
        left = max(global_left, sigma_seed - width)
        right = min(global_right, sigma_seed + width)
    if not left < right:
        raise ValueError("Could not create a sigma bracket from Beyn seed.")
    return float(left), float(right)


def global_sigma_fallback_seed(
    epsilon: float,
    config: Config,
) -> tuple[DiscoverySeed, tuple[float, float]] | None:
    """Theory-independent fallback over the target odd spectral band."""
    smin, smax = sigma_target_band_limits(config)
    n = max(int(config.local_sigma_scan_points), 64)
    log_grid = np.geomspace(smin, smax, n)
    linear_grid = np.linspace(max(smin, 1.0e-3 * smax), smax, max(n // 2, 32))
    sigmas = np.unique(np.concatenate((log_grid, linear_grid)))
    M = int(config.bem_M)
    values = np.array(
        [absolute_sv_from_sigma(float(s), epsilon, M, config) for s in sigmas],
        dtype=float,
    )
    minima = [
        i for i in range(1, len(sigmas) - 1)
        if values[i] < values[i - 1] and values[i] < values[i + 1]
    ]
    if not minima:
        return None
    idx = min(minima, key=lambda i: values[i])
    seed_sigma = float(sigmas[idx])
    seed = DiscoverySeed(kb_from_sigma(seed_sigma, config), 0.0, "global-sigma-scan")
    return seed, (float(sigmas[idx - 1]), float(sigmas[idx + 1]))


def certify_odd_mode(kb,epsilon,M,interior,drop,config):
    full=assemble_full_matrix(complex(kb),epsilon,M,config);Q=odd_lift_matrix(M)
    if not np.all(np.isfinite(full)):raise ValueError("Nonfinite mode matrix")
    A=Q.conj().T@full@Q;_,sv,vh=np.linalg.svd(A);v=vh.conj().T[:,-1]
    rows=np.linalg.norm(A,axis=1)
    if sv[0]<=0 or np.any(rows<=np.finfo(float).tiny):raise ValueError("Degenerate mode matrix")
    B=A/rows[:,None];scaled=np.linalg.svd(B,compute_uv=False)
    residual=float(np.linalg.norm(A@v)/(sv[0]*np.linalg.norm(v)))
    eq=float(np.linalg.norm(B@v)/(scaled[0]*np.linalg.norm(v)))
    vf=Q@v;fullnorm=np.linalg.norm(full,2)
    fullres=float(np.linalg.norm(full@vf)/(fullnorm*np.linalg.norm(vf)))
    leakage=float(np.linalg.norm((np.eye(M)-Q@Q.conj().T)@full@vf)/(fullnorm*np.linalg.norm(vf)))
    values=[residual,eq,scaled[-1],fullres,leakage,drop]
    supported=bool(np.all(np.isfinite(values)) and interior and drop>=config.minimum_drop_factor
                   and sv[-1]/sv[0]<=config.relative_near_singular_tolerance
                   and residual<=config.singular_vector_residual_tolerance
                   and eq<=config.equilibrated_residual_tolerance
                   and scaled[-1]<=config.near_singular_tolerance
                   and fullres<=config.full_vector_residual_tolerance
                   and leakage<=config.parity_operator_tolerance)
    return dict(supported=supported,normalized_residual=residual,equilibrated_residual=eq,
                equilibrated_sigma_min=float(scaled[-1]),full_vector_residual=fullres,parity_leakage_residual=leakage)


def run_candidate_refinement(
    epsilon: float,
    candidate_index: int,
    seed: DiscoverySeed,
    sigma_bracket: tuple[float, float],
    config: Config,
) -> tuple[list[ModeRefinementRow], ModeResult]:
    rows: list[ModeRefinementRow] = []
    left, right = sigma_bracket

    for M in (config.bem_M,):
        kb, sigma, smin, smax, rel, drop, interior = refine_candidate_for_M(
            epsilon, int(M), left, right, config
        )
        change = math.nan  # Compatibility field; not applicable at fixed M.
        rows.append(
            ModeRefinementRow(
                epsilon=float(epsilon),
                candidate_index=int(candidate_index),
                seed_source=seed.source,
                full_M=int(M),
                reduced_dimension=int(M) // 2,
                kb=float(kb),
                sigma_bem=float(sigma),
                sigma_min=float(smin),
                sigma_max=float(smax),
                relative_singular_value=float(rel),
                drop_factor=float(drop),
                minimum_is_interior=bool(interior),
                relative_sigma_change_from_previous=float(change),
            )
        )

    final = rows[-1]
    final_change = math.nan  # Not applicable: certification uses one fixed M.
    certificate=certify_odd_mode(final.kb,epsilon,config.bem_M,final.minimum_is_interior,final.drop_factor,config)
    resolved=certificate['supported']
    mode = ModeResult(
        epsilon=float(epsilon),
        candidate_index=int(candidate_index),
        seed_source=seed.source,
        kb_numerical=final.kb,
        sigma_numerical=final.sigma_bem,
        sigma_min_final=final.sigma_min,
        sigma_max_final=final.sigma_max,
        relative_singular_value_final=final.relative_singular_value,
        final_drop_factor=final.drop_factor,
        final_relative_mesh_change=float(final_change),
        resolved=resolved,
        normalized_residual=certificate['normalized_residual'],
        equilibrated_residual=certificate['equilibrated_residual'],
        equilibrated_sigma_min=certificate['equilibrated_sigma_min'],
        full_vector_residual=certificate['full_vector_residual'],
        parity_leakage_residual=certificate['parity_leakage_residual'],
    )
    return rows, mode


# ---------------------------------------------------------------------------
# Physical diagnostics for the embedded odd mode
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
    config: Config,
    G: Callable[..., complex],
) -> complex:
    xi, eta, xi_p, eta_p, *_ = obstacle_geometry(theta, epsilon, config)
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
    config: Config,
    G: Callable[..., complex],
) -> np.ndarray:
    dtheta = 2.0 * PI / len(theta)
    out = np.empty(len(points), dtype=np.complex128)
    for i, (x, y) in enumerate(points):
        kernel = np.array(
            [
                weighted_kernel_at_field_point(x, y, float(t), epsilon, config, G)
                for t in theta
            ],
            dtype=np.complex128,
        )
        out[i] = dtheta * np.dot(kernel, boundary_vector)
    return out


def periodic_fourier_interpolate(
    theta_nodes: np.ndarray,
    values: np.ndarray,
    theta_targets: np.ndarray,
) -> np.ndarray:
    M = len(theta_nodes)
    theta0 = float(theta_nodes[0])
    coeffs = np.fft.fft(values) / M
    frequencies = np.fft.fftfreq(M, d=1.0 / M)
    phase = np.exp(1j * np.outer(np.asarray(theta_targets) - theta0, frequencies))
    return phase @ coeffs


def offgrid_boundary_integral_residual(
    epsilon: float,
    kb: float,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    config: Config,
) -> float:
    sample_count = max(config.physical_boundary_residual_samples, 2 * len(theta))
    targets = (
        (np.arange(sample_count, dtype=float) + 0.371)
        * (2.0 * PI / sample_count)
    ) % (2.0 * PI)
    u_targets = periodic_fourier_interpolate(theta, boundary_vector, targets)
    G, Greg = make_green_functions(complex(kb, 0.0), config)
    dtheta = 2.0 * PI / len(theta)

    residuals = np.empty(sample_count, dtype=np.complex128)
    rhs_values = np.empty(sample_count, dtype=np.complex128)
    for i, psi in enumerate(targets):
        kernel = np.array(
            [
                weighted_normal_kernel(
                    float(psi), float(src), epsilon, config, G, Greg
                )
                for src in theta
            ],
            dtype=np.complex128,
        )
        rhs = dtheta * np.dot(kernel, boundary_vector)
        rhs_values[i] = rhs
        residuals[i] = 0.5 * u_targets[i] - rhs

    scale = max(
        float(np.max(np.abs(0.5 * u_targets))),
        float(np.max(np.abs(rhs_values))),
        float(np.max(np.abs(boundary_vector))),
        1.0e-30,
    )
    return float(np.max(np.abs(residuals)) / scale)


def physical_mode_diagnostics(
    epsilon: float,
    kb: float,
    M: int,
    config: Config,
) -> PhysicalDiagnostics:
    _, _, _, _, boundary_vector = odd_singular_pair(kb, epsilon, M, config)
    theta = boundary_nodes(M)
    G = make_extended_field_green(kb, config)
    b = config.b

    nodes, weights = leggauss(config.physical_quadrature_points)
    y_values = b * nodes
    y_weights = b * weights
    distances = np.asarray(config.physical_x_over_b, dtype=float) * b

    left_norms: list[float] = []
    right_norms: list[float] = []
    open_ratios: list[float] = []
    cross_peak = 0.0

    phi1 = np.cos(PI * y_values / (2.0 * b))
    phi1_norm = math.sqrt(float(np.sum(y_weights * phi1**2)))

    for x in distances:
        for sign, target in [(-1.0, left_norms), (1.0, right_norms)]:
            pts = [(float(sign * x), float(y)) for y in y_values]
            field = reconstruct_field_at_points(
                pts, boundary_vector, theta, epsilon, config, G
            )
            norm = math.sqrt(float(np.sum(y_weights * np.abs(field) ** 2)))
            target.append(norm)
            cross_peak = max(cross_peak, float(np.max(np.abs(field))))
            overlap = abs(np.sum(y_weights * np.conjugate(phi1) * field))
            open_ratio = float(overlap / max(phi1_norm * norm, 1.0e-30))
            open_ratios.append(open_ratio)

    left_arr = np.asarray(left_norms)
    right_arr = np.asarray(right_norms)
    monotone_left = bool(np.all(np.diff(left_arr) <= 1.0e-10 * max(left_arr[0], 1.0)))
    monotone_right = bool(np.all(np.diff(right_arr) <= 1.0e-10 * max(right_arr[0], 1.0)))

    def fit_decay(norms: np.ndarray) -> float:
        slope, _ = np.polyfit(distances, np.log(np.maximum(norms, np.finfo(float).tiny)), 1)
        return float(-slope)

    decay_left = fit_decay(left_arr)
    decay_right = fit_decay(right_arr)
    expected = sigma_from_kb(kb, config)
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
        wall_points, boundary_vector, theta, epsilon, config, G
    )

    # Centerline values: exact odd symmetry demands u(x,0)=0.
    centerline_points = [(float(x), 0.0) for x in (-x0, -0.5*x0, 0.5*x0, x0)]
    centerline_values = reconstruct_field_at_points(
        centerline_points, boundary_vector, theta, epsilon, config, G
    )

    # Direct parity test away from the obstacle and away from y=0.
    parity_points: list[tuple[float, float]] = []
    parity_pairs: list[tuple[int, int]] = []
    for x in (-x0, x0):
        for frac in (0.20, 0.45, 0.70):
            y = frac * b
            i = len(parity_points)
            parity_points.extend([(float(x), y), (float(x), -y)])
            parity_pairs.append((i, i + 1))
    parity_values = reconstruct_field_at_points(
        parity_points, boundary_vector, theta, epsilon, config, G
    )

    amplitude_scale = max(
        float(np.max(np.abs(boundary_vector))),
        cross_peak,
        1.0e-30,
    )
    wall_rel = float(np.max(np.abs(wall_values)) / amplitude_scale)
    centerline_rel = float(np.max(np.abs(centerline_values)) / amplitude_scale)
    parity_rel = 0.0
    for i, j in parity_pairs:
        parity_rel = max(
            parity_rel,
            float(abs(parity_values[i] + parity_values[j]) / amplitude_scale),
        )

    open_channel_rel = float(max(open_ratios) if open_ratios else math.nan)
    bie_rel = offgrid_boundary_integral_residual(
        epsilon, kb, boundary_vector, theta, config
    )

    walls_ok = wall_rel <= config.physical_wall_relative_tolerance
    center_ok = centerline_rel <= config.physical_centerline_relative_tolerance
    parity_ok = parity_rel <= config.physical_parity_relative_tolerance
    open_ok = open_channel_rel <= config.physical_open_channel_relative_tolerance
    bie_ok = bie_rel <= config.physical_boundary_residual_tolerance
    decay_ok = bool(
        monotone_left
        and monotone_right
        and decay_error_left <= config.physical_decay_relative_tolerance
        and decay_error_right <= config.physical_decay_relative_tolerance
    )

    return PhysicalDiagnostics(
        epsilon=float(epsilon),
        kb=float(kb),
        full_M=int(M),
        wall_relative_residual=wall_rel,
        centerline_relative_residual=centerline_rel,
        odd_parity_relative_residual=float(parity_rel),
        first_open_channel_relative_amplitude=open_channel_rel,
        boundary_integral_relative_residual=bie_rel,
        decay_rate_left=decay_left,
        decay_rate_right=decay_right,
        expected_decay_rate=expected,
        decay_relative_error_left=float(decay_error_left),
        decay_relative_error_right=float(decay_error_right),
        monotone_decay_left=monotone_left,
        monotone_decay_right=monotone_right,
        walls_verified=walls_ok,
        centerline_verified=center_ok,
        parity_verified=parity_ok,
        open_channel_suppressed=open_ok,
        boundary_integral_verified=bie_ok,
        decay_verified=decay_ok,
    )


# ---------------------------------------------------------------------------
# One-epsilon theorem validation
# ---------------------------------------------------------------------------


def validate_epsilon(
    epsilon: float,
    mu: float,
    nu: float,
    config: Config,
) -> tuple[
    ValidationResult,
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
    list[ModeRefinementRow],
    list[ModeResult],
    PhysicalDiagnostics | None,
]:
    ok, reason = geometry_admissibility(epsilon, config)
    if not ok:
        raise RuntimeError(f"Invalid geometry at epsilon={epsilon}: {reason}")

    kb_asym, sigma_asym = asymptotic_prediction(epsilon, mu, config)
    coefficient = asymptotic_coefficient(mu, config)
    cutoff_margin = effective_cutoff_margin(epsilon, mu, config)

    (
        final_eigs,
        final_s,
        final_diag,
        final_rows,
        seeds,
        diagnostics,
        all_eigen_rows,
        convergence_rows,
        aitken,
    ) = run_beyn_convergence_study(epsilon, mu, config)

    rank_stable = bool(
        len(convergence_rows) >= 2
        and convergence_rows[-1].estimated_rank == convergence_rows[-2].estimated_rank
    )
    final_s0_change = (
        convergence_rows[-1].s0_relative_change_from_previous
        if convergence_rows
        else math.nan
    )

    candidate_failures=[]
    sigma_tasks: list[tuple[DiscoverySeed, tuple[float, float]]] = []
    for seed in seeds:
        try:
            sigma_tasks.append((seed, sigma_bracket_from_seed(seed, config)))
        except ValueError as exc:
            candidate_failures.append({"seed":asdict(seed),"error":str(exc)})
            continue

    if not sigma_tasks and final_diag.estimated_rank == 1 and rank_stable:
        fallback = global_sigma_fallback_seed(epsilon, config)
        if fallback is not None:
            sigma_tasks.append(fallback)

    refinement_rows: list[ModeRefinementRow] = []
    modes: list[ModeResult] = []
    for idx, (seed, bracket) in enumerate(sigma_tasks, start=1):
        try:
            rows, mode = run_candidate_refinement(epsilon, idx, seed, bracket, config)
        except Exception as exc:
            candidate_failures.append({"seed":asdict(seed),"error":f"{type(exc).__name__}: {exc}"})
            continue
        refinement_rows.extend(rows)
        modes.append(mode)

    resolved=[]
    for mode in modes:
        if mode.resolved and not any(abs(mode.sigma_numerical-old.sigma_numerical)
                <=1e-5*max(mode.sigma_numerical,old.sigma_numerical) for old in resolved):
            resolved.append(mode)

    # Tighter upper-threshold contour check.
    tighter_consistent: bool | None = None
    if config.beyn_check_tighter_cutoff_margin:
        tighter_margin = max(
            config.beyn_min_cutoff_margin,
            cutoff_margin * config.beyn_tighter_cutoff_margin_factor,
        )
        if tighter_margin < cutoff_margin * (1.0 - 1.0e-12):
            tight_eigs, _, tight_diag, tight_rows, _ = beyn_discover(
                epsilon,
                final_diag.quadrature_points,
                mu,
                config,
                cutoff_margin_override=tighter_margin,
            )
            tight_candidates,_=candidate_values(tight_eigs,tight_rows,config)
            original_candidates,_=candidate_values(final_eigs,final_rows,config)
            tighter_consistent = bool(tight_diag.estimated_rank==tight_diag.threshold_rank==tight_diag.gap_rank==final_diag.estimated_rank
                                      and tight_diag.max_linear_solve_relative_residual<=config.beyn_linear_solve_tolerance
                                      and candidate_set_shift(original_candidates,tight_candidates)<=config.beyn_candidate_convergence_tolerance)

    one_mode_supported = bool(
        rank_stable
        and final_diag.estimated_rank == 1
        and len(resolved) == 1
        and bool(tighter_consistent)
        and final_diag.moments_converged and final_diag.candidate_stable
        and final_diag.threshold_rank==final_diag.gap_rank==final_diag.estimated_rank
        and final_diag.max_linear_solve_relative_residual<=config.beyn_linear_solve_tolerance
        and not candidate_failures and all(m.resolved for m in modes)
    )

    if resolved:
        best = resolved[0]
        kb_num = best.kb_numerical
        sigma_num = best.sigma_numerical
        k2_num = (kb_num / config.b) ** 2
        rel_sigma_error = abs(sigma_num - sigma_asym) / max(abs(sigma_asym), 1.0e-30)
        scaled_remainder = abs(sigma_num - sigma_asym) / max(
            epsilon**3 * abs(math.log(epsilon)), 1.0e-30
        )
        parity_comm = parity_commutator_residual(
            kb_num, epsilon, config.bem_M, config
        )
        physical = (
            physical_mode_diagnostics(
                epsilon, kb_num, config.bem_M, config
            )
            if config.run_physical_diagnostics
            else None
        )
        final_mode = best
    else:
        kb_num = sigma_num = k2_num = rel_sigma_error = scaled_remainder = math.nan
        parity_comm = math.nan
        physical = None
        final_mode = None

    physical_verified: bool | None
    if physical is None:
        physical_verified = None
    else:
        physical_verified = bool(
            physical.walls_verified
            and physical.centerline_verified
            and physical.parity_verified
            and physical.open_channel_suppressed
            and physical.boundary_integral_verified
            and physical.decay_verified
        )

    result = ValidationResult(
        epsilon=float(epsilon),
        beta=float(config.shape_beta),
        mu=float(mu),
        nu=float(nu),
        asymptotic_coefficient=float(coefficient),
        lambda_1=lambda_1(config),
        lambda_2=lambda_2(config),
        kb_lower_threshold=kb_1(config),
        kb_upper_threshold=kb_2(config),
        effective_cutoff_margin=float(cutoff_margin),
        kb_asymptotic=float(kb_asym),
        sigma_asymptotic=float(sigma_asym),
        beyn_final_quadrature_points=int(final_diag.quadrature_points),
        beyn_estimated_rank=int(final_diag.estimated_rank),
        beyn_rank_stable=rank_stable,
        beyn_final_s0_relative_change=float(final_s0_change),
        tighter_margin_rank_consistent=tighter_consistent,
        local_refinement_seed_count=len(sigma_tasks),
        resolved_mode_count=len(resolved),
        kb_numerical=float(kb_num),
        k2_numerical=float(k2_num),
        sigma_numerical=float(sigma_num),
        sigma_over_epsilon_squared=(
            float(sigma_num / epsilon**2) if np.isfinite(sigma_num) else math.nan
        ),
        scaled_asymptotic_remainder=float(scaled_remainder),
        relative_error_sigma=float(rel_sigma_error),
        sigma_min_final=(
            float(final_mode.sigma_min_final) if final_mode is not None else math.nan
        ),
        relative_singular_value_final=(
            float(final_mode.relative_singular_value_final)
            if final_mode is not None
            else math.nan
        ),
        final_drop_factor=(
            float(final_mode.final_drop_factor) if final_mode is not None else math.nan
        ),
        final_relative_mesh_change=(
            float(final_mode.final_relative_mesh_change)
            if final_mode is not None
            else math.nan
        ),
        operator_parity_commutator_residual=float(parity_comm),
        unique_embedded_mode_supported=bool(
            one_mode_supported
            and np.isfinite(k2_num)
            and lambda_1(config) <= k2_num < lambda_2(config)
            and bool(physical_verified) and parity_comm<=config.parity_operator_tolerance
        ),
        wall_relative_residual=(
            physical.wall_relative_residual if physical is not None else math.nan
        ),
        centerline_relative_residual=(
            physical.centerline_relative_residual if physical is not None else math.nan
        ),
        odd_parity_relative_residual=(
            physical.odd_parity_relative_residual if physical is not None else math.nan
        ),
        first_open_channel_relative_amplitude=(
            physical.first_open_channel_relative_amplitude
            if physical is not None
            else math.nan
        ),
        boundary_integral_relative_residual=(
            physical.boundary_integral_relative_residual
            if physical is not None
            else math.nan
        ),
        decay_rate_left=(physical.decay_rate_left if physical is not None else math.nan),
        decay_rate_right=(physical.decay_rate_right if physical is not None else math.nan),
        physical_embedded_mode_verified=physical_verified,
        mode_exists_supported=bool(resolved and physical_verified and parity_comm<=config.parity_operator_tolerance),
        odd_uniqueness_supported=bool(one_mode_supported and physical_verified and parity_comm<=config.parity_operator_tolerance),
        full_problem_uniqueness_supported=False,
        diagnostic_complete=bool(not candidate_failures and all(m.resolved for m in modes)
                                 and final_diag.moments_converged),
        beyn_threshold_rank=final_diag.threshold_rank,beyn_gap_rank=final_diag.gap_rank,
        beyn_final_s1_relative_change=final_diag.s1_relative_change,
        beyn_moments_converged=final_diag.moments_converged,
        normalized_residual=final_mode.normalized_residual if final_mode else math.nan,
        equilibrated_residual=final_mode.equilibrated_residual if final_mode else math.nan,
        equilibrated_sigma_min=final_mode.equilibrated_sigma_min if final_mode else math.nan,
        full_vector_residual=final_mode.full_vector_residual if final_mode else math.nan,
        parity_leakage_residual=final_mode.parity_leakage_residual if final_mode else math.nan,
    )
    directory=Path(config.output_directory)/'geometry_diagnostics'/f'epsilon_{epsilon:.12g}'
    directory.mkdir(parents=True,exist_ok=True)
    harness.atomic_json(directory/'candidate_failures.json',candidate_failures)
    harness.atomic_json(directory/'certificate.json',asdict(result))
    write_dataclass_csv(directory/'beyn_levels.csv',diagnostics)
    write_dataclass_csv(directory/'candidate_modes.csv',modes)
    write_dataclass_csv(directory/'physical.csv',[physical] if physical is not None else [])
    return (
        result,
        diagnostics,
        all_eigen_rows,
        convergence_rows,
        refinement_rows,
        modes,
        physical,
    )


# ---------------------------------------------------------------------------
# Numerical convergence study
# ---------------------------------------------------------------------------


def run_internal_convergence_study(
    results: list[ValidationResult],
    config: Config,
) -> list[InternalConvergenceRow]:
    if not config.run_internal_convergence_study:
        return []
    finite = [r for r in results if np.isfinite(r.sigma_numerical)]
    if not finite:
        return []
    target = min(
        finite,
        key=lambda r: abs(r.epsilon - config.internal_convergence_epsilon),
    )
    eps = float(target.epsilon)
    factor = max(config.internal_sigma_bracket_factor, 1.001)
    global_left, global_right = sigma_target_band_limits(config)
    left = max(global_left, target.sigma_numerical / factor)
    right = min(global_right, target.sigma_numerical * factor)
    if not left < right:
        return []

    reference_M = int(config.bem_M)
    ref_kb, ref_sigma, _, _, ref_rel, _, _ = refine_candidate_for_M(
        eps, reference_M, left, right, config
    )

    rows: list[InternalConvergenceRow] = []

    def add(parameter: str, value: float, M: int, cfg: Config) -> None:
        kb, sigma, _, _, rel, _, _ = refine_candidate_for_M(
            eps, int(M), left, right, cfg
        )
        rows.append(
            InternalConvergenceRow(
                epsilon=eps,
                parameter=parameter,
                value=float(value),
                full_M=int(M),
                kb=float(kb),
                sigma_bem=float(sigma),
                relative_singular_value=float(rel),
                relative_kb_shift_from_reference=abs(kb - ref_kb) / max(abs(ref_kb), 1.0e-30),
                relative_sigma_shift_from_reference=abs(sigma - ref_sigma) / max(abs(ref_sigma), 1.0e-30),
            )
        )

    for h in config.finite_difference_steps_test:
        add(
            "finite_difference_step",
            float(h),
            reference_M,
            replace(config, finite_difference_step=float(h)),
        )
    for terms in config.lattice_terms_test:
        add(
            "lattice_terms",
            float(terms),
            reference_M,
            replace(config, lattice_terms=int(terms)),
        )
    for order in config.harmonic_orders_test:
        add(
            "harmonic_order",
            float(order),
            reference_M,
            replace(config, harmonic_order=int(order)),
        )
    return rows


# ---------------------------------------------------------------------------
# Reporting and plots
# ---------------------------------------------------------------------------


def write_dataclass_csv(path: Path, rows: Sequence[Any]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0].keys()))
        writer.writeheader()
        writer.writerows(data)


def plot_shape(config: Config, output: Path) -> None:
    t = np.linspace(0.0, 2.0 * PI, 1000)
    X, Y, *_ = reference_obstacle_geometry(t, config)
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.plot(X, Y, lw=1.8)
    ax.axhline(0.0, lw=0.8)
    ax.axvline(0.0, lw=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"Inflated obstacle, beta={config.shape_beta:g}")
    fig.tight_layout()
    fig.savefig(output / "inflated_obstacle.png", dpi=config.plot_dpi)
    plt.close(fig)


def plot_mu_convergence(rows: list[MuMFSRow], output: Path, config: Config) -> None:
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot([r.order for r in rows], [abs(r.mu-rows[-1].mu) for r in rows], marker="o")
    ax.set_yscale("symlog",linthresh=1e-13)
    ax.set_xlabel("MFS order")
    ax.set_ylabel("|mu_N - mu_last| (zero retained)")
    ax.set_title("Independent dipole-strength convergence")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "mu_mfs_convergence.png", dpi=config.plot_dpi)
    plt.close(fig)


def plot_summary(results: list[ValidationResult], output: Path, config: Config) -> None:
    finite = [r for r in results if np.isfinite(r.sigma_numerical)]
    if not finite:
        return
    eps = np.array([r.epsilon for r in finite])
    kb_num = np.array([r.kb_numerical for r in finite])
    kb_asym = np.array([r.kb_asymptotic for r in finite])
    ratio = np.array([r.sigma_over_epsilon_squared for r in finite])
    C = finite[0].asymptotic_coefficient
    Q = np.array([r.scaled_asymptotic_remainder for r in finite])

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(eps, kb_num, marker="o", label="BEM/SVD")
    ax.plot(eps, kb_asym, linestyle="--", label="leading asymptotic")
    ax.axhline(kb_2(config), linestyle=":", label="sqrt(Lambda_2)b")
    ax.set_xlabel("epsilon")
    ax.set_ylabel("kb")
    ax.set_title("Theorem 2.3(iii): embedded eigenvalue branch")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "kb_vs_epsilon.png", dpi=config.plot_dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(eps, ratio, marker="o", label="sigma_BEM / epsilon^2")
    ax.axhline(C, linestyle="--", label="pi^3 mu / b^3")
    ax.set_xlabel("epsilon")
    ax.set_ylabel("sigma / epsilon^2")
    ax.set_title("Leading-order scaling")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "sigma_over_epsilon_squared.png", dpi=config.plot_dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(eps, Q, marker="o")
    ax.set_xlabel("epsilon")
    ax.set_ylabel("|sigma-C eps^2| / (eps^3 |log eps|)")
    ax.set_title("Scaled asymptotic remainder")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "scaled_remainder.png", dpi=config.plot_dpi)
    plt.close(fig)

    small = [r for r in finite if r.epsilon <= config.publication_asymptotic_epsilon_max]
    if small:
        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        ax.plot(
            [r.epsilon for r in small],
            [r.sigma_over_epsilon_squared for r in small],
            marker="o",
        )
        ax.axhline(C, linestyle="--")
        ax.set_xlabel("epsilon")
        ax.set_ylabel("sigma / epsilon^2")
        ax.set_title("Small-epsilon asymptotic window")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(output / "publication_small_epsilon_scaling.png", dpi=config.plot_dpi)
        plt.close(fig)


def plot_internal_convergence(
    rows: list[InternalConvergenceRow],
    output: Path,
    config: Config,
) -> None:
    parameters = sorted({r.parameter for r in rows})
    for parameter in parameters:
        subset = [r for r in rows if r.parameter == parameter]
        if not subset:
            continue
        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        ax.plot(
            [r.value for r in subset],
            [r.relative_sigma_shift_from_reference for r in subset],
            marker="o",
        )
        if all(r.value > 0 for r in subset):
            ax.set_xscale("log" if parameter != "M" else "linear")
        ax.set_yscale("symlog", linthresh=1e-12)
        ax.set_xlabel(parameter)
        ax.set_ylabel("relative shift in sigma (zero retained)")
        ax.set_title(f"Numerical convergence: {parameter}")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            output / f"internal_convergence_{parameter}.png",
            dpi=config.plot_dpi,
        )
        plt.close(fig)


def build_publication_claims(
    shape: ShapeDiagnostics,
    mu: float,
    nu: float,
    results: list[ValidationResult],
    convergence: list[InternalConvergenceRow],
    config: Config,
) -> list[PublicationClaimRow]:
    claims: list[PublicationClaimRow] = []

    symmetry_ok = bool(
        abs(config.a) <= 100 * np.finfo(float).eps
        and shape.x_even_defect <= 1.0e-12
        and shape.y_odd_defect <= 1.0e-12
        and abs(nu) <= 1.0e-7 * max(mu, 1.0)
    )
    claims.append(
        PublicationClaimRow(
            claim="Theorem 2.3(iii) symmetry assumptions (a=0, X even, Y odd, nu=0)",
            status="supported" if symmetry_ok else "not-established",
            evidence=(
                f"a={config.a:.3e}, X-even defect={shape.x_even_defect:.3e}, "
                f"Y-odd defect={shape.y_odd_defect:.3e}, MFS nu={nu:.3e}, mu={mu:.8g}."
            ),
            limitation="nu=0 is checked numerically in addition to the exact parametrization symmetry.",
        )
    )

    small = [
        r for r in results
        if r.epsilon <= config.publication_asymptotic_epsilon_max
    ]
    unique = sum(r.odd_uniqueness_supported for r in small)
    physical = sum(bool(r.physical_embedded_mode_verified) for r in small)
    existence_status = (
        "supported"
        if len(small) >= config.publication_min_asymptotic_points
        and unique == len(small) and physical==len(small)
        else "partial"
    )
    claims.append(
        PublicationClaimRow(
            claim="One physical mode; uniqueness supported only in the odd sector and sampled contour",
            status=existence_status,
            evidence=(
                f"{unique}/{len(small)} small-epsilon cases have stable odd-sector Beyn rank one "
                f"and exactly one resolved fixed-M SVD mode; {physical}/{len(small)} pass the "
                "embedded-field diagnostics."
            ),
            limitation="Finite sampling supports but does not prove theorem-wide uniqueness.",
        )
    )

    finite_small = [r for r in small if np.isfinite(r.sigma_numerical)]
    if len(finite_small) >= config.publication_min_asymptotic_points:
        ordered = sorted(finite_small, key=lambda r: r.epsilon)
        C = ordered[0].asymptotic_coefficient
        errors = [abs(r.sigma_over_epsilon_squared - C) for r in ordered]
        q_values = [r.scaled_asymptotic_remainder for r in ordered]
        trend = errors[0] <= errors[-1]
        bounded_q = all(np.isfinite(q_values))
        asym_status = "consistent-with-leading-term" if trend and bounded_q and all(r.mode_exists_supported for r in ordered) else "not-established"
        evidence = (
            f"At the smallest sampled epsilon, sigma/epsilon^2={ordered[0].sigma_over_epsilon_squared:.8g} "
            f"versus C=pi^3 mu/b^3={C:.8g}; scaled remainders are finite over the reporting window."
        )
    else:
        asym_status = "not-established"
        evidence = "Insufficient resolved small-epsilon points."
    claims.append(
        PublicationClaimRow(
            claim="sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)",
            status=asym_status,
            evidence=evidence,
            limitation="Bounded finite-sample scaled remainder is numerical evidence, not a proof of Big-O.",
        )
    )

    claims.append(
        PublicationClaimRow(
            claim="Numerical implementation sensitivity at fixed M",
            status="documented" if convergence else "not-run",
            evidence=(
                f"Representative one-at-a-time sensitivity at fixed M={config.bem_M} in "
                "finite-difference step, lattice truncation, and harmonic order is saved to internal_convergence.csv."
                if convergence
                else "Convergence study disabled or unavailable."
            ),
            limitation="Representative convergence is not repeated at every epsilon.",
        )
    )

    claims.append(
        PublicationClaimRow(
            claim="Analyticity in epsilon and epsilon log epsilon",
            status="theoretical-only",
            evidence="The paper proves analyticity; finite numerical samples cannot establish analyticity.",
            limitation="Not a numerically provable statement from a finite sweep.",
        )
    )
    return claims


def write_publication_report(
    output: Path,
    shape: ShapeDiagnostics,
    mu_rows: list[MuMFSRow],
    circle: MuMFSRow,
    results: list[ValidationResult],
    claims: list[PublicationClaimRow],
    config: Config,
) -> None:
    lines: list[str] = []
    lines.append("# Numerical validation report — Theorem 2.3(iii)")
    lines.append("")
    lines.append("This run uses the exact x-axis symmetry reduction (odd in y) described in Remark 2.5. Full-problem uniqueness remains not established: the even sector and omitted thresholds are not counted. Finite-sample asymptotic trends are not a proof of Big-O.")
    lines.append("")
    lines.append("## Geometry and dipole strength")
    lines.append("")
    lines.append(f"- beta = {config.shape_beta:g}")
    lines.append(f"- a = {config.a:g}")
    lines.append(f"- fixed physical boundary order M = {config.bem_M}; Beyn/SVD use the odd reduced operator")
    lines.append(f"- X-even defect = {shape.x_even_defect:.3e}")
    lines.append(f"- Y-odd defect = {shape.y_odd_defect:.3e}")
    lines.append(f"- minimum reference speed = {shape.min_reference_speed:.6e}")
    lines.append(f"- circle calibration mu = {circle.mu:.12g}, nu = {circle.nu:.3e}")
    if mu_rows:
        lines.append(f"- final shape mu = {mu_rows[-1].mu:.12g}, nu = {mu_rows[-1].nu:.3e}")
        lines.append(f"- C = pi^3 mu/b^3 = {asymptotic_coefficient(mu_rows[-1].mu, config):.12g}")
    lines.append("")
    lines.append("## Small-epsilon branch")
    lines.append("")
    lines.append("| epsilon | odd Beyn rank | resolved modes | kb | sigma/eps^2 | Q(eps) | physical |")
    lines.append("|---:|---:|---:|---:|---:|---:|:---:|")
    for r in results:
        if r.epsilon <= config.publication_asymptotic_epsilon_max:
            lines.append(
                f"| {r.epsilon:.4f} | {r.beyn_estimated_rank} | {r.resolved_mode_count} | "
                f"{r.kb_numerical:.12g} | {r.sigma_over_epsilon_squared:.8g} | "
                f"{r.scaled_asymptotic_remainder:.6g} | "
                f"{'PASS' if r.physical_embedded_mode_verified else 'FAIL'} |"
            )
    lines.append("")
    lines.append("## Claim audit")
    lines.append("")
    for c in claims:
        lines.append(f"- **{c.status}** — {c.claim}: {c.evidence} Limitation: {c.limitation}")
    lines.append("")
    lines.append("## Interpretation boundary")
    lines.append("")
    lines.append(
        "The odd-parity reduction is essential: it removes the first even propagating channel, "
        "turning the embedded full-space eigenvalue into a discrete eigenvalue of the odd sector. "
        "The full reconstructed field is then checked to be odd, to suppress the first open channel, "
        "and to decay exponentially."
    )
    lines.append(
        "Points above the configured small-epsilon reporting window are finite-size exploration and "
        "are not used to define the formal asymptotic claim."
    )
    (output / "publication_validation_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def print_result(result: ValidationResult) -> None:
    print("\n  --- Theorem 2.3(iii) validation result ---")
    print(f"  epsilon                         = {result.epsilon:.6f}")
    print(f"  odd-sector Beyn rank            = {result.beyn_estimated_rank}")
    print(f"  Beyn rank stable                = {'YES' if result.beyn_rank_stable else 'no'}")
    print(f"  final Beyn Nq                   = {result.beyn_final_quadrature_points}")
    print(f"  tighter-cutoff rank check       = {result.tighter_margin_rank_consistent}")
    print(f"  resolved odd modes              = {result.resolved_mode_count}")
    print(f"  kb asymptotic                   = {result.kb_asymptotic:.12f}")
    print(f"  kb BEM                          = {result.kb_numerical:.12f}")
    print(f"  sigma asym                      = {result.sigma_asymptotic:.8e}")
    print(f"  sigma BEM                       = {result.sigma_numerical:.8e}")
    print(f"  sigma/epsilon^2                 = {result.sigma_over_epsilon_squared:.8e}")
    print(f"  C=pi^3 mu/b^3                   = {result.asymptotic_coefficient:.8e}")
    print(f"  scaled asymptotic remainder     = {result.scaled_asymptotic_remainder:.8e}")
    print(f"  relative sigma error            = {result.relative_error_sigma:.3%}")
    print(f"  relative singular value         = {result.relative_singular_value_final:.3e}")
    print(f"  minimum drop                    = {result.final_drop_factor:.3e}")
    print("  cross-M change                      = N/A (single fixed M)")
    print(f"  parity commutator residual      = {result.operator_parity_commutator_residual:.3e}")
    print(f"  one embedded mode supported     = {'PASS' if result.unique_embedded_mode_supported else 'FAIL'}")
    if np.isfinite(result.wall_relative_residual):
        print(f"  wall residual                   = {result.wall_relative_residual:.3e}")
        print(f"  centerline residual             = {result.centerline_relative_residual:.3e}")
        print(f"  odd-parity residual             = {result.odd_parity_relative_residual:.3e}")
        print(f"  first open-channel amplitude    = {result.first_open_channel_relative_amplitude:.3e}")
        print(f"  off-grid BIE residual           = {result.boundary_integral_relative_residual:.3e}")
        print(f"  decay left/right                = {result.decay_rate_left:.6e}/{result.decay_rate_right:.6e}")
        print(
            f"  physical embedded-mode checks   = "
            f"{'PASS' if result.physical_embedded_mode_verified else 'FAIL'}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _run_validation(config: Config) -> None:
    output = Path(config.output_directory)
    output.mkdir(parents=True, exist_ok=True)

    if abs(config.a) > 100.0 * np.finfo(float).eps:
        raise RuntimeError("Theorem 2.3(iii) requires a=0 exactly.")
    if abs(config.shape_beta) >= 1.0:
        raise RuntimeError(
            "For this validation script use |shape_beta|<1 so the supplied smooth "
            "parametrization remains safely regular/star-shaped for the MFS sources."
        )
    for M in (config.bem_M,):
        if int(M) % 2:
            raise RuntimeError("Every full boundary order M must be even.")

    print("=== Theorem 2.3(iii) numerical validation: odd-sector Beyn + BEM + sigma-SVD ===")
    print(f"b = {config.b}")
    print(f"a = {config.a} (required exactly by statement iii)")
    print(f"shape beta = {config.shape_beta}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"Lambda_2 = {lambda_2(config):.12f}")
    print(f"sqrt(Lambda_1)b = {kb_1(config):.12f}")
    print(f"sqrt(Lambda_2)b = {kb_2(config):.12f}")

    shape = compute_shape_diagnostics(config)
    write_dataclass_csv(output / "shape_diagnostics.csv", [shape])
    plot_shape(config, output)
    print("\n=== SHAPE / SYMMETRY CHECK ===")
    print(f"reference area = {shape.area:.12g}")
    print(f"min |r'(t)| = {shape.min_reference_speed:.6e}")
    print(f"X-even defect = {shape.x_even_defect:.3e}")
    print(f"Y-odd defect = {shape.y_odd_defect:.3e}")

    print("\n=== INDEPENDENT DIPOLE-STRENGTH mu (MFS) ===")
    mu, nu, mu_rows, circle = run_mu_convergence(config)
    write_dataclass_csv(output / "mu_mfs_convergence.csv", mu_rows)
    write_dataclass_csv(output / "mu_circle_calibration.csv", [circle])
    plot_mu_convergence(mu_rows, output, config)
    print(f"circle calibration: mu={circle.mu:.12g}, nu={circle.nu:.3e}")
    for row in mu_rows:
        print(
            f"  N={row.order:>4d}: mu={row.mu:.12g}, nu={row.nu:.3e}, "
            f"boundary_res={row.boundary_relative_residual:.3e}"
        )
    print(f"selected mu = {mu:.12g}")
    print(f"selected nu = {nu:.3e} (must vanish by x-axis symmetry)")
    print(f"C = pi^3 mu/b^3 = {asymptotic_coefficient(mu, config):.12g}")

    # Complex-kb smoke test on the reduced odd operator.
    eps0 = config.epsilon_values[0]
    left, right, center, _, ry, _ = contour_geometry(eps0, mu, config)
    ztest = complex(center, 0.37 * ry)
    Atest = assemble_odd_matrix(ztest, eps0, config.bem_M, config)
    if not np.all(np.isfinite(Atest)):
        raise RuntimeError("Complex-kb odd-sector BEM smoke test produced non-finite values.")
    print("\ncomplex-kb odd-sector assembly: PASS")
    print(f"target real band ~ [{left:.12f}, {right:.12f}]")

    summaries: list[ValidationResult] = []
    all_beyn_diags: list[BeynDiagnostics] = []
    all_eig_rows: list[BeynEigenvalueRow] = []
    all_beyn_conv: list[BeynConvergenceRow] = []
    all_refinement: list[ModeRefinementRow] = []
    all_modes: list[ModeResult] = []
    all_physical: list[PhysicalDiagnostics] = []

    print("\n=== EMBEDDED-MODE EPSILON SWEEP ===")
    print(f"epsilons = {config.epsilon_values}")
    print(
        f"Beyn full M={config.bem_M} -> odd reduced dimension={config.bem_M//2}; "
        f"Nq base={config.beyn_quadrature_levels}, adaptive={config.beyn_adaptive_quadrature_levels}"
    )
    print(f"fixed full BEM/SVD order M = {config.bem_M}")

    module = sys.modules[__name__]
    cases = harness.independent_cases(module, config, mu, nu)
    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.8f} ===")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, mu, config)
        margin = effective_cutoff_margin(epsilon, mu, config)
        print(f"predicted kb = {kb_asym:.12f}")
        print(f"predicted sigma = {sigma_asym:.8e}")
        print(f"predicted Lambda_2 cutoff gap in kb = {kb_2(config)-kb_asym:.3e}")
        print(f"effective Beyn upper margin = {margin:.3e}")

        (
            result,
            diags,
            eigrows,
            convrows,
            refrows,
            modes,
            physical,
        ) = cases[epsilon]
        summaries.append(result)
        all_beyn_diags.extend(diags)
        all_eig_rows.extend(eigrows)
        all_beyn_conv.extend(convrows)
        all_refinement.extend(refrows)
        all_modes.extend(modes)
        if physical is not None:
            all_physical.append(physical)
        print_result(result)

    write_dataclass_csv(output / "summary.csv", summaries)
    write_dataclass_csv(output / "beyn_diagnostics.csv", all_beyn_diags)
    write_dataclass_csv(output / "beyn_raw_eigenvalues.csv", all_eig_rows)
    write_dataclass_csv(output / "beyn_convergence.csv", all_beyn_conv)
    write_dataclass_csv(output / "mode_refinement.csv", all_refinement)
    write_dataclass_csv(output / "mode_results.csv", all_modes)
    write_dataclass_csv(output / "physical_diagnostics.csv", all_physical)
    plot_summary(summaries, output, config)

    convergence = run_internal_convergence_study(summaries, config)
    write_dataclass_csv(output / "internal_convergence.csv", convergence)
    plot_internal_convergence(convergence, output, config)

    small = [
        r for r in summaries
        if r.epsilon <= config.publication_asymptotic_epsilon_max
    ]
    write_dataclass_csv(output / "publication_asymptotic_window.csv", small)

    claims = build_publication_claims(shape, mu, nu, summaries, convergence, config)
    write_dataclass_csv(output / "publication_claims.csv", claims)
    write_publication_report(output, shape, mu_rows, circle, summaries, claims, config)

    print("\n=== FINAL SUMMARY ===")
    unique = sum(r.unique_embedded_mode_supported for r in summaries)
    physical_ok = sum(bool(r.physical_embedded_mode_verified) for r in summaries)
    print(f"stable odd-sector rank=1 + exactly one resolved mode: {unique}/{len(summaries)}")
    print(f"full embedded-field physical diagnostics: {physical_ok}/{len(summaries)}")
    print("Publication claim audit:")
    for claim in claims:
        print(f"  - {claim.status:>16}: {claim.claim}")
    print(f"\nFiles written to: {output.resolve()}")


def main() -> None:
    harness.run_entrypoint(sys.modules[__name__], CONFIG, _run_validation)


if __name__ == "__main__":
    main()
