from __future__ import annotations

import csv
import math
import os
import sys
from collections.abc import Callable
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

# Same project import convention as the existing scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice

PI = np.pi


@dataclass(frozen=True)
class Config:
    # Geometry.
    b: float = 1.0
    a: float = 0.6
    epsilon_values: tuple[float, ...] = (
        0.01, 0.03111111, 0.05222222, 0.07333333, 0.09444444,
        0.11555556, 0.13666667, 0.15777778, 0.17888889, 0.2,
    )

    # Green function / BEM.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Internal numerical-convergence diagnostics. These do NOT alter the main
    # mathematical discretization; they probe sensitivity to implementation
    # parameters at one representative epsilon.
    run_internal_convergence_study: bool = True
    internal_convergence_epsilon: float = 0.09
    internal_convergence_M: int = 32
    internal_refine_half_width: float = 7.5e-5
    finite_difference_steps_test: tuple[float, ...] = (3.0e-6, 1.0e-6, 3.0e-7)
    lattice_terms_test: tuple[int, ...] = (100, 200, 300)
    harmonic_orders_test: tuple[int, ...] = (12, 20, 28)

    # ------------------------------------------------------------------
    # Beyn global discovery.
    # ------------------------------------------------------------------
    # BEM order used only for contour discovery. Local candidates are then
    # validated at the higher refinement_M values below.
    beyn_M: int = 24

    # Base quadrature-convergence study. These levels are always run so the
    # rank and raw Beyn eigenvalue sequence can be diagnosed consistently.
    beyn_quadrature_levels: tuple[int, ...] = (96, 192, 384)

    # Adaptive continuation. These expensive levels are used only when the base
    # study has not produced a strict in-contour real candidate or the final
    # Beyn rank is not yet stable. This is especially useful for modes extremely
    # close to a cutoff, where contour quadrature converges slowly.
    beyn_adaptive_quadrature_levels: tuple[int, ...] = (768, 1536)

    # Number of probing vectors. Must exceed the number of eigenvalues enclosed
    # by the contour (counting algebraic multiplicity). If the detected rank is
    # close to this value, the script warns that probe_dim should be increased.
    beyn_probe_dim: int = 8
    beyn_random_seed: int = 1729

    # Rank(S0) is diagnosed in two complementary ways:
    #   (a) a conventional relative singular-value threshold, and
    #   (b) the dominant spectral gap among those provisionally retained
    #       singular values.  The gap estimate prevents quadrature noise near a
    #       cutoff from being automatically interpreted as an extra mode.
    beyn_rank_relative_tolerance: float = 1.0e-8
    beyn_rank_absolute_tolerance: float = 1.0e-12
    beyn_rank_gap_threshold: float = 1.0e3

    # Critical fix: permit an EMPTY contour.  If every singular value of S0 is
    # at the numerical-noise level, the estimated rank is zero rather than
    # being forcibly promoted to one.  This is essential for testing the
    # non-existence side of Theorem 2.1.
    beyn_empty_s0_tolerance: float = 1.0e-9

    # The contour is an ellipse enclosing the real discrete band. It stays a
    # small positive distance away from k=0 and from the first cutoff.
    beyn_low_k_margin: float = 1.0e-4
    beyn_cutoff_margin: float = 1.0e-5
    beyn_min_cutoff_margin: float = 1.0e-10
    # The effective right-end margin is reduced when the asymptotic branch is
    # closer to the cutoff.  The contour still discovers globally; the
    # asymptotic formula only prevents us from accidentally chopping off the
    # near-cutoff region we intend to test.
    beyn_cutoff_margin_fraction_of_predicted_gap: float = 0.25
    # Independent contour-geometry diagnostic: repeat the FINAL Beyn solve with
    # a tighter cutoff margin and check that the enclosed rank does not change.
    beyn_check_tighter_cutoff_margin: bool = True
    beyn_tighter_cutoff_margin_factor: float = 0.20
    beyn_imag_half_height: float = 3.0e-2

    # A true discrete trapped mode is real. Beyn may return a small imaginary
    # part from quadrature/discretization error. Candidates below this threshold
    # are projected to the real axis and locally verified with SVD.
    beyn_real_axis_tolerance: float = 5.0e-5

    # Cluster duplicate near-real Beyn eigenvalues if numerical multiplicity
    # produces nearly identical values.
    beyn_cluster_tolerance: float = 5.0e-5

    # A raw Beyn eigenvalue that is real to numerical precision but lies only a
    # tiny distance outside the real span of the contour is retained as a
    # *fallback seed only*. It is never accepted as an eigenvalue by itself; the
    # local SVD stage must still resolve a genuine interior minimum.
    beyn_near_contour_real_tolerance: float = 5.0e-5

    # If adaptive quadrature still leaves no strict candidate, an Aitken
    # Delta^2 extrapolation of the last three raw rank-one Beyn estimates may be
    # used as a local-SVD seed. This remains a diagnostic/extrapolation, never a
    # certified eigenvalue.
    beyn_enable_aitken_seed: bool = True
    beyn_aitken_require_contracting_shifts: bool = True

    # Only a discovery diagnostic: how much the final Beyn candidate is allowed
    # to move between the last two quadrature levels before we flag slow contour
    # convergence. Local SVD refinement is still the final position estimate.
    beyn_candidate_convergence_tolerance: float = 2.0e-4

    # ------------------------------------------------------------------
    # Local SVD verification of each Beyn-discovered candidate.
    # ------------------------------------------------------------------
    refinement_M: tuple[int, ...] = (16, 24, 32, 40, 48)
    local_refine_half_width: float = 2.0e-3
    minimizer_xatol: float = 2.0e-10

    minimum_drop_factor: float = 100.0
    # Retained as a descriptive diagnostic only.  Certification is based on
    # the scale-invariant ratio sigma_min/sigma_max below.
    near_singular_tolerance: float = 1.0e-4
    relative_near_singular_tolerance: float = 1.0e-4
    mesh_sigma_relative_tolerance: float = 0.03

    # Asymptotic accuracy is intentionally separate from mode existence.
    relative_sigma_error_tolerance: float = 0.05

    # Physical trapped-mode diagnostics.  These are independent of the SVD
    # search: reconstruct the field, verify wall Dirichlet values, estimate the
    # obstacle Neumann residual, and check exponential decay away from the
    # obstacle.  They are reported separately from spectral certification.
    run_physical_diagnostics: bool = True
    physical_quadrature_points: int = 80
    physical_x_over_b: tuple[float, ...] = (1.5, 2.0, 2.5, 3.0)
    physical_neumann_samples: int = 24
    physical_neumann_offset_fraction_of_epsilon: float = 0.04
    physical_wall_relative_tolerance: float = 1.0e-6
    physical_neumann_relative_tolerance: float = 5.0e-2
    physical_decay_relative_tolerance: float = 0.20

    # Full theorem-side diagnostic for the critical height a*.  This is costly
    # because each bisection point launches a complete Beyn+SVD mode-count
    # calculation.  It is implemented here but disabled by default; enable it
    # for the final publication-quality validation run.
    run_critical_height_study: bool = False
    critical_height_epsilon_values: tuple[float, ...] = (0.05, 0.07, 0.09)
    critical_height_half_width: float = 0.12
    critical_height_bisection_iterations: int = 5

    output_directory: str = "theorem_2_1_beyn_v4_full_validation"


@dataclass
class BeynDiagnostics:
    epsilon: float
    M: int
    quadrature_points: int
    probe_dim: int
    threshold_rank: int
    gap_rank: int
    estimated_rank: int
    selected_gap_ratio: float
    leading_s0_singular_value: float
    trailing_kept_s0_singular_value: float
    first_discarded_s0_singular_value: float
    max_linear_solve_relative_residual: float
    raw_eigenvalues: int
    near_real_candidates: int
    near_contour_seeds: int


@dataclass
class BeynEigenvalueRow:
    epsilon: float
    quadrature_points: int
    raw_index: int
    real_part: float
    imag_part: float
    inside_contour: bool
    inside_real_band: bool
    near_real: bool
    real_span_overrun: float
    near_contour_seed: bool
    accepted_for_real_refinement: bool


@dataclass
class DiscoverySeed:
    real: float
    imag: float
    source: str  # strict-beyn | near-contour | aitken


@dataclass
class BeynConvergenceRow:
    epsilon: float
    quadrature_points: int
    estimated_rank: int
    threshold_rank: int
    gap_rank: int
    selected_gap_ratio: float
    s0_sv1: float
    s0_sv2: float
    s0_sv3: float
    s0_relative_change_from_previous: float
    accepted_candidates: int
    near_contour_seeds: int
    primary_candidate_real: float
    primary_candidate_imag: float
    primary_raw_real: float
    primary_raw_imag: float
    max_candidate_shift_from_previous: float
    rank_stable_from_previous: bool


@dataclass
class ModeRefinementRow:
    epsilon: float
    candidate_index: int
    seed_source: str
    beyn_kb_real: float
    beyn_kb_imag: float
    M: int
    kb: float
    sigma_bem: float
    sigma_min: float
    sigma_max: float
    relative_singular_value: float
    left_value: float
    right_value: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_previous: float


@dataclass
class ModeResult:
    epsilon: float
    candidate_index: int
    seed_source: str
    beyn_kb_real: float
    beyn_kb_imag: float
    kb_numerical: float
    sigma_numerical: float
    sigma_min_final: float
    sigma_max_final: float
    relative_singular_value_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    resolved: bool


@dataclass
class PhysicalDiagnostics:
    epsilon: float
    kb: float
    M: int
    wall_relative_residual: float
    neumann_relative_residual: float
    decay_rate_left: float
    decay_rate_right: float
    expected_decay_rate: float
    decay_relative_error_left: float
    decay_relative_error_right: float
    monotone_decay_left: bool
    monotone_decay_right: bool
    walls_verified: bool
    neumann_verified: bool
    decay_verified: bool


@dataclass
class InternalConvergenceRow:
    epsilon: float
    parameter: str
    value: float
    kb: float
    sigma_bem: float
    relative_singular_value: float
    relative_kb_shift_from_baseline: float
    relative_sigma_shift_from_baseline: float


@dataclass
class CriticalHeightRow:
    epsilon: float
    a_lower_zero_mode: float
    a_upper_one_mode: float
    a_critical_estimate: float
    bracket_width: float
    a0_star: float
    normalized_shift_over_epsilon: float
    status: str


@dataclass
class ValidationResult:
    epsilon: float
    a: float
    lambda_1: float
    kb_cutoff: float
    effective_beyn_cutoff_margin: float
    tighter_margin_rank_consistent: bool | None
    kb_asymptotic: float
    sigma_asymptotic: float
    asymptotic_coefficient: float
    beyn_final_quadrature_points: int
    beyn_estimated_rank: int
    beyn_rank_stable: bool
    beyn_near_real_candidates: int
    local_refinement_seed_count: int
    local_refinement_seed_source: str
    aitken_estimate: float
    resolved_mode_count: int
    kb_numerical: float
    sigma_numerical: float
    sigma_over_epsilon_squared: float
    scaled_asymptotic_remainder: float
    sigma_min_final: float
    sigma_max_final: float
    relative_singular_value_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    relative_error_kb: float
    relative_error_sigma: float
    unique_mode_verified: bool
    asymptotic_agreement_verified: bool | None
    wall_relative_residual: float
    neumann_relative_residual: float
    decay_rate_left: float
    decay_rate_right: float
    physical_decay_verified: bool | None


CONFIG = Config()


# ---------------------------------------------------------------------------
# Analytic / asymptotic quantities
# ---------------------------------------------------------------------------


def lambda_1(config: Config) -> float:
    return (PI / (2.0 * config.b)) ** 2


def kb_cutoff(config: Config) -> float:
    return math.sqrt(lambda_1(config)) * config.b


def critical_height_leading_order(config: Config) -> float:
    """Leading a_0* for a circle with R0=1, mu=1, S=pi."""
    mu = 1.0
    area = PI
    return (2.0 * config.b / PI) * math.atan(math.sqrt(area / (2.0 * PI * mu)))


def asymptotic_coefficient(config: Config) -> float:
    """C(a) in sigma = C(a) epsilon^2 + O(epsilon^3 log epsilon)."""
    alpha = PI * config.a / config.b
    mu = 1.0
    area = PI
    bracket = (
        PI * mu * math.sin(alpha / 2.0) ** 2
        - 0.5 * area * math.cos(alpha / 2.0) ** 2
    )
    return PI**2 / (4.0 * config.b**3) * bracket


def effective_cutoff_margin(epsilon: float, config: Config) -> float:
    """
    Pick a cutoff margin small enough not to exclude the near-threshold branch.

    This does not select an eigenvalue; it only chooses how close the global
    contour gets to sqrt(Lambda_1)b.  If the leading coefficient is non-positive
    (the non-existence side of the theorem), use the configured minimum margin.
    """
    coefficient = asymptotic_coefficient(config)
    if coefficient <= 0.0:
        return float(config.beyn_min_cutoff_margin)

    sigma = coefficient * epsilon**2
    k_squared = lambda_1(config) - sigma**2
    if k_squared <= 0.0:
        return float(config.beyn_min_cutoff_margin)

    kb_asym = config.b * math.sqrt(k_squared)
    gap = max(kb_cutoff(config) - kb_asym, config.beyn_min_cutoff_margin)
    margin = config.beyn_cutoff_margin_fraction_of_predicted_gap * gap
    return float(
        min(
            config.beyn_cutoff_margin,
            max(config.beyn_min_cutoff_margin, margin),
        )
    )


def config_for_epsilon(epsilon: float, config: Config) -> Config:
    """Return a frozen Config copy with an epsilon-appropriate cutoff margin."""
    return replace(config, beyn_cutoff_margin=effective_cutoff_margin(epsilon, config))


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    """Return (kb_asymptotic, sigma_asymptotic)."""
    sigma = asymptotic_coefficient(config) * epsilon**2

    if sigma <= 0.0:
        raise ValueError(
            "The leading-order sigma is not positive. The selected geometry "
            "does not satisfy the predicted discrete-mode regime."
        )

    k_squared = lambda_1(config) - sigma**2
    if k_squared <= 0.0:
        raise ValueError("The asymptotic formula produced non-positive k^2.")

    return math.sqrt(k_squared) * config.b, sigma


def sigma_from_kb(kb: float, config: Config) -> float:
    k = kb / config.b
    gap = lambda_1(config) - k**2
    return math.sqrt(max(gap, 0.0))


# ---------------------------------------------------------------------------
# BEM assembly — same mathematical discretization as the current script.
# ---------------------------------------------------------------------------


def boundary_nodes(M: int) -> np.ndarray:
    return (np.arange(M) + 0.5) * (2.0 * PI / M)


def circle_geometry(
    t: np.ndarray | float,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    t = np.asarray(t)
    X = epsilon * np.cos(t)
    Y = epsilon * np.sin(t)
    Xp = -epsilon * np.sin(t)
    Yp = epsilon * np.cos(t)
    Xpp = -epsilon * np.cos(t)
    Ypp = -epsilon * np.sin(t)
    return X, Y, Xp, Yp, Xpp, Ypp


def make_green_functions(
    kb: complex,
    config: Config,
) -> tuple[Callable[..., complex], Callable[..., complex]]:
    """
    Return the same Dirichlet waveguide Green functions, now allowing complex kb.

    Beyn requires analytic continuation of A(kb) away from the real axis. This
    function therefore intentionally does not cast kb or k to float.
    """
    b = config.b
    d = 2.0 * b
    k = kb / b

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
            y + b + config.a,
            xi,
            eta + b + config.a,
            lattice_coefficients,
            k,
            d,
        )

    def green_regularized(x: float, y: float, xi: float, eta: float) -> complex:
        return lattice.greens_dirichlet_reg(
            x,
            y + b + config.a,
            xi,
            eta + b + config.a,
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
    x, y, _, _, _, _ = circle_geometry(psi, epsilon)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = circle_geometry(theta, epsilon)

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
    kb: complex,
    epsilon: float,
    M: int,
    config: Config,
) -> np.ndarray:
    G, G_regularized = make_green_functions(kb, config)
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

    # Full contour [0, 2pi]:
    #   1/2 u_i = (2pi/M) sum_j K^w_ij u_j
    # hence A = I - (4pi/M) K^w.
    return np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted


def singular_metrics(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float]:
    """Return sigma_min(A), sigma_max(A), and the scale-invariant ratio."""
    A = assemble_matrix(complex(kb, 0.0), epsilon, M, config)
    singular_values = np.linalg.svd(A, compute_uv=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    relative = sigma_min / max(sigma_max, np.finfo(float).tiny)
    return sigma_min, sigma_max, float(relative)


def smallest_singular_value(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> float:
    return singular_metrics(kb, epsilon, M, config)[0]


def smallest_singular_pair(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, np.ndarray]:
    """Return singular diagnostics and normalized right singular vector."""
    A = assemble_matrix(complex(kb, 0.0), epsilon, M, config)
    _, singular_values, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= max(float(np.linalg.norm(vector)), 1.0e-30)
    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    relative = sigma_min / max(sigma_max, np.finfo(float).tiny)
    return sigma_min, sigma_max, float(relative), vector


# ---------------------------------------------------------------------------
# Beyn contour method
# ---------------------------------------------------------------------------


def contour_geometry(config: Config) -> tuple[float, float, float, float, float]:
    """Return (left, right, center, radius_x, radius_y) in kb-plane."""
    left = config.beyn_low_k_margin
    right = kb_cutoff(config) - config.beyn_cutoff_margin
    if not (0.0 < left < right < kb_cutoff(config)):
        raise ValueError("Invalid Beyn contour margins.")

    center = 0.5 * (left + right)
    radius_x = 0.5 * (right - left)
    radius_y = config.beyn_imag_half_height
    return left, right, center, radius_x, radius_y


def ellipse_point(theta: float, config: Config) -> tuple[complex, complex]:
    """Return z(theta) and dz/dtheta for the counter-clockwise ellipse."""
    _, _, center, rx, ry = contour_geometry(config)
    z = center + rx * math.cos(theta) + 1j * ry * math.sin(theta)
    dz_dtheta = -rx * math.sin(theta) + 1j * ry * math.cos(theta)
    return complex(z), complex(dz_dtheta)


def point_inside_ellipse(z: complex, config: Config, tolerance: float = 1.0e-9) -> bool:
    _, _, center, rx, ry = contour_geometry(config)
    q = ((z.real - center) / rx) ** 2 + (z.imag / ry) ** 2
    return bool(q <= 1.0 + tolerance)


def complex_support_smoke_test(epsilon: float, config: Config) -> None:
    """Check that A(z) can be assembled at a genuinely complex kb."""
    _, _, center, _, ry = contour_geometry(config)
    z_test = complex(center, 0.37 * ry)
    try:
        A = assemble_matrix(z_test, epsilon, min(8, config.beyn_M), config)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "Beyn requires A(kb) for complex kb, but lattice_sums / Green "
            f"failed at kb={z_test!r}. Original error: "
            f"{type(exc).__name__}: {exc}"
        ) from exc

    if not np.all(np.isfinite(A.real)) or not np.all(np.isfinite(A.imag)):
        raise RuntimeError(
            "Complex-kb smoke test produced non-finite entries in A(kb)."
        )


def probing_matrix(M: int, config: Config) -> np.ndarray:
    if config.beyn_probe_dim >= M:
        raise ValueError("beyn_probe_dim must be strictly smaller than BEM size M.")

    rng = np.random.default_rng(config.beyn_random_seed)
    V = rng.standard_normal((M, config.beyn_probe_dim)) + 1j * rng.standard_normal(
        (M, config.beyn_probe_dim)
    )
    Q, _ = np.linalg.qr(V)
    return Q[:, : config.beyn_probe_dim]


def estimate_beyn_rank(
    singular_values: np.ndarray,
    config: Config,
) -> tuple[int, int, int, float]:
    """
    Estimate rank(S0) while explicitly allowing rank zero.

    The previous implementation forced threshold_rank >= 1, which made an
    empty contour impossible to diagnose numerically.  Here an S0 whose leading
    singular value is at the configured noise floor is classified as rank zero.
    Above that floor we retain the threshold + dominant-gap logic.
    """
    s = np.asarray(singular_values, dtype=float)
    if len(s) == 0 or not np.isfinite(s[0]) or s[0] <= 0.0:
        return 0, 0, 0, math.nan

    if s[0] <= config.beyn_empty_s0_tolerance:
        return 0, 0, 0, math.nan

    threshold = max(
        config.beyn_rank_absolute_tolerance,
        config.beyn_rank_relative_tolerance * s[0],
    )
    threshold_rank = int(np.sum(s > threshold))
    if threshold_rank == 0:
        return 0, 0, 0, math.nan
    threshold_rank = min(threshold_rank, len(s))

    if len(s) == 1:
        return 1, threshold_rank, 1, math.inf

    n_gaps = min(threshold_rank, len(s) - 1)
    if n_gaps <= 0:
        return threshold_rank, threshold_rank, threshold_rank, math.nan

    floor = max(config.beyn_rank_absolute_tolerance, np.finfo(float).tiny)
    ratios = s[:n_gaps] / np.maximum(s[1 : n_gaps + 1], floor)
    best_index = int(np.argmax(ratios))
    selected_gap = float(ratios[best_index])
    gap_rank = best_index + 1

    rank = gap_rank if selected_gap >= config.beyn_rank_gap_threshold else threshold_rank
    return int(rank), int(threshold_rank), int(gap_rank), selected_gap


def beyn_discover(
    epsilon: float,
    quadrature_points: int,
    config: Config,
) -> tuple[
    np.ndarray,
    np.ndarray,
    BeynDiagnostics,
    list[BeynEigenvalueRow],
    np.ndarray,
]:
    """Run one Beyn contour solve at a prescribed quadrature level."""
    M = config.beyn_M
    Nq = int(quadrature_points)
    V = probing_matrix(M, config)

    S0 = np.zeros((M, config.beyn_probe_dim), dtype=np.complex128)
    S1 = np.zeros_like(S0)
    max_solve_residual = 0.0

    print(
        f"  Beyn contour discovery: M={M}, quadrature={Nq}, "
        f"probe_dim={config.beyn_probe_dim}"
    )

    for j in range(Nq):
        # Midpoint-shifted periodic nodes avoid the real-axis ellipse endpoints.
        theta = 2.0 * PI * (j + 0.5) / Nq
        z, dz_dtheta = ellipse_point(theta, config)
        A = assemble_matrix(z, epsilon, M, config)

        try:
            X = np.linalg.solve(A, V)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                f"Linear solve failed on the Beyn contour at kb={z}. "
                "Move the contour away from the spectrum / threshold or "
                "increase numerical resolution."
            ) from exc

        rel_res = np.linalg.norm(A @ X - V) / max(np.linalg.norm(V), 1.0e-30)
        max_solve_residual = max(max_solve_residual, float(rel_res))

        weight = dz_dtheta / (1j * Nq)
        S0 += weight * X
        S1 += weight * z * X

        progress_stride = max(16, Nq // 6)
        if (j + 1) % progress_stride == 0 or j + 1 == Nq:
            print(f"    contour node {j + 1:>3}/{Nq}")

    U, s, Vh = np.linalg.svd(S0, full_matrices=False)
    rank, threshold_rank, gap_rank, selected_gap = estimate_beyn_rank(s, config)

    if rank == 0:
        raw_eigenvalues = np.array([], dtype=np.complex128)
    else:
        Ur = U[:, :rank]
        Wr = Vh[:rank, :].conj().T
        sr = s[:rank]
        B = Ur.conj().T @ S1 @ Wr
        B = B @ np.diag(1.0 / sr)
        raw_eigenvalues = np.linalg.eigvals(B)

    left, right, _, _, _ = contour_geometry(config)
    cutoff = kb_cutoff(config)
    rows: list[BeynEigenvalueRow] = []
    accepted_count = 0

    near_contour_count = 0
    for idx, eig in enumerate(raw_eigenvalues):
        eig = complex(eig)
        inside = point_inside_ellipse(eig, config)
        inside_real_band = bool(0.0 < eig.real < cutoff)
        inside_contour_real_span = bool(left < eig.real < right)
        near_real = bool(abs(eig.imag) <= config.beyn_real_axis_tolerance)
        real_span_overrun = max(left - eig.real, eig.real - right, 0.0)

        accepted = bool(
            inside and inside_real_band and inside_contour_real_span and near_real
        )
        near_contour_seed = bool(
            not accepted
            and inside_real_band
            and near_real
            and real_span_overrun <= config.beyn_near_contour_real_tolerance
        )

        accepted_count += int(accepted)
        near_contour_count += int(near_contour_seed)
        rows.append(
            BeynEigenvalueRow(
                epsilon=epsilon,
                quadrature_points=Nq,
                raw_index=idx,
                real_part=float(eig.real),
                imag_part=float(eig.imag),
                inside_contour=inside,
                inside_real_band=inside_real_band,
                near_real=near_real,
                real_span_overrun=float(real_span_overrun),
                near_contour_seed=near_contour_seed,
                accepted_for_real_refinement=accepted,
            )
        )

    kept_last = float(s[rank - 1]) if rank > 0 else math.nan
    discarded_first = float(s[rank]) if rank < len(s) else math.nan
    diagnostics = BeynDiagnostics(
        epsilon=epsilon,
        M=M,
        quadrature_points=Nq,
        probe_dim=config.beyn_probe_dim,
        threshold_rank=threshold_rank,
        gap_rank=gap_rank,
        estimated_rank=rank,
        selected_gap_ratio=selected_gap,
        leading_s0_singular_value=float(s[0]) if len(s) else math.nan,
        trailing_kept_s0_singular_value=kept_last,
        first_discarded_s0_singular_value=discarded_first,
        max_linear_solve_relative_residual=max_solve_residual,
        raw_eigenvalues=len(raw_eigenvalues),
        near_real_candidates=accepted_count,
        near_contour_seeds=near_contour_count,
    )

    print(
        f"  rank diagnostics: threshold-rank={threshold_rank}, "
        f"gap-rank={gap_rank}, selected gap={selected_gap:.3e}"
    )
    print(f"  Beyn estimated enclosed rank = {rank}")
    if len(s):
        shown = ", ".join(f"{x:.3e}" for x in s[: min(len(s), 8)])
        print(f"  singular values of S0: [{shown}]")
    print(f"  max contour linear-solve residual = {max_solve_residual:.3e}")

    if rank >= config.beyn_probe_dim - 1:
        print(
            "  WARNING: detected Beyn rank is close to probe_dim. Increase "
            "beyn_probe_dim before interpreting the eigenvalue count."
        )

    if len(raw_eigenvalues):
        print("  raw Beyn eigenvalues:")
        for idx, eig in enumerate(raw_eigenvalues):
            row = rows[idx]
            if row.accepted_for_real_refinement:
                tag = "real-candidate"
            elif row.near_contour_seed:
                tag = "near-contour-seed"
            elif not row.inside_contour or not (left < eig.real < right):
                tag = "outside-contour/band"
            elif not row.near_real:
                tag = "complex"
            else:
                tag = "other"
            print(f"    [{idx}] {eig.real:+.12f} {eig.imag:+.3e}i  {tag}")
    else:
        print("  raw Beyn eigenvalues: none")

    return raw_eigenvalues, s, diagnostics, rows, S0


def cluster_real_candidates(
    eigenvalues: np.ndarray,
    rows: list[BeynEigenvalueRow],
    config: Config,
) -> list[complex]:
    accepted = [
        complex(eigenvalues[row.raw_index])
        for row in rows
        if row.accepted_for_real_refinement
    ]
    if not accepted:
        return []

    accepted.sort(key=lambda z: z.real)
    clusters: list[list[complex]] = [[accepted[0]]]
    for z in accepted[1:]:
        if abs(z.real - clusters[-1][-1].real) <= config.beyn_cluster_tolerance:
            clusters[-1].append(z)
        else:
            clusters.append([z])

    return [sum(cluster) / len(cluster) for cluster in clusters]


def cluster_near_contour_seeds(
    eigenvalues: np.ndarray,
    rows: list[BeynEigenvalueRow],
    config: Config,
) -> list[complex]:
    values = [
        complex(eigenvalues[row.raw_index]) for row in rows if row.near_contour_seed
    ]
    if not values:
        return []

    values.sort(key=lambda z: z.real)
    clusters: list[list[complex]] = [[values[0]]]
    for z in values[1:]:
        if abs(z.real - clusters[-1][-1].real) <= config.beyn_cluster_tolerance:
            clusters[-1].append(z)
        else:
            clusters.append([z])
    return [sum(cluster) / len(cluster) for cluster in clusters]


def primary_raw_eigenvalue(eigenvalues: np.ndarray, config: Config) -> complex | None:
    """Return one raw near-real eigenvalue for quadrature-convergence diagnostics."""
    if len(eigenvalues) == 0:
        return None

    values = [complex(z) for z in eigenvalues]
    near_real = [
        z for z in values if abs(z.imag) <= max(config.beyn_real_axis_tolerance, 1.0e-4)
    ]
    pool = near_real if near_real else values
    cutoff = kb_cutoff(config)
    return min(pool, key=lambda z: (abs(z.imag), abs(z.real - cutoff)))


def aitken_delta_squared(values: list[float]) -> float:
    """Aitken Delta^2 extrapolation using the last three scalar iterates."""
    if len(values) < 3:
        return math.nan
    x0, x1, x2 = map(float, values[-3:])
    denominator = x2 - 2.0 * x1 + x0
    scale = max(abs(x0), abs(x1), abs(x2), 1.0)
    if abs(denominator) <= 100.0 * np.finfo(float).eps * scale:
        return math.nan
    return float(x0 - (x1 - x0) ** 2 / denominator)


def aitken_seed_from_history(
    raw_history: list[complex],
    config: Config,
) -> DiscoverySeed | None:
    """
    Build a fallback local-SVD seed from the raw Beyn sequence.

    This is intentionally *not* an eigenvalue acceptance criterion.  It is only
    used after all adaptive contour quadrature levels have been exhausted and a
    stable rank-one Beyn count has failed to produce a strict candidate.
    """
    if not config.beyn_enable_aitken_seed or len(raw_history) < 3:
        return None

    last = raw_history[-3:]
    reals = [z.real for z in last]
    if config.beyn_aitken_require_contracting_shifts:
        d1 = abs(reals[1] - reals[0])
        d2 = abs(reals[2] - reals[1])
        if not (d2 < d1):
            return None

    estimate = aitken_delta_squared(reals)
    if not np.isfinite(estimate):
        return None

    left, right, _, _, _ = contour_geometry(config)
    cutoff = kb_cutoff(config)
    if not (0.0 < estimate < cutoff):
        return None

    # Prefer an estimate inside the intended contour real span.  A tiny overrun
    # is allowed only because this value is a seed that must still pass SVD.
    overrun = max(left - estimate, estimate - right, 0.0)
    if overrun > config.beyn_near_contour_real_tolerance:
        return None

    return DiscoverySeed(real=float(estimate), imag=0.0, source="aitken")


def candidate_set_shift(previous: list[complex], current: list[complex]) -> float:
    """Maximum real-part shift for equally-sized sorted candidate sets."""
    if len(previous) != len(current) or not current:
        return math.nan
    p = sorted(previous, key=lambda z: z.real)
    c = sorted(current, key=lambda z: z.real)
    return float(max(abs(z1.real - z0.real) for z0, z1 in zip(p, c, strict=True)))


def run_beyn_convergence_study(
    epsilon: float,
    config: Config,
) -> tuple[
    np.ndarray,
    np.ndarray,
    BeynDiagnostics,
    list[BeynEigenvalueRow],
    list[complex],
    list[DiscoverySeed],
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
    float,
]:
    """
    Run the base Beyn quadrature study and adaptively increase Nq only when
    discovery still needs it.

    The base levels are always run.  After that, expensive adaptive levels are
    used only if (a) no strict real in-contour candidate exists or (b) the Beyn
    rank has not stabilized.  If all levels are exhausted without a strict
    candidate, a near-contour raw estimate or Aitken Delta^2 estimate may be
    passed to the local SVD stage strictly as a *seed*.
    """
    all_diagnostics: list[BeynDiagnostics] = []
    all_eigen_rows: list[BeynEigenvalueRow] = []
    convergence_rows: list[BeynConvergenceRow] = []

    previous_S0: np.ndarray | None = None
    previous_candidates: list[complex] = []
    previous_rank: int | None = None
    raw_primary_history: list[complex] = []

    final_raw = np.array([], dtype=np.complex128)
    final_s = np.array([], dtype=float)
    final_diag: BeynDiagnostics | None = None
    final_rows: list[BeynEigenvalueRow] = []
    final_candidates: list[complex] = []

    base_levels = list(config.beyn_quadrature_levels)
    adaptive_levels = [
        int(nq)
        for nq in config.beyn_adaptive_quadrature_levels
        if int(nq) > max(base_levels, default=0)
    ]
    levels = base_levels + adaptive_levels
    base_count = len(base_levels)

    print("  Beyn quadrature-convergence study:")

    for level_index, Nq in enumerate(levels):
        if level_index >= base_count:
            print(f"  adaptive Beyn escalation -> Nq={Nq}")

        raw, svals, diag, rows, S0 = beyn_discover(epsilon, Nq, config)
        candidates = cluster_real_candidates(raw, rows, config)
        near_seeds = cluster_near_contour_seeds(raw, rows, config)
        primary_raw = primary_raw_eigenvalue(raw, config)
        if primary_raw is not None:
            raw_primary_history.append(primary_raw)

        if previous_S0 is None:
            s0_change = math.nan
        else:
            s0_change = float(
                np.linalg.norm(S0 - previous_S0) / max(np.linalg.norm(S0), 1.0e-30)
            )

        shift = candidate_set_shift(previous_candidates, candidates)
        rank_stable = bool(
            previous_rank is not None and previous_rank == diag.estimated_rank
        )

        primary = (
            min(candidates, key=lambda z: abs(z.imag))
            if candidates
            else complex(math.nan, math.nan)
        )
        raw_for_row = (
            primary_raw if primary_raw is not None else complex(math.nan, math.nan)
        )

        convergence_rows.append(
            BeynConvergenceRow(
                epsilon=epsilon,
                quadrature_points=Nq,
                estimated_rank=diag.estimated_rank,
                threshold_rank=diag.threshold_rank,
                gap_rank=diag.gap_rank,
                selected_gap_ratio=diag.selected_gap_ratio,
                s0_sv1=float(svals[0]) if len(svals) > 0 else math.nan,
                s0_sv2=float(svals[1]) if len(svals) > 1 else math.nan,
                s0_sv3=float(svals[2]) if len(svals) > 2 else math.nan,
                s0_relative_change_from_previous=s0_change,
                accepted_candidates=len(candidates),
                near_contour_seeds=len(near_seeds),
                primary_candidate_real=float(primary.real),
                primary_candidate_imag=float(primary.imag),
                primary_raw_real=float(raw_for_row.real),
                primary_raw_imag=float(raw_for_row.imag),
                max_candidate_shift_from_previous=shift,
                rank_stable_from_previous=rank_stable,
            )
        )

        change_text = f"{s0_change:.3e}" if np.isfinite(s0_change) else "--"
        shift_text = f"{shift:.3e}" if np.isfinite(shift) else "--"
        raw_text = (
            f"{raw_for_row.real:.12f}{raw_for_row.imag:+.2e}i"
            if np.isfinite(raw_for_row.real)
            else "--"
        )
        print(
            f"    Nq={Nq}: rank={diag.estimated_rank}, "
            f"strict={len(candidates)}, near-contour={len(near_seeds)}, "
            f"S0-change={change_text}, candidate-shift={shift_text}, raw={raw_text}"
        )

        all_diagnostics.append(diag)
        all_eigen_rows.extend(rows)

        previous_S0 = S0
        previous_candidates = candidates
        previous_rank = diag.estimated_rank

        final_raw = raw
        final_s = svals
        final_diag = diag
        final_rows = rows
        final_candidates = candidates

        # Always complete the configured base convergence study first.
        if level_index + 1 < base_count:
            continue

        # After the base study, stop as soon as global discovery is usable:
        # a strict candidate exists and the enclosed rank is stable.
        if candidates and rank_stable:
            break

        # Otherwise continue through the adaptive levels, if any remain.
        if level_index + 1 < len(levels):
            reason = []
            if not candidates:
                reason.append("no strict candidate")
            if not rank_stable:
                reason.append("rank not yet stable")
            print("    continuing adaptive quadrature: " + ", ".join(reason))

    assert final_diag is not None

    # Aitken is recorded as a diagnostic even when a strict candidate exists.
    aitken_estimate = aitken_delta_squared([z.real for z in raw_primary_history])

    if final_candidates:
        local_seeds = [
            DiscoverySeed(real=float(z.real), imag=float(z.imag), source="strict-beyn")
            for z in final_candidates
        ]
    else:
        local_seeds: list[DiscoverySeed] = []

        # Fallbacks are allowed only when the global count itself is stable and
        # indicates one enclosed mode.  They merely initialize local SVD.
        final_rank_stable = bool(
            len(convergence_rows) >= 2
            and convergence_rows[-1].estimated_rank
            == convergence_rows[-2].estimated_rank
        )
        if final_diag.estimated_rank == 1 and final_rank_stable:
            aitken_seed = aitken_seed_from_history(raw_primary_history, config)
            if aitken_seed is not None:
                local_seeds = [aitken_seed]
                print(
                    "  no strict final candidate; using Aitken Delta^2 only as "
                    f"local-SVD seed: kb={aitken_seed.real:.12f}"
                )
            else:
                near_values = cluster_near_contour_seeds(final_raw, final_rows, config)
                if near_values:
                    local_seeds = [
                        DiscoverySeed(
                            real=float(z.real),
                            imag=float(z.imag),
                            source="near-contour",
                        )
                        for z in near_values
                    ]
                    print(
                        "  no strict final candidate; using near-contour raw Beyn "
                        "estimate(s) only as local-SVD seed(s)."
                    )

    local_seeds.sort(key=lambda seed: seed.real)
    return (
        final_raw,
        final_s,
        final_diag,
        final_rows,
        final_candidates,
        local_seeds,
        all_diagnostics,
        all_eigen_rows,
        convergence_rows,
        float(aitken_estimate),
    )


# ---------------------------------------------------------------------------
# Local real-axis verification of Beyn candidates
# ---------------------------------------------------------------------------


def candidate_brackets(
    seeds: list[DiscoverySeed],
    config: Config,
) -> list[tuple[float, float]]:
    """Build non-overlapping real-axis SVD brackets around discovery seeds."""
    if not seeds:
        return []

    left_band = config.beyn_low_k_margin
    right_band = kb_cutoff(config) - config.beyn_cutoff_margin
    tiny = max(1.0e-12, 1.0e-8 * (right_band - left_band))

    xs_raw = np.array([seed.real for seed in seeds], dtype=float)
    # A near-contour fallback can sit microscopically outside the contour real
    # span because of quadrature error. Clamp only the *bracket center*; the raw
    # value remains reported unchanged and SVD must certify the actual minimum.
    xs = np.clip(xs_raw, left_band + tiny, right_band - tiny)

    brackets: list[tuple[float, float]] = []
    for i, x in enumerate(xs):
        neighbor_left = 0.5 * (xs[i - 1] + x) if i > 0 else left_band
        neighbor_right = 0.5 * (x + xs[i + 1]) if i + 1 < len(xs) else right_band

        left = max(left_band, neighbor_left, x - config.local_refine_half_width)
        right = min(right_band, neighbor_right, x + config.local_refine_half_width)

        if not left < right:
            raise RuntimeError(
                f"Could not build local bracket around discovery seed kb={xs_raw[i]}."
            )
        brackets.append((float(left), float(right)))

    return brackets


def refine_candidate_for_M(
    epsilon: float,
    M: int,
    left: float,
    right: float,
    config: Config,
) -> tuple[float, float, float, float, float, float, float, bool]:
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
            f"Local SVD refinement failed for epsilon={epsilon}, M={M}: {result.message}"
        )

    kb = float(result.x)
    sigma_min, sigma_max, relative_sv = singular_metrics(kb, epsilon, M, config)
    drop = min(left_value, right_value) / max(sigma_min, tiny)

    width = right - left
    edge_margin = 0.005 * width
    interior = (kb > left + edge_margin) and (kb < right - edge_margin)
    return (
        kb,
        sigma_min,
        sigma_max,
        relative_sv,
        left_value,
        right_value,
        drop,
        interior,
    )


def run_candidate_refinement(
    epsilon: float,
    candidate_index: int,
    seed: DiscoverySeed,
    bracket: tuple[float, float],
    config: Config,
) -> tuple[list[ModeRefinementRow], ModeResult]:
    rows: list[ModeRefinementRow] = []
    previous_sigma = math.nan
    left, right = bracket

    print(
        f"  candidate {candidate_index} [{seed.source}]: seed kb={seed.real:.12f} "
        f"{seed.imag:+.3e}i, local bracket=[{left:.12f}, {right:.12f}]"
    )

    for M in config.refinement_M:
        (
            kb,
            sigma_min,
            sigma_max,
            relative_sv,
            left_value,
            right_value,
            drop,
            interior,
        ) = refine_candidate_for_M(epsilon, M, left, right, config)
        sigma_bem = sigma_from_kb(kb, config)
        change = (
            abs(sigma_bem - previous_sigma) / max(abs(sigma_bem), 1.0e-30)
            if np.isfinite(previous_sigma)
            else math.nan
        )

        row = ModeRefinementRow(
            epsilon=epsilon,
            candidate_index=candidate_index,
            seed_source=seed.source,
            beyn_kb_real=float(seed.real),
            beyn_kb_imag=float(seed.imag),
            M=M,
            kb=kb,
            sigma_bem=sigma_bem,
            sigma_min=sigma_min,
            sigma_max=sigma_max,
            relative_singular_value=relative_sv,
            left_value=left_value,
            right_value=right_value,
            drop_factor=drop,
            minimum_is_interior=interior,
            relative_sigma_change_from_previous=change,
        )
        rows.append(row)
        previous_sigma = sigma_bem

        change_text = f"{change:.3%}" if np.isfinite(change) else "--"
        print(
            f"    M={M:>2}: kb={kb:.12f}, sigma_BEM={sigma_bem:.8e}, "
            f"sv_min={sigma_min:.3e}, sv_min/sv_max={relative_sv:.3e}, "
            f"drop={drop:.2e}, mesh_change={change_text}, "
            f"interior={'yes' if interior else 'no'}"
        )

    final = rows[-1]
    previous = rows[-2] if len(rows) >= 2 else rows[-1]
    final_change = abs(final.sigma_bem - previous.sigma_bem) / max(
        abs(final.sigma_bem), 1.0e-30
    )

    # Certification uses the scale-invariant singular-value ratio rather than
    # an absolute sigma_min threshold.  The absolute value is still reported.
    resolved = bool(
        final.minimum_is_interior
        and final.relative_singular_value <= config.relative_near_singular_tolerance
        and final.drop_factor >= config.minimum_drop_factor
        and final_change <= config.mesh_sigma_relative_tolerance
    )

    result = ModeResult(
        epsilon=epsilon,
        candidate_index=candidate_index,
        seed_source=seed.source,
        beyn_kb_real=float(seed.real),
        beyn_kb_imag=float(seed.imag),
        kb_numerical=final.kb,
        sigma_numerical=final.sigma_bem,
        sigma_min_final=final.sigma_min,
        sigma_max_final=final.sigma_max,
        relative_singular_value_final=final.relative_singular_value,
        final_drop_factor=final.drop_factor,
        final_relative_mesh_change=final_change,
        resolved=resolved,
    )
    print(f"    resolved={'YES' if resolved else 'no'}")
    return rows, result


# ---------------------------------------------------------------------------
# Independent physical trapped-mode diagnostics
# ---------------------------------------------------------------------------


def make_extended_field_green(
    kb: float,
    config: Config,
) -> Callable[..., complex]:
    """
    Dirichlet waveguide Green function for off-boundary field reconstruction.

    The lattice-sum representation used by greens_periodic requires the polar
    radius to stay below one period.  For beta=0 the periodic Green function is
    exactly periodic in Y, so we wrap each image separation into the nearest
    periodic representative before evaluating it.  This extends the usable x
    range for the independent decay diagnostics without changing the Green
    function mathematically.
    """
    b = config.b
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
        y_field = y + b + config.a
        y_source = eta + b + config.a
        X = x - xi
        Y1 = wrap_periodic_y(y_field - y_source)
        Y2 = wrap_periodic_y(y_field + y_source)
        term1 = lattice.greens_periodic(X, Y1, coefficients, k, period)
        term2 = lattice.greens_periodic(X, Y2, coefficients, k, period)
        return term1 - term2

    return green


def weighted_kernel_at_field_point(
    x: float,
    y: float,
    theta: float,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
) -> complex:
    """Weighted double-layer kernel at a field point away from the boundary."""
    xi, eta, xi_p, eta_p, _, _ = circle_geometry(theta, epsilon)
    G_xi, G_eta = source_derivatives(
        G,
        x,
        y,
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
    """Reconstruct u(P) = integral_gamma u(q) dG/dn_q ds_q off the boundary."""
    delta_theta = 2.0 * PI / len(theta)
    values = np.empty(len(points), dtype=np.complex128)
    for i, (x, y) in enumerate(points):
        kernel = np.array(
            [
                weighted_kernel_at_field_point(
                    float(x), float(y), float(t), epsilon, config, G
                )
                for t in theta
            ],
            dtype=np.complex128,
        )
        values[i] = delta_theta * np.dot(kernel, boundary_vector)
    return values


def cross_section_field(
    x: float,
    y_values: np.ndarray,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    epsilon: float,
    config: Config,
    G: Callable[..., complex],
) -> np.ndarray:
    points = [(float(x), float(y)) for y in y_values]
    return reconstruct_field_at_points(
        points, boundary_vector, theta, epsilon, config, G
    )


def physical_mode_diagnostics(
    epsilon: float,
    kb: float,
    M: int,
    config: Config,
) -> PhysicalDiagnostics:
    """
    Independently inspect the reconstructed trapped mode.

    Checks:
      * Dirichlet walls y=+-b (in obstacle-centred coordinates y=-b-a,b-a),
      * approximate homogeneous Neumann data on the obstacle,
      * exponential L2 decay of cross sections as |x| grows.

    The Neumann test is deliberately a diagnostic rather than part of spectral
    certification because evaluating a double-layer field very close to its
    boundary is numerically delicate.
    """
    _, _, _, boundary_vector = smallest_singular_pair(kb, epsilon, M, config)
    theta = boundary_nodes(M)
    G = make_extended_field_green(kb, config)

    b = config.b
    y_lower = -b - config.a
    y_upper = b - config.a

    # Gauss-Legendre rule for cross-section L2 norms.
    nodes, weights = leggauss(config.physical_quadrature_points)
    y_values = 0.5 * (y_upper - y_lower) * nodes + 0.5 * (y_upper + y_lower)
    y_weights = 0.5 * (y_upper - y_lower) * weights

    distances = np.array(config.physical_x_over_b, dtype=float) * b
    left_norms: list[float] = []
    right_norms: list[float] = []
    cross_section_peak = 0.0
    for distance in distances:
        field_left = cross_section_field(
            -float(distance), y_values, boundary_vector, theta, epsilon, config, G
        )
        field_right = cross_section_field(
            float(distance), y_values, boundary_vector, theta, epsilon, config, G
        )
        left_norm = float(math.sqrt(np.sum(y_weights * np.abs(field_left) ** 2)))
        right_norm = float(math.sqrt(np.sum(y_weights * np.abs(field_right) ** 2)))
        left_norms.append(left_norm)
        right_norms.append(right_norm)
        cross_section_peak = max(
            cross_section_peak,
            float(np.max(np.abs(field_left))),
            float(np.max(np.abs(field_right))),
        )

    left_arr = np.asarray(left_norms, dtype=float)
    right_arr = np.asarray(right_norms, dtype=float)
    monotone_left = bool(np.all(np.diff(left_arr) <= 1.0e-10 * max(left_arr[0], 1.0)))
    monotone_right = bool(
        np.all(np.diff(right_arr) <= 1.0e-10 * max(right_arr[0], 1.0))
    )

    def fit_decay(norms: np.ndarray) -> float:
        positive = np.maximum(norms, np.finfo(float).tiny)
        slope, _ = np.polyfit(distances, np.log(positive), 1)
        return float(-slope)

    decay_left = fit_decay(left_arr)
    decay_right = fit_decay(right_arr)
    expected_decay = sigma_from_kb(kb, config)
    decay_error_left = abs(decay_left - expected_decay) / max(expected_decay, 1.0e-30)
    decay_error_right = abs(decay_right - expected_decay) / max(expected_decay, 1.0e-30)

    # Wall values at three longitudinal locations.  Because the waveguide Green
    # function already satisfies Dirichlet data, these should be near roundoff.
    x0 = float(distances[0]) if len(distances) else 2.0 * b
    wall_points = [
        (-x0, y_lower), (-x0, y_upper),
        (0.0, y_lower), (0.0, y_upper),
        (x0, y_lower), (x0, y_upper),
    ]
    wall_values = reconstruct_field_at_points(
        wall_points, boundary_vector, theta, epsilon, config, G
    )
    amplitude_scale = max(
        float(np.max(np.abs(boundary_vector))), cross_section_peak, 1.0e-30
    )
    wall_relative = float(np.max(np.abs(wall_values)) / amplitude_scale)

    # One-sided second-order normal derivative from the fluid side of a circle.
    sample_t = np.linspace(
        0.0, 2.0 * PI, config.physical_neumann_samples, endpoint=False
    )
    offset = max(
        config.physical_neumann_offset_fraction_of_epsilon * epsilon,
        2.0e-5 * b,
    )
    periodic_theta = np.r_[theta[-1] - 2.0 * PI, theta, theta[0] + 2.0 * PI]
    periodic_values = np.r_[boundary_vector[-1], boundary_vector, boundary_vector[0]]
    boundary_samples = np.interp(sample_t, periodic_theta, periodic_values)
    normal_derivatives: list[complex] = []
    for t, u0 in zip(sample_t, boundary_samples, strict=True):
        x0q = epsilon * math.cos(float(t))
        y0q = epsilon * math.sin(float(t))
        nx = math.cos(float(t))
        ny = math.sin(float(t))
        points = [
            (x0q + offset * nx, y0q + offset * ny),
            (x0q + 2.0 * offset * nx, y0q + 2.0 * offset * ny),
        ]
        u1, u2 = reconstruct_field_at_points(
            points, boundary_vector, theta, epsilon, config, G
        )
        derivative = (-3.0 * u0 + 4.0 * u1 - u2) / (2.0 * offset)
        normal_derivatives.append(complex(derivative))

    neumann_scale = max(float(np.max(np.abs(boundary_vector))), 1.0e-30)
    neumann_relative = float(
        epsilon * np.max(np.abs(np.asarray(normal_derivatives))) / neumann_scale
    )

    walls_ok = wall_relative <= config.physical_wall_relative_tolerance
    neumann_ok = neumann_relative <= config.physical_neumann_relative_tolerance
    decay_ok = bool(
        monotone_left
        and monotone_right
        and decay_error_left <= config.physical_decay_relative_tolerance
        and decay_error_right <= config.physical_decay_relative_tolerance
    )

    return PhysicalDiagnostics(
        epsilon=epsilon,
        kb=kb,
        M=M,
        wall_relative_residual=wall_relative,
        neumann_relative_residual=neumann_relative,
        decay_rate_left=decay_left,
        decay_rate_right=decay_right,
        expected_decay_rate=expected_decay,
        decay_relative_error_left=float(decay_error_left),
        decay_relative_error_right=float(decay_error_right),
        monotone_decay_left=monotone_left,
        monotone_decay_right=monotone_right,
        walls_verified=walls_ok,
        neumann_verified=neumann_ok,
        decay_verified=decay_ok,
    )


# ---------------------------------------------------------------------------
# Validation / reporting
# ---------------------------------------------------------------------------


def validate_epsilon(
    epsilon: float,
    base_config: Config,
) -> tuple[
    ValidationResult,
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
    list[ModeRefinementRow],
    list[ModeResult],
    np.ndarray,
    list[BeynEigenvalueRow],
    float,
    PhysicalDiagnostics | None,
]:
    """Validate one epsilon using an epsilon-adapted near-cutoff contour."""
    config = config_for_epsilon(epsilon, base_config)
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
    coefficient = asymptotic_coefficient(config)

    (
        final_raw,
        final_s,
        final_diag,
        final_eigen_rows,
        strict_candidates,
        local_seeds,
        all_beyn_diagnostics,
        all_eigen_rows,
        convergence_rows,
        aitken_estimate,
    ) = run_beyn_convergence_study(epsilon, config)

    brackets = candidate_brackets(local_seeds, config)
    print(f"  effective cutoff margin               = {config.beyn_cutoff_margin:.3e}")
    print(f"  final strict near-real Beyn candidates = {len(strict_candidates)}")
    print(f"  local refinement seeds                 = {len(local_seeds)}")

    all_refinement: list[ModeRefinementRow] = []
    mode_results: list[ModeResult] = []
    for i, (seed, bracket) in enumerate(
        zip(local_seeds, brackets, strict=True), start=1
    ):
        rows, mode = run_candidate_refinement(epsilon, i, seed, bracket, config)
        all_refinement.extend(rows)
        mode_results.append(mode)

    resolved = [mode for mode in mode_results if mode.resolved]

    rank_stable = bool(
        len(convergence_rows) >= 2
        and convergence_rows[-1].estimated_rank == convergence_rows[-2].estimated_rank
    )
    candidate_stable = bool(
        len(convergence_rows) >= 2
        and np.isfinite(convergence_rows[-1].max_candidate_shift_from_previous)
        and convergence_rows[-1].max_candidate_shift_from_previous
        <= config.beyn_candidate_convergence_tolerance
    )

    # Independent contour-geometry check: move the right endpoint closer to the
    # cutoff and verify that the enclosed rank is unchanged.  This specifically
    # guards against accidentally truncating a mode that is O(epsilon^4) from
    # Lambda_1.
    tighter_margin_consistent: bool | None = None
    if config.beyn_check_tighter_cutoff_margin:
        tighter_margin = max(
            config.beyn_min_cutoff_margin,
            config.beyn_cutoff_margin * config.beyn_tighter_cutoff_margin_factor,
        )
        if tighter_margin < config.beyn_cutoff_margin * (1.0 - 1.0e-12):
            tight_config = replace(config, beyn_cutoff_margin=tighter_margin)
            print(
                f"  tighter-cutoff diagnostic: margin {config.beyn_cutoff_margin:.3e} "
                f"-> {tighter_margin:.3e}"
            )
            _, _, tight_diag, _, _ = beyn_discover(
                epsilon, final_diag.quadrature_points, tight_config
            )
            tighter_margin_consistent = bool(
                tight_diag.estimated_rank == final_diag.estimated_rank
            )
            print(
                "  tighter-cutoff rank consistency = "
                f"{'PASS' if tighter_margin_consistent else 'FAIL'} "
                f"({final_diag.estimated_rank} -> {tight_diag.estimated_rank})"
            )
        else:
            tighter_margin_consistent = True

    margin_ok = tighter_margin_consistent is not False
    one_mode_count_supported = bool(
        final_diag.estimated_rank == 1
        and rank_stable
        and len(resolved) == 1
        and margin_ok
    )

    # IMPORTANT: no theory-driven selection among several numerical modes.
    # If more than one mode is resolved, uniqueness is not supported and the
    # asymptotic comparison is intentionally left N/A.
    physical: PhysicalDiagnostics | None = None
    if len(resolved) == 1:
        mode = resolved[0]
        kb_num = mode.kb_numerical
        sigma_num = mode.sigma_numerical
        sv_final = mode.sigma_min_final
        sv_max_final = mode.sigma_max_final
        relative_sv_final = mode.relative_singular_value_final
        drop_final = mode.final_drop_factor
        mesh_change = mode.final_relative_mesh_change
        error_kb = abs(kb_num - kb_asym) / abs(kb_asym)
        error_sigma = abs(sigma_num - sigma_asym) / abs(sigma_asym)
        asymptotic_ok: bool | None = (
            error_sigma <= config.relative_sigma_error_tolerance
        )
        sigma_over_eps2 = sigma_num / epsilon**2
        remainder_denom = epsilon**3 * abs(math.log(epsilon))
        scaled_remainder = abs(sigma_num - coefficient * epsilon**2) / max(
            remainder_denom, 1.0e-30
        )

        if config.run_physical_diagnostics:
            print("  reconstructing field for independent physical diagnostics...")
            physical = physical_mode_diagnostics(
                epsilon, kb_num, config.refinement_M[-1], config
            )
            print(
                f"    wall residual={physical.wall_relative_residual:.3e}, "
                f"Neumann residual={physical.neumann_relative_residual:.3e}"
            )
            print(
                f"    decay rates left/right="
                f"{physical.decay_rate_left:.6e}/{physical.decay_rate_right:.6e}, "
                f"expected sigma={physical.expected_decay_rate:.6e}"
            )
    else:
        kb_num = math.nan
        sigma_num = math.nan
        sv_final = math.nan
        sv_max_final = math.nan
        relative_sv_final = math.nan
        drop_final = math.nan
        mesh_change = math.nan
        error_kb = math.nan
        error_sigma = math.nan
        sigma_over_eps2 = math.nan
        scaled_remainder = math.nan
        asymptotic_ok = None

    if not candidate_stable and strict_candidates:
        print(
            "  NOTE: final strict Beyn candidate position is still moving by more "
            f"than {config.beyn_candidate_convergence_tolerance:.1e} between the "
            "last two Nq levels. Local SVD refinement remains the trusted position."
        )

    seed_source = (
        ",".join(sorted({mode.seed_source for mode in mode_results})) or "none"
    )
    result = ValidationResult(
        epsilon=epsilon,
        a=config.a,
        lambda_1=lambda_1(config),
        kb_cutoff=kb_cutoff(config),
        effective_beyn_cutoff_margin=config.beyn_cutoff_margin,
        tighter_margin_rank_consistent=tighter_margin_consistent,
        kb_asymptotic=kb_asym,
        sigma_asymptotic=sigma_asym,
        asymptotic_coefficient=coefficient,
        beyn_final_quadrature_points=final_diag.quadrature_points,
        beyn_estimated_rank=final_diag.estimated_rank,
        beyn_rank_stable=rank_stable,
        beyn_near_real_candidates=len(strict_candidates),
        local_refinement_seed_count=len(local_seeds),
        local_refinement_seed_source=seed_source,
        aitken_estimate=aitken_estimate,
        resolved_mode_count=len(resolved),
        kb_numerical=kb_num,
        sigma_numerical=sigma_num,
        sigma_over_epsilon_squared=sigma_over_eps2,
        scaled_asymptotic_remainder=scaled_remainder,
        sigma_min_final=sv_final,
        sigma_max_final=sv_max_final,
        relative_singular_value_final=relative_sv_final,
        final_drop_factor=drop_final,
        final_relative_mesh_change=mesh_change,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        unique_mode_verified=one_mode_count_supported,
        asymptotic_agreement_verified=asymptotic_ok,
        wall_relative_residual=(
            physical.wall_relative_residual if physical is not None else math.nan
        ),
        neumann_relative_residual=(
            physical.neumann_relative_residual if physical is not None else math.nan
        ),
        decay_rate_left=(physical.decay_rate_left if physical is not None else math.nan),
        decay_rate_right=(physical.decay_rate_right if physical is not None else math.nan),
        physical_decay_verified=(physical.decay_verified if physical is not None else None),
    )

    return (
        result,
        all_beyn_diagnostics,
        all_eigen_rows,
        convergence_rows,
        all_refinement,
        mode_results,
        final_s,
        final_eigen_rows,
        aitken_estimate,
        physical,
    )


def write_dataclass_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def run_internal_convergence_study(
    results: list[ValidationResult],
    base_config: Config,
) -> list[InternalConvergenceRow]:
    """One-at-a-time sensitivity to h, lattice truncation, and harmonic order."""
    if not base_config.run_internal_convergence_study:
        return []

    finite = [row for row in results if np.isfinite(row.kb_numerical)]
    if not finite:
        return []

    target = min(
        finite,
        key=lambda row: abs(row.epsilon - base_config.internal_convergence_epsilon),
    )
    epsilon = target.epsilon
    baseline_kb = target.kb_numerical
    baseline_sigma = target.sigma_numerical
    config = config_for_epsilon(epsilon, base_config)
    M = base_config.internal_convergence_M
    half = base_config.internal_refine_half_width
    left = max(config.beyn_low_k_margin, baseline_kb - half)
    right = min(kb_cutoff(config) - config.beyn_cutoff_margin, baseline_kb + half)
    if not left < right:
        return []

    variants: list[tuple[str, float, Config]] = []
    for h in base_config.finite_difference_steps_test:
        variants.append(("finite_difference_step", float(h), replace(config, finite_difference_step=float(h))))
    for terms in base_config.lattice_terms_test:
        variants.append(("lattice_terms", float(terms), replace(config, lattice_terms=int(terms))))
    for order in base_config.harmonic_orders_test:
        variants.append(("harmonic_order", float(order), replace(config, harmonic_order=int(order))))

    rows: list[InternalConvergenceRow] = []
    print(f"\n=== INTERNAL GREEN / DERIVATIVE CONVERGENCE at epsilon={epsilon:.3f} ===")
    for parameter, value, variant_config in variants:
        (
            kb,
            _,
            _,
            relative_sv,
            _,
            _,
            _,
            _,
        ) = refine_candidate_for_M(epsilon, M, left, right, variant_config)
        sigma = sigma_from_kb(kb, variant_config)
        kb_shift = abs(kb - baseline_kb) / max(abs(baseline_kb), 1.0e-30)
        sigma_shift = abs(sigma - baseline_sigma) / max(abs(baseline_sigma), 1.0e-30)
        rows.append(
            InternalConvergenceRow(
                epsilon=epsilon,
                parameter=parameter,
                value=value,
                kb=kb,
                sigma_bem=sigma,
                relative_singular_value=relative_sv,
                relative_kb_shift_from_baseline=kb_shift,
                relative_sigma_shift_from_baseline=sigma_shift,
            )
        )
        print(
            f"  {parameter}={value:g}: kb={kb:.12f}, sigma={sigma:.8e}, "
            f"rel_sv={relative_sv:.3e}, delta_sigma={sigma_shift:.3%}"
        )
    return rows


def classify_mode_count_for_geometry(
    epsilon: float,
    a: float,
    base_config: Config,
) -> tuple[str, int, int]:
    """Return ('zero'|'one'|'ambiguous', Beyn rank, resolved count)."""
    trial = replace(base_config, a=float(a), run_physical_diagnostics=False)
    config = config_for_epsilon(epsilon, trial)
    (
        _,
        _,
        final_diag,
        _,
        _,
        local_seeds,
        _,
        _,
        convergence_rows,
        _,
    ) = run_beyn_convergence_study(epsilon, config)

    brackets = candidate_brackets(local_seeds, config)
    resolved_count = 0
    for i, (seed, bracket) in enumerate(zip(local_seeds, brackets, strict=True), start=1):
        _, mode = run_candidate_refinement(epsilon, i, seed, bracket, config)
        resolved_count += int(mode.resolved)

    rank_stable = bool(
        len(convergence_rows) >= 2
        and convergence_rows[-1].estimated_rank == convergence_rows[-2].estimated_rank
    )

    # Near-cutoff safeguard for the critical-height test.
    tighter_consistent = True
    tighter_margin = max(
        config.beyn_min_cutoff_margin,
        config.beyn_cutoff_margin * config.beyn_tighter_cutoff_margin_factor,
    )
    if tighter_margin < config.beyn_cutoff_margin * (1.0 - 1.0e-12):
        tight_config = replace(config, beyn_cutoff_margin=tighter_margin)
        _, _, tight_diag, _, _ = beyn_discover(
            epsilon, final_diag.quadrature_points, tight_config
        )
        tighter_consistent = tight_diag.estimated_rank == final_diag.estimated_rank

    if rank_stable and tighter_consistent and final_diag.estimated_rank == 0 and resolved_count == 0:
        return "zero", final_diag.estimated_rank, resolved_count
    if rank_stable and tighter_consistent and final_diag.estimated_rank == 1 and resolved_count == 1:
        return "one", final_diag.estimated_rank, resolved_count
    return "ambiguous", final_diag.estimated_rank, resolved_count


def run_critical_height_study(base_config: Config) -> list[CriticalHeightRow]:
    """Numerically bracket a_c(epsilon) and test a_c(epsilon)->a0*."""
    if not base_config.run_critical_height_study:
        return []

    a0 = critical_height_leading_order(base_config)
    rows: list[CriticalHeightRow] = []
    cache: dict[tuple[float, float], tuple[str, int, int]] = {}

    def classify(epsilon: float, a: float) -> tuple[str, int, int]:
        key = (round(float(epsilon), 12), round(float(a), 12))
        if key not in cache:
            print(f"\n  critical-height probe: epsilon={epsilon:.4f}, a={a:.8f}")
            cache[key] = classify_mode_count_for_geometry(epsilon, a, base_config)
            print(f"    classification={cache[key][0]}, rank={cache[key][1]}, resolved={cache[key][2]}")
        return cache[key]

    print("\n=== CRITICAL-HEIGHT STUDY a*(epsilon) ===")
    for epsilon in base_config.critical_height_epsilon_values:
        geom_limit = base_config.b - epsilon - 1.0e-5
        lower = max(-geom_limit, a0 - base_config.critical_height_half_width)
        upper = min(geom_limit, a0 + base_config.critical_height_half_width)
        lower_status, _, _ = classify(epsilon, lower)
        upper_status, _, _ = classify(epsilon, upper)

        if lower_status != "zero" or upper_status != "one":
            rows.append(
                CriticalHeightRow(
                    epsilon=epsilon,
                    a_lower_zero_mode=lower if lower_status == "zero" else math.nan,
                    a_upper_one_mode=upper if upper_status == "one" else math.nan,
                    a_critical_estimate=math.nan,
                    bracket_width=math.nan,
                    a0_star=a0,
                    normalized_shift_over_epsilon=math.nan,
                    status="unbracketed-or-ambiguous",
                )
            )
            continue

        for _ in range(base_config.critical_height_bisection_iterations):
            mid = 0.5 * (lower + upper)
            status, _, _ = classify(epsilon, mid)
            if status == "zero":
                lower = mid
            elif status == "one":
                upper = mid
            else:
                break

        estimate = 0.5 * (lower + upper)
        width = upper - lower
        rows.append(
            CriticalHeightRow(
                epsilon=epsilon,
                a_lower_zero_mode=lower,
                a_upper_one_mode=upper,
                a_critical_estimate=estimate,
                bracket_width=width,
                a0_star=a0,
                normalized_shift_over_epsilon=(estimate - a0) / epsilon,
                status="bracketed",
            )
        )
        print(
            f"  epsilon={epsilon:.4f}: a_c~{estimate:.10f}, "
            f"width={width:.3e}, (a_c-a0*)/epsilon={(estimate-a0)/epsilon:.6e}"
        )
    return rows


def plot_beyn_contour(
    epsilon: float,
    eigen_rows: list[BeynEigenvalueRow],
    output_directory: Path,
    config: Config,
) -> None:
    theta = np.linspace(0.0, 2.0 * PI, 500)
    contour = np.array([ellipse_point(float(t), config)[0] for t in theta])

    plt.figure(figsize=(8, 5))
    plt.plot(contour.real, contour.imag, label="Beyn contour")
    if eigen_rows:
        raw = np.array([complex(row.real_part, row.imag_part) for row in eigen_rows])
        plt.scatter(raw.real, raw.imag, marker="x", label="final raw Beyn eigenvalues")
    plt.axhline(0.0, linewidth=0.8)
    plt.axvline(kb_cutoff(config), linestyle="--", label=r"$\sqrt{\Lambda_1}b$")
    plt.xlabel(r"$\operatorname{Re}(kb)$")
    plt.ylabel(r"$\operatorname{Im}(kb)$")
    plt.title(rf"Beyn contour spectrum, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"beyn_contour_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_s0_singular_values(
    epsilon: float,
    singular_values: np.ndarray,
    output_directory: Path,
) -> None:
    if len(singular_values) == 0:
        return
    indices = np.arange(1, len(singular_values) + 1)
    plt.figure(figsize=(7, 4.5))
    plt.semilogy(indices, singular_values, "o-")
    plt.xlabel("index")
    plt.ylabel(r"singular value of $S_0$")
    plt.title(rf"Final Beyn rank diagnostic, $\varepsilon={epsilon:.2f}$")
    plt.tight_layout()
    plt.savefig(output_directory / f"beyn_rank_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_beyn_convergence(
    epsilon: float,
    rows: list[BeynConvergenceRow],
    output_directory: Path,
    config: Config,
    aitken_estimate: float = math.nan,
) -> None:
    if not rows:
        return
    nq = np.array([row.quadrature_points for row in rows], dtype=int)

    plt.figure(figsize=(7, 4.5))
    for idx, attr in enumerate(("s0_sv1", "s0_sv2", "s0_sv3"), start=1):
        vals = np.array([getattr(row, attr) for row in rows], dtype=float)
        if np.any(np.isfinite(vals)):
            plt.semilogy(nq, vals, "o-", label=rf"$s_{idx}(S_0)$")
    plt.xlabel(r"Beyn quadrature points $N_q$")
    plt.ylabel(r"singular values of $S_0$")
    plt.title(rf"Beyn moment convergence, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        output_directory / f"beyn_s0_convergence_epsilon_{epsilon:.2f}.png", dpi=180
    )
    plt.close()

    raw_values = np.array([row.primary_raw_real for row in rows], dtype=float)
    finite = np.isfinite(raw_values)
    if np.any(finite):
        plt.figure(figsize=(7, 4.5))
        plt.plot(nq[finite], raw_values[finite], "o-", label="raw Beyn estimate")
        if np.isfinite(aitken_estimate):
            plt.axhline(aitken_estimate, linestyle=":", label=r"Aitken $\Delta^2$")
        _, right, _, _, _ = contour_geometry(config)
        plt.axhline(right, linestyle="-.", label="contour right endpoint")
        plt.axhline(kb_cutoff(config), linestyle="--", label=r"$\sqrt{\Lambda_1}b$")
        plt.xlabel(r"Beyn quadrature points $N_q$")
        plt.ylabel(r"raw Beyn $kb$")
        plt.title(rf"Beyn raw-eigenvalue convergence, $\varepsilon={epsilon:.2f}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            output_directory / f"beyn_candidate_convergence_epsilon_{epsilon:.2f}.png",
            dpi=180,
        )
        plt.close()


def plot_candidate_refinement(
    epsilon: float,
    rows: list[ModeRefinementRow],
    output_directory: Path,
) -> None:
    if not rows:
        return

    candidate_ids = sorted(set(row.candidate_index for row in rows))
    for candidate_index in candidate_ids:
        subset = [row for row in rows if row.candidate_index == candidate_index]
        M = np.array([row.M for row in subset])
        sv = np.array([row.sigma_min for row in subset])
        relative_sv = np.array([row.relative_singular_value for row in subset])
        kb = np.array([row.kb for row in subset])

        plt.figure(figsize=(7, 4.5))
        plt.semilogy(M, sv, "o-")
        plt.xlabel(r"$M$")
        plt.ylabel(r"$\sigma_{\min}(A(k_*))$")
        plt.title(
            rf"Beyn candidate {candidate_index}, $\varepsilon={epsilon:.2f}$: absolute SVD"
        )
        plt.tight_layout()
        plt.savefig(
            output_directory
            / f"candidate_{candidate_index}_svd_epsilon_{epsilon:.2f}.png",
            dpi=180,
        )
        plt.close()

        plt.figure(figsize=(7, 4.5))
        plt.semilogy(M, relative_sv, "o-")
        plt.xlabel(r"$M$")
        plt.ylabel(r"$\sigma_{\min}(A)/\sigma_{\max}(A)$")
        plt.title(
            rf"Beyn candidate {candidate_index}, $\varepsilon={epsilon:.2f}$: relative singularity"
        )
        plt.tight_layout()
        plt.savefig(
            output_directory
            / f"candidate_{candidate_index}_relative_svd_epsilon_{epsilon:.2f}.png",
            dpi=180,
        )
        plt.close()

        plt.figure(figsize=(7, 4.5))
        plt.plot(M, kb, "o-")
        plt.xlabel(r"$M$")
        plt.ylabel(r"$kb$")
        plt.title(
            rf"Beyn candidate {candidate_index}, $\varepsilon={epsilon:.2f}$: $kb$ convergence"
        )
        plt.tight_layout()
        plt.savefig(
            output_directory
            / f"candidate_{candidate_index}_kb_epsilon_{epsilon:.2f}.png",
            dpi=180,
        )
        plt.close()


def plot_summary(
    results: list[ValidationResult], output_directory: Path, config: Config
) -> None:
    eps = np.array([row.epsilon for row in results])
    counts = np.array([row.resolved_mode_count for row in results])
    finite_kb = [row for row in results if np.isfinite(row.kb_numerical)]

    if finite_kb:
        eps_kb = np.array([row.epsilon for row in finite_kb])
        kb_num = np.array([row.kb_numerical for row in finite_kb])
        kb_asym_finite = np.array([row.kb_asymptotic for row in finite_kb])
        plt.figure(figsize=(7, 4.5))
        plt.plot(eps_kb, kb_num, "o-", label=r"$kb_{\mathrm{BEM}}$")
        plt.plot(eps_kb, kb_asym_finite, "s--", label=r"$kb_{\mathrm{asym}}$")
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"$kb$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_directory / "kb_vs_epsilon.png", dpi=180)
        plt.close()

        sigma_over_eps2 = np.array(
            [row.sigma_over_epsilon_squared for row in finite_kb], dtype=float
        )
        coefficient = np.array(
            [row.asymptotic_coefficient for row in finite_kb], dtype=float
        )
        plt.figure(figsize=(7, 4.5))
        plt.plot(eps_kb, sigma_over_eps2, "o-", label=r"$\sigma_{\mathrm{BEM}}/\varepsilon^2$")
        plt.plot(eps_kb, coefficient, "--", label=r"$C(a)$")
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"$\sigma/\varepsilon^2$")
        plt.title(r"Direct test of $\sigma=C(a)\varepsilon^2+O(\varepsilon^3\log\varepsilon)$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_directory / "sigma_over_epsilon_squared.png", dpi=180)
        plt.close()

        scaled = np.array(
            [row.scaled_asymptotic_remainder for row in finite_kb], dtype=float
        )
        plt.figure(figsize=(7, 4.5))
        plt.plot(eps_kb, scaled, "o-")
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(
            r"$|\sigma_{\mathrm{BEM}}-C\varepsilon^2|/(\varepsilon^3|\log\varepsilon|)$"
        )
        plt.title("Scaled asymptotic remainder")
        plt.tight_layout()
        plt.savefig(output_directory / "scaled_asymptotic_remainder.png", dpi=180)
        plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, counts, "o-")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel("locally resolved discrete modes")
    plt.title("Beyn discovery + BEM local certification")
    plt.tight_layout()
    plt.savefig(output_directory / "resolved_mode_count.png", dpi=180)
    plt.close()

    finite = [row for row in results if np.isfinite(row.relative_error_sigma)]
    if finite:
        eps_f = np.array([row.epsilon for row in finite])
        error = np.array([row.relative_error_sigma for row in finite])
        plt.figure(figsize=(7, 4.5))
        plt.plot(eps_f, error, "o-")
        plt.axhline(
            config.relative_sigma_error_tolerance, linestyle="--", label="reporting tolerance"
        )
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"relative error in $\sigma$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_directory / "relative_error_sigma.png", dpi=180)
        plt.close()


def plot_internal_convergence(
    rows: list[InternalConvergenceRow], output_directory: Path
) -> None:
    if not rows:
        return
    parameters = sorted(set(row.parameter for row in rows))
    for parameter in parameters:
        subset = [row for row in rows if row.parameter == parameter]
        x = np.array([row.value for row in subset], dtype=float)
        y = np.array([row.relative_sigma_shift_from_baseline for row in subset], dtype=float)
        order = np.argsort(x)
        plt.figure(figsize=(7, 4.5))
        plt.semilogy(x[order], np.maximum(y[order], 1.0e-16), "o-")
        plt.xlabel(parameter)
        plt.ylabel("relative shift in sigma from baseline")
        plt.title(f"Internal numerical sensitivity: {parameter}")
        plt.tight_layout()
        plt.savefig(output_directory / f"internal_convergence_{parameter}.png", dpi=180)
        plt.close()


def plot_critical_height(
    rows: list[CriticalHeightRow], output_directory: Path
) -> None:
    finite = [row for row in rows if np.isfinite(row.a_critical_estimate)]
    if not finite:
        return
    eps = np.array([row.epsilon for row in finite])
    ac = np.array([row.a_critical_estimate for row in finite])
    a0 = finite[0].a0_star
    plt.figure(figsize=(7, 4.5))
    plt.plot(eps, ac, "o-", label=r"$a_c(\varepsilon)$ numerical")
    plt.axhline(a0, linestyle="--", label=r"$a_0^*$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"critical height $a$")
    plt.title(r"Critical-height test: $a_c(\varepsilon)\to a_0^*$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "critical_height_vs_epsilon.png", dpi=180)
    plt.close()


def print_result(result: ValidationResult, config: Config) -> None:
    print("\n  --- Beyn-v4 + local-SVD validation result ---")
    print(f"  effective cutoff margin        = {result.effective_beyn_cutoff_margin:.3e}")
    print(f"  final Beyn quadrature Nq       = {result.beyn_final_quadrature_points}")
    print(f"  final Beyn estimated rank      = {result.beyn_estimated_rank}")
    print(
        f"  Beyn rank stable (last 2 Nq)   = {'YES' if result.beyn_rank_stable else 'no'}"
    )
    if result.tighter_margin_rank_consistent is None:
        tighter_text = "N/A"
    else:
        tighter_text = "PASS" if result.tighter_margin_rank_consistent else "FAIL"
    print(f"  tighter-cutoff rank check      = {tighter_text}")
    print(f"  final strict Beyn candidates   = {result.beyn_near_real_candidates}")
    print(f"  local refinement seeds         = {result.local_refinement_seed_count}")
    print(f"  local seed source              = {result.local_refinement_seed_source}")
    if np.isfinite(result.aitken_estimate):
        print(f"  Aitken diagnostic kb           = {result.aitken_estimate:.12f}")
    else:
        print("  Aitken diagnostic kb           = --")
    print(f"  locally resolved modes         = {result.resolved_mode_count}")
    print(
        "  one-mode count supported      = "
        f"{'PASS' if result.unique_mode_verified else 'FAIL'}"
    )
    print(f"  kb asymptotic                  = {result.kb_asymptotic:.12f}")
    if np.isfinite(result.kb_numerical):
        print(f"  kb BEM                         = {result.kb_numerical:.12f}")
        print(f"  sigma asym                     = {result.sigma_asymptotic:.8e}")
        print(f"  sigma BEM                      = {result.sigma_numerical:.8e}")
        print(f"  sigma/epsilon^2                = {result.sigma_over_epsilon_squared:.8e}")
        print(f"  asymptotic coefficient C(a)    = {result.asymptotic_coefficient:.8e}")
        print(f"  scaled asymptotic remainder    = {result.scaled_asymptotic_remainder:.8e}")
        print(f"  final sigma_min(A)             = {result.sigma_min_final:.3e}")
        print(f"  final sigma_max(A)             = {result.sigma_max_final:.3e}")
        print(
            f"  final sigma_min/sigma_max      = {result.relative_singular_value_final:.3e}"
        )
        print(f"  final minimum drop             = {result.final_drop_factor:.2e}")
        print(
            f"  final mesh change in sigma     = {result.final_relative_mesh_change:.3%}"
        )
        print(f"  relative sigma error           = {result.relative_error_sigma:.3%}")
        if np.isfinite(result.wall_relative_residual):
            print(f"  wall Dirichlet residual        = {result.wall_relative_residual:.3e}")
            print(f"  obstacle Neumann residual      = {result.neumann_relative_residual:.3e}")
            print(
                f"  decay rate left/right          = "
                f"{result.decay_rate_left:.6e}/{result.decay_rate_right:.6e}"
            )
            decay_text = (
                "PASS" if result.physical_decay_verified is True else "FAIL"
            )
            print(f"  physical decay check           = {decay_text}")
    else:
        print("  kb BEM                         = --")
        print(f"  sigma asym                     = {result.sigma_asymptotic:.8e}")
        print("  sigma BEM                      = --")
        print("  relative singular value        = --")
        print("  relative sigma error           = --")

    if result.asymptotic_agreement_verified is None:
        asym_text = "N/A (requires exactly one certified BEM mode)"
    else:
        asym_text = "PASS" if result.asymptotic_agreement_verified else "FAIL"
    print(
        f"  asymptotic accuracy <= {config.relative_sigma_error_tolerance:.1%}: "
        f"{asym_text}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    a0_star = critical_height_leading_order(config)

    print("=== Theorem 2.1 numerical validation: Beyn + BEM + SVD (v4) ===")
    print(f"b = {config.b}")
    print(f"a = {config.a}")
    print(f"leading-order a0* = {a0_star:.12f}")
    print(f"geometric leading condition a > a0*: {'PASS' if config.a > a0_star else 'FAIL'}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"sqrt(Lambda_1)b = {kb_cutoff(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")
    print(
        f"Beyn global discovery uses fixed M={config.beyn_M} "
        "(no multi-M global Beyn sweep, by design)."
    )
    print(
        f"Beyn Nq base={config.beyn_quadrature_levels}, "
        f"adaptive={config.beyn_adaptive_quadrature_levels}, "
        f"probe_dim={config.beyn_probe_dim}"
    )
    print(
        f"rank diagnostics: empty-S0={config.beyn_empty_s0_tolerance:.1e}, "
        f"rel_tol={config.beyn_rank_relative_tolerance:.1e}, "
        f"gap_threshold={config.beyn_rank_gap_threshold:.1e}"
    )
    print(f"local BEM refinement M = {config.refinement_M}")

    if config.a <= a0_star:
        raise RuntimeError(
            "The main epsilon-scaling experiment is configured on the existence "
            "side and requires a > a0*. Use run_critical_height_study for the "
            "zero/one-mode transition around a*."
        )

    first_config = config_for_epsilon(config.epsilon_values[0], config)
    print("\nChecking complex-kb support required by Beyn...")
    complex_support_smoke_test(config.epsilon_values[0], first_config)
    print("  complex-kb assembly: PASS")

    summaries: list[ValidationResult] = []
    all_beyn_diagnostics: list[BeynDiagnostics] = []
    all_eigen_rows: list[BeynEigenvalueRow] = []
    all_convergence_rows: list[BeynConvergenceRow] = []
    all_refinement_rows: list[ModeRefinementRow] = []
    all_mode_results: list[ModeResult] = []
    all_physical: list[PhysicalDiagnostics] = []

    for epsilon in config.epsilon_values:
        work_config = config_for_epsilon(epsilon, config)
        left, right, center, rx, ry = contour_geometry(work_config)
        print(f"\n=== epsilon={epsilon:.8f} ===")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, work_config)
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted sigma = {sigma_asym:.8e}")
        print(f"  predicted cutoff gap = {kb_cutoff(work_config) - kb_asym:.3e}")
        print(
            "  Beyn ellipse: "
            f"real=[{left:.12f}, {right:.12f}], "
            f"center={center:.12f}, rx={rx:.12f}, ry={ry:.3e}"
        )

        (
            result,
            beyn_diags,
            eigen_rows,
            convergence_rows,
            refinement_rows,
            mode_results,
            final_s0_singular_values,
            final_eigen_rows,
            aitken_estimate,
            physical,
        ) = validate_epsilon(epsilon, config)

        summaries.append(result)
        all_beyn_diagnostics.extend(beyn_diags)
        all_eigen_rows.extend(eigen_rows)
        all_convergence_rows.extend(convergence_rows)
        all_refinement_rows.extend(refinement_rows)
        all_mode_results.extend(mode_results)
        if physical is not None:
            all_physical.append(physical)

        print_result(result, work_config)
        plot_beyn_contour(epsilon, final_eigen_rows, output_directory, work_config)
        plot_s0_singular_values(epsilon, final_s0_singular_values, output_directory)
        plot_beyn_convergence(
            epsilon, convergence_rows, output_directory, work_config, aitken_estimate
        )
        plot_candidate_refinement(epsilon, refinement_rows, output_directory)

    write_dataclass_csv(output_directory / "summary.csv", summaries)
    write_dataclass_csv(output_directory / "beyn_diagnostics.csv", all_beyn_diagnostics)
    write_dataclass_csv(output_directory / "beyn_convergence.csv", all_convergence_rows)
    write_dataclass_csv(output_directory / "beyn_raw_eigenvalues.csv", all_eigen_rows)
    write_dataclass_csv(output_directory / "mode_refinement.csv", all_refinement_rows)
    write_dataclass_csv(output_directory / "mode_results.csv", all_mode_results)
    write_dataclass_csv(output_directory / "physical_diagnostics.csv", all_physical)
    plot_summary(summaries, output_directory, config)

    internal_rows = run_internal_convergence_study(summaries, config)
    write_dataclass_csv(output_directory / "internal_convergence.csv", internal_rows)
    plot_internal_convergence(internal_rows, output_directory)

    critical_rows = run_critical_height_study(config)
    write_dataclass_csv(output_directory / "critical_height.csv", critical_rows)
    plot_critical_height(critical_rows, output_directory)

    one_mode_supported = sum(row.unique_mode_verified for row in summaries)
    asymptotic_evaluated = [
        row for row in summaries if row.asymptotic_agreement_verified is not None
    ]
    asymptotic_passed = sum(
        row.asymptotic_agreement_verified is True for row in asymptotic_evaluated
    )
    physical_evaluated = [row for row in summaries if row.physical_decay_verified is not None]
    physical_passed = sum(row.physical_decay_verified is True for row in physical_evaluated)

    print("\n=== FINAL SUMMARY ===")
    print("Global discovery: adaptive Beyn contour method with Nq convergence study")
    print("Global Beyn BEM order: one fixed M, as requested")
    print("Local certification: relative SVD singularity + BEM mesh refinement")
    print("Near-cutoff safeguard: epsilon-adaptive margin + tighter-margin rank check")
    print("Asymptotics: direct sigma/epsilon^2 and scaled O(epsilon^3 log epsilon) diagnostics")
    print(
        f"Stable rank=1 + exactly one locally resolved mode: "
        f"{one_mode_supported}/{len(summaries)}"
    )
    print(
        f"Leading asymptotic sigma within {config.relative_sigma_error_tolerance:.1%}: "
        f"{asymptotic_passed}/{len(asymptotic_evaluated)}"
    )
    if physical_evaluated:
        print(
            f"Independent physical decay checks: "
            f"{physical_passed}/{len(physical_evaluated)}"
        )
    if not config.run_critical_height_study:
        print(
            "Critical-height a*(epsilon) study: IMPLEMENTED but disabled. "
            "Set run_critical_height_study=True for the full existence/non-existence transition test."
        )

    for row in summaries:
        sigma_error_text = (
            f"{row.relative_error_sigma:.3%}"
            if np.isfinite(row.relative_error_sigma)
            else "N/A"
        )
        asym_text = (
            "N/A"
            if row.asymptotic_agreement_verified is None
            else ("PASS" if row.asymptotic_agreement_verified else "FAIL")
        )
        margin_text = (
            "N/A"
            if row.tighter_margin_rank_consistent is None
            else ("PASS" if row.tighter_margin_rank_consistent else "FAIL")
        )
        print(
            f"  epsilon={row.epsilon:.5f}: "
            f"margin={row.effective_beyn_cutoff_margin:.2e}, "
            f"Nq={row.beyn_final_quadrature_points}, "
            f"rank={row.beyn_estimated_rank}, "
            f"rank_stable={'yes' if row.beyn_rank_stable else 'no'}, "
            f"margin_check={margin_text}, "
            f"resolved={row.resolved_mode_count}, "
            f"one-mode={'PASS' if row.unique_mode_verified else 'FAIL'}, "
            f"asymptotic={asym_text}, sigma_error={sigma_error_text}"
        )

    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
