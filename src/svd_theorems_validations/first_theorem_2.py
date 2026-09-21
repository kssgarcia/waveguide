"""
Focused numerical validation of Theorem 2.1: one discrete trapped mode.

The script is intentionally limited to the circular obstacle and the interval

    0 < k^2 < Lambda_1.

It separates three different numerical questions:

1. EXPECTED MODE:
   Does the local minimum predicted by the asymptotic formula persist and
   converge when the BEM boundary discretization M is refined?

2. UNIQUENESS SCREEN:
   Does a whole-band scan reveal any additional resolved near-singular
   minima away from the expected mode?

3. ASYMPTOTIC ACCURACY:
   How close is the converged numerical spectral distance sigma_BEM to the
   leading asymptotic prediction sigma_asym?

Important methodological point
------------------------------
Existence is NOT decided by a single hard test such as sigma_min < 1e-6.
A numerical mode is considered resolved only when:

- sigma_min(A(k)) has a genuine local dip,
- the matrix is sufficiently close to singular at the refined minimum,
- the location of the mode stabilizes under mesh refinement.

The 5% asymptotic-error criterion is reported separately. Therefore a mode
may be numerically resolved even when the leading asymptotic approximation
has already lost 5% accuracy.

This is a numerical validation, not a mathematical proof of uniqueness.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

# Same project import convention used by the previous scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice

PI = np.pi


@dataclass(frozen=True)
class Config:
    # Geometry used in the paper experiment.
    b: float = 1.0
    a: float = 0.6

    # Keep the first experiment small and interpretable.
    epsilon_values: tuple[float, ...] = (
        0.01, 0.03111111, 0.05222222, 0.07333333, 0.09444444,
        0.11555556, 0.13666667, 0.15777778, 0.17888889, 0.2,
    )

    # BEM / Green function.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Mesh refinement used to validate the expected branch.
    refinement_M: tuple[int, ...] = (16, 24, 32, 40, 48)

    # Expected-mode search window in spectral distance
    # delta = Lambda_1 - k^2 = sigma^2.
    expected_delta_lower_factor: float = 0.20
    expected_delta_upper_factor: float = 5.00
    cutoff_delta_floor: float = 1.0e-9
    minimizer_xatol: float = 2.0e-10

    # A resolved minimum should be substantially lower than the values at
    # the two ends of its search bracket. 100 means two orders of magnitude.
    minimum_drop_factor: float = 100.0

    # This is deliberately looser than the old 1e-6 binary threshold.
    # It is only one of several diagnostics, not the sole existence test.
    near_singular_tolerance: float = 1.0e-4

    # Mesh-convergence criterion.  We compare sigma_BEM rather than kb because
    # kb is extremely close to the cutoff and can hide meaningful differences.
    mesh_sigma_relative_tolerance: float = 0.03

    # Accuracy of the leading asymptotic approximation; reported separately.
    relative_sigma_error_tolerance: float = 0.05

    # Whole-band uniqueness screen. This is cheap discovery with a low-order
    # BEM, followed by limited refinement only for sampled local minima.
    uniqueness_scan_M: int = 16
    uniqueness_refine_M: int = 32
    uniqueness_linear_points: int = 32
    uniqueness_log_points: int = 32
    low_k_margin: float = 1.0e-5

    # Two-resolution persistence test for additional candidates.
    additional_candidate_kb_tolerance: float = 2.0e-3

    output_directory: str = "theorem_2_1_refined_validation"


@dataclass
class RefinementRow:
    epsilon: float
    M: int
    kb: float
    sigma_bem: float
    sigma_min: float
    left_value: float
    right_value: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_previous: float


@dataclass
class AdditionalCandidate:
    epsilon: float
    bracket_left: float
    bracket_right: float
    kb_coarse: float
    sigma_min_coarse: float
    kb_fine: float
    sigma_min_fine: float
    fine_drop_factor: float
    kb_shift: float
    persistent: bool


@dataclass
class ValidationResult:
    epsilon: float
    a: float
    lambda_1: float
    kb_cutoff: float
    kb_asymptotic: float
    kb_numerical: float
    sigma_asymptotic: float
    sigma_numerical: float
    sigma_min_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    relative_error_kb: float
    relative_error_sigma: float
    expected_mode_resolved: bool
    additional_resolved_modes: int
    uniqueness_screen_passed: bool
    unique_mode_verified: bool
    asymptotic_agreement_verified: bool


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


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    """Return (kb_asymptotic, sigma_asymptotic)."""
    alpha = PI * config.a / config.b
    mu = 1.0
    area = PI

    bracket = (
        PI * mu * math.sin(alpha / 2.0) ** 2 - 0.5 * area * math.cos(alpha / 2.0) ** 2
    )
    sigma = epsilon**2 * PI**2 / (4.0 * config.b**3) * bracket

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


def kb_from_delta(delta: float, config: Config) -> float:
    """delta = Lambda_1-k^2; return kb inside the discrete band."""
    lam1 = lambda_1(config)
    delta = min(max(delta, config.cutoff_delta_floor), lam1 * (1.0 - 1.0e-12))
    return config.b * math.sqrt(max(lam1 - delta, 0.0))


# ---------------------------------------------------------------------------
# BEM assembly
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
    kb: float,
    config: Config,
) -> tuple[Callable[..., complex], Callable[..., complex]]:
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


def assemble_matrix(kb: float, epsilon: float, M: int, config: Config) -> np.ndarray:
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

    # Full contour [0,2pi]:
    #   1/2 u_i = (2pi/M) sum_j K^w_ij u_j
    # therefore A = I - (4pi/M) K^w.
    return np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted


def smallest_singular_value(kb: float, epsilon: float, M: int, config: Config) -> float:
    A = assemble_matrix(kb, epsilon, M, config)
    singular_values = np.linalg.svd(A, compute_uv=False)
    return float(singular_values[-1])


# ---------------------------------------------------------------------------
# Expected branch: targeted mesh-refinement experiment
# ---------------------------------------------------------------------------


def expected_mode_bracket(epsilon: float, config: Config) -> tuple[float, float]:
    """Bracket the predicted mode in delta = Lambda_1-k^2 coordinates."""
    _, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2
    lam1 = lambda_1(config)

    # Larger delta means smaller k, hence left/right are reversed in delta.
    delta_far = min(delta_asym * config.expected_delta_upper_factor, lam1 * 0.95)
    delta_near = max(
        delta_asym * config.expected_delta_lower_factor, config.cutoff_delta_floor
    )

    left = kb_from_delta(delta_far, config)
    right = kb_from_delta(delta_near, config)

    if not left < right:
        raise RuntimeError("Invalid expected-mode search bracket.")
    return left, right


def refine_expected_mode_for_M(
    epsilon: float,
    M: int,
    config: Config,
) -> tuple[float, float, float, float, float, bool]:
    """
    Return:
      kb, sigma_min, left_value, right_value, drop_factor, minimum_is_interior
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
            f"Expected-mode minimization failed for M={M}: {result.message}"
        )

    kb = float(result.x)
    sigma_min = smallest_singular_value(kb, epsilon, M, config)

    endpoint_level = min(left_value, right_value)
    drop_factor = endpoint_level / max(sigma_min, tiny)

    # Reject a minimum that simply lands on a bracket endpoint.
    width = right - left
    edge_margin = 0.01 * width
    interior = (kb > left + edge_margin) and (kb < right - edge_margin)

    return kb, sigma_min, left_value, right_value, drop_factor, interior


def run_expected_branch_refinement(
    epsilon: float,
    config: Config,
) -> list[RefinementRow]:
    rows: list[RefinementRow] = []
    previous_sigma = math.nan

    print("  expected-mode mesh refinement:")
    for M in config.refinement_M:
        kb, sigma_min, left_value, right_value, drop_factor, interior = (
            refine_expected_mode_for_M(epsilon, M, config)
        )
        sigma_bem = sigma_from_kb(kb, config)

        change = (
            abs(sigma_bem - previous_sigma) / max(abs(sigma_bem), 1.0e-30)
            if np.isfinite(previous_sigma)
            else math.nan
        )

        row = RefinementRow(
            epsilon=epsilon,
            M=M,
            kb=kb,
            sigma_bem=sigma_bem,
            sigma_min=sigma_min,
            left_value=left_value,
            right_value=right_value,
            drop_factor=drop_factor,
            minimum_is_interior=interior,
            relative_sigma_change_from_previous=change,
        )
        rows.append(row)
        previous_sigma = sigma_bem

        change_text = f"{change:.3%}" if np.isfinite(change) else "--"
        print(
            f"    M={M:>2}: kb={kb:.12f}, "
            f"sigma_BEM={sigma_bem:.8e}, sigma_min={sigma_min:.3e}, "
            f"drop={drop_factor:.2e}, mesh_change={change_text}, "
            f"interior={'yes' if interior else 'no'}"
        )

    return rows


def expected_mode_is_resolved(rows: list[RefinementRow], config: Config) -> bool:
    if len(rows) < 2:
        return False

    final = rows[-1]
    previous = rows[-2]

    final_change = abs(final.sigma_bem - previous.sigma_bem) / max(
        abs(final.sigma_bem), 1.0e-30
    )

    near_singular = final.sigma_min <= config.near_singular_tolerance
    deep_minimum = final.drop_factor >= config.minimum_drop_factor
    converged = final_change <= config.mesh_sigma_relative_tolerance

    return bool(
        final.minimum_is_interior and near_singular and deep_minimum and converged
    )


# ---------------------------------------------------------------------------
# Whole-band uniqueness screen
# ---------------------------------------------------------------------------


def build_uniqueness_grid(epsilon: float, config: Config) -> np.ndarray:
    lam1 = lambda_1(config)
    kb1 = kb_cutoff(config)
    _, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2

    linear = np.linspace(
        config.low_k_margin,
        kb1 - config.low_k_margin,
        config.uniqueness_linear_points,
    )

    delta_max = lam1 * (1.0 - 1.0e-8)
    deltas = np.geomspace(
        config.cutoff_delta_floor,
        delta_max,
        config.uniqueness_log_points,
    )
    logarithmic = config.b * np.sqrt(np.maximum(lam1 - deltas, 0.0))

    # Explicit points around the predicted mode so that the expected dip is
    # represented in the full-band diagnostic plot.
    local_factors = np.array([0.15, 0.25, 0.5, 0.75, 1.0, 1.5, 2.5, 5.0])
    local_delta = delta_asym * local_factors
    local_delta = local_delta[
        (local_delta > config.cutoff_delta_floor) & (local_delta < lam1)
    ]
    local = config.b * np.sqrt(lam1 - local_delta)

    grid = np.concatenate([linear, logarithmic, local])
    grid = grid[(grid > 0.0) & (grid < kb1)]
    return np.unique(np.sort(grid))


def sampled_local_minimum_brackets(
    scan: np.ndarray,
    values: np.ndarray,
) -> list[tuple[float, float]]:
    """Only interior sampled minima; open band edges are not treated as modes."""
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
) -> tuple[float, float, float]:
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
        return math.nan, math.nan, math.nan

    kb = float(result.x)
    sigma_min = smallest_singular_value(kb, epsilon, M, config)
    drop_factor = min(left_value, right_value) / max(sigma_min, tiny)
    return kb, sigma_min, drop_factor


def whole_band_uniqueness_screen(
    epsilon: float,
    config: Config,
) -> tuple[np.ndarray, np.ndarray, list[AdditionalCandidate]]:
    scan = build_uniqueness_grid(epsilon, config)
    values = np.empty(len(scan), dtype=float)

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
        f"{len(other_brackets)} outside expected-mode window"
    )

    additional: list[AdditionalCandidate] = []
    for j, (left, right) in enumerate(other_brackets, start=1):
        print(f"    checking additional minimum {j}/{len(other_brackets)}")

        kb_c, sv_c, _ = refine_bracket_once(
            left, right, epsilon, config.uniqueness_scan_M, config
        )
        kb_f, sv_f, drop_f = refine_bracket_once(
            left, right, epsilon, config.uniqueness_refine_M, config
        )

        if np.isfinite(kb_c) and np.isfinite(kb_f):
            kb_shift = abs(kb_f - kb_c)
        else:
            kb_shift = math.nan

        persistent = bool(
            np.isfinite(kb_f)
            and np.isfinite(sv_f)
            and sv_f <= config.near_singular_tolerance
            and drop_f >= config.minimum_drop_factor
            and np.isfinite(kb_shift)
            and kb_shift <= config.additional_candidate_kb_tolerance
        )

        additional.append(
            AdditionalCandidate(
                epsilon=epsilon,
                bracket_left=left,
                bracket_right=right,
                kb_coarse=kb_c,
                sigma_min_coarse=sv_c,
                kb_fine=kb_f,
                sigma_min_fine=sv_f,
                fine_drop_factor=drop_f,
                kb_shift=kb_shift,
                persistent=persistent,
            )
        )

        print(
            f"      kb(M={config.uniqueness_scan_M})={kb_c:.12f}, sigma_min={sv_c:.3e}"
        )
        print(
            f"      kb(M={config.uniqueness_refine_M})={kb_f:.12f}, "
            f"sigma_min={sv_f:.3e}, drop={drop_f:.2e}, "
            f"shift={kb_shift:.3e}, persistent={'YES' if persistent else 'no'}"
        )

    return scan, values, additional


# ---------------------------------------------------------------------------
# Validation and reporting
# ---------------------------------------------------------------------------


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
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)

    refinement = run_expected_branch_refinement(epsilon, config)
    expected_resolved = expected_mode_is_resolved(refinement, config)

    final = refinement[-1]
    kb_num = final.kb
    sigma_num = final.sigma_bem
    final_mesh_change = final.relative_sigma_change_from_previous

    error_kb = abs(kb_num - kb_asym) / abs(kb_asym)
    error_sigma = abs(sigma_num - sigma_asym) / abs(sigma_asym)
    asymptotic_ok = error_sigma <= config.relative_sigma_error_tolerance

    scan, scan_values, additional = whole_band_uniqueness_screen(epsilon, config)
    additional_resolved = [
        candidate for candidate in additional if candidate.persistent
    ]
    uniqueness_screen_passed = len(additional_resolved) == 0
    unique_mode_verified = expected_resolved and uniqueness_screen_passed

    result = ValidationResult(
        epsilon=epsilon,
        a=config.a,
        lambda_1=lambda_1(config),
        kb_cutoff=kb_cutoff(config),
        kb_asymptotic=kb_asym,
        kb_numerical=kb_num,
        sigma_asymptotic=sigma_asym,
        sigma_numerical=sigma_num,
        sigma_min_final=final.sigma_min,
        final_drop_factor=final.drop_factor,
        final_relative_mesh_change=final_mesh_change,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        expected_mode_resolved=expected_resolved,
        additional_resolved_modes=len(additional_resolved),
        uniqueness_screen_passed=uniqueness_screen_passed,
        unique_mode_verified=unique_mode_verified,
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
    plt.semilogy(scan, values, "o-", markersize=3, label=r"coarse $\sigma_{\min}(A)$")
    plt.axvline(result.kb_asymptotic, linestyle="--", label=r"$kb_{\mathrm{asym}}$")
    plt.axvline(result.kb_numerical, linestyle=":", label=r"$kb_{\mathrm{BEM}}$")
    plt.xlabel(r"$kb$")
    plt.ylabel(r"$\sigma_{\min}(A(k))$")
    plt.title(rf"Theorem 2.1 full-band screen, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"full_band_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_refinement(
    epsilon: float,
    rows: list[RefinementRow],
    sigma_asym: float,
    output_directory: Path,
) -> None:
    M = np.array([row.M for row in rows])
    sigma_bem = np.array([row.sigma_bem for row in rows])
    sigma_min = np.array([row.sigma_min for row in rows])

    plt.figure(figsize=(7, 4.5))
    plt.plot(M, sigma_bem, "o-", label=r"$\sigma_{\mathrm{BEM}}$")
    plt.axhline(sigma_asym, linestyle="--", label=r"$\sigma_{\mathrm{asym}}$")
    plt.xlabel(r"$M$")
    plt.ylabel(r"spectral distance $\sigma$")
    plt.title(rf"Mesh refinement, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        output_directory / f"refinement_sigma_epsilon_{epsilon:.2f}.png", dpi=180
    )
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(M, sigma_min, "o-")
    plt.axhline(
        CONFIG.near_singular_tolerance, linestyle="--", label="diagnostic tolerance"
    )
    plt.xlabel(r"$M$")
    plt.ylabel(r"$\sigma_{\min}(A(k_*))$")
    plt.title(rf"Near-singularity under refinement, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"refinement_svd_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_summary(results: list[ValidationResult], output_directory: Path) -> None:
    eps = np.array([row.epsilon for row in results])
    sigma_asym = np.array([row.sigma_asymptotic for row in results])
    sigma_num = np.array([row.sigma_numerical for row in results])
    kb_asym = np.array([row.kb_asymptotic for row in results])
    kb_num = np.array([row.kb_numerical for row in results])
    error = np.array([row.relative_error_sigma for row in results])

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
        CONFIG.relative_sigma_error_tolerance, linestyle="--", label="5% tolerance"
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"relative error in $\sigma$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "relative_error_sigma.png", dpi=180)
    plt.close()


def print_result(result: ValidationResult, config: Config) -> None:
    print("\n  --- validation result ---")
    print(
        f"  expected mode resolved: "
        f"{'PASS' if result.expected_mode_resolved else 'FAIL'}"
    )
    print(f"  additional resolved modes: {result.additional_resolved_modes}")
    print(
        f"  uniqueness screen: {'PASS' if result.uniqueness_screen_passed else 'FAIL'}"
    )
    print(
        f"  exactly one resolved discrete mode: "
        f"{'PASS' if result.unique_mode_verified else 'FAIL'}"
    )

    print(f"  kb asymptotic = {result.kb_asymptotic:.12f}")
    print(f"  kb BEM        = {result.kb_numerical:.12f}")
    print(f"  sigma asym    = {result.sigma_asymptotic:.8e}")
    print(f"  sigma BEM     = {result.sigma_numerical:.8e}")
    print(f"  final sigma_min(A) = {result.sigma_min_final:.3e}")
    print(f"  final minimum drop = {result.final_drop_factor:.2e}")
    print(f"  final mesh change in sigma = {result.final_relative_mesh_change:.3%}")
    print(f"  relative sigma error = {result.relative_error_sigma:.3%}")
    print(
        f"  asymptotic accuracy <= {config.relative_sigma_error_tolerance:.1%}: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )

    # Deliberately separate theorem existence/uniqueness from approximation accuracy.
    print(
        f"  EXISTENCE/UNIQUENESS CHECK: "
        f"{'PASS' if result.unique_mode_verified else 'FAIL'}"
    )
    print(
        f"  ASYMPTOTIC APPROXIMATION CHECK: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    a0_star = critical_height_leading_order(config)

    print("=== Focused validation of Theorem 2.1 ===")
    print(f"b = {config.b}")
    print(f"a = {config.a}")
    print(f"leading-order a0* = {a0_star:.12f}")
    print(f"geometric condition a > a0*: {'PASS' if config.a > a0_star else 'FAIL'}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"sqrt(Lambda_1) b = {kb_cutoff(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")
    print(f"refinement M = {config.refinement_M}")

    if config.a <= a0_star:
        raise RuntimeError(
            "The selected geometry does not satisfy the leading-order condition a > a0*."
        )

    summary: list[ValidationResult] = []
    all_refinement: list[RefinementRow] = []
    all_additional: list[AdditionalCandidate] = []

    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.3f} ===")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted sigma = {sigma_asym:.8e}")
        print(f"  predicted kb cutoff gap = {kb_cutoff(config) - kb_asym:.3e}")

        result, refinement, scan, scan_values, additional = validate_epsilon(
            epsilon, config
        )

        summary.append(result)
        all_refinement.extend(refinement)
        all_additional.extend(additional)

        print_result(result, config)
        plot_full_band_scan(epsilon, scan, scan_values, result, output_directory)
        plot_refinement(epsilon, refinement, sigma_asym, output_directory)

    write_dataclass_csv(output_directory / "summary.csv", summary)
    write_dataclass_csv(output_directory / "mesh_refinement.csv", all_refinement)
    write_dataclass_csv(output_directory / "additional_candidates.csv", all_additional)
    plot_summary(summary, output_directory)

    uniqueness_passed = sum(row.unique_mode_verified for row in summary)
    asymptotic_passed = sum(row.asymptotic_agreement_verified for row in summary)

    print("\n=== FINAL SUMMARY ===")
    print("Interval screened: 0 < k^2 < Lambda_1")
    print(
        f"Existence + no additional resolved modes: {uniqueness_passed}/{len(summary)}"
    )
    print(
        f"Leading asymptotic sigma within "
        f"{config.relative_sigma_error_tolerance:.1%}: "
        f"{asymptotic_passed}/{len(summary)}"
    )

    for row in summary:
        print(
            f"  epsilon={row.epsilon:.3f}: "
            f"unique={'PASS' if row.unique_mode_verified else 'FAIL'}, "
            f"asymptotic={'PASS' if row.asymptotic_agreement_verified else 'FAIL'}, "
            f"sigma_min={row.sigma_min_final:.3e}, "
            f"mesh_change={row.final_relative_mesh_change:.3%}, "
            f"sigma_error={row.relative_error_sigma:.3%}"
        )

    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
