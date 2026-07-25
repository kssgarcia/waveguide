"""
Reviewer-oriented validation for trapped modes and BICs.

This script extends the original determinant-only experiment so that the
conference paper can report:

1. Quantitative BEM-versus-asymptotic comparisons.
2. Smallest-singular-value diagnostics.
3. Mesh-refinement stability with respect to the number of boundary points.
4. Cancellation of the first propagating waveguide component for BICs.
5. Separate distances to Lambda_1 (discrete modes) and Lambda_2 (BICs).
6. CSV tables and LaTeX-ready tables.

It assumes that the project modules `dipole_theorem.py` and `lattice_sums.py`
are available through the same relative path used by the original script.

Important:
- A small singular value is a candidate criterion, not by itself a proof that
  the mode is physical.
- For a BIC, the candidate is complemented by the projection of the reconstructed
  field onto the first open transverse mode.
- The tolerances below are numerical reporting criteria. They should be stated
  explicitly in the paper and, if necessary, agreed with the authors.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Literal

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

import dipole_theorem as dip
import lattice_sums as lattice


PI = np.pi
CaseName = Literal["discrete_circle", "bic_x", "bic_y"]


@dataclass(frozen=True)
class NumericalConfig:
    b: float = 1.0
    shape_beta: float = 0.1

    # BEM / Green-function discretization
    boundary_points: int = 32
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Spectral search
    scan_points: int = 25
    search_half_width: float = 0.02
    minimizer_xatol: float = 1.0e-8

    # Reporting / diagnostic tolerances
    relative_error_tolerance: float = 0.05  # 5%; state this explicitly
    singular_value_tolerance: float = 1.0e-6
    propagating_ratio_tolerance: float = 1.0e-3

    # Far-section projection used for BIC diagnostics
    projection_x_over_b: float = 2.0
    projection_quadrature_points: int = 80

    # Representative mesh-refinement experiment
    refinement_M: tuple[int, ...] = (16, 24, 32, 48, 64)
    refinement_epsilon: float = 0.10

    # Output
    output_directory: str = "reviewer_validation_output"


@dataclass
class SpectralResult:
    case: str
    epsilon: float
    a: float
    M: int
    kb_asymptotic: float
    kb_numerical: float
    sigma_asymptotic: float
    sigma_numerical: float
    relative_error_k: float
    relative_error_sigma: float
    sigma_min: float
    svd_score: float
    propagating_ratio_left: float
    propagating_ratio_right: float
    first_mode_coefficient_left: float
    first_mode_coefficient_right: float
    field_norm_left: float
    field_norm_right: float
    passes_singular_value: bool
    passes_propagating_test: bool
    accepted_local_candidate: bool


@dataclass
class RefinementResult:
    case: str
    epsilon: float
    M: int
    kb_numerical: float
    sigma_min: float
    relative_change_from_previous: float
    propagating_ratio_left: float
    propagating_ratio_right: float


def obstacle_geometry(
    t: np.ndarray | float,
    epsilon: float,
    shape_beta: float,
    case: CaseName,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return X, Y and their first and second derivatives."""
    t = np.asarray(t)

    if case == "discrete_circle":
        X = epsilon * np.cos(t)
        Y = epsilon * np.sin(t)
        Xp = -epsilon * np.sin(t)
        Yp = epsilon * np.cos(t)
        Xpp = -epsilon * np.cos(t)
        Ypp = -epsilon * np.sin(t)
        return X, Y, Xp, Yp, Xpp, Ypp

    if case == "bic_x":
        X = epsilon * (np.cos(t) - 0.5 * shape_beta * np.cos(2.0 * t))
        Y = epsilon * (np.sin(t) - 0.5 * shape_beta * np.sin(2.0 * t))
        Xp = epsilon * (-np.sin(t) + shape_beta * np.sin(2.0 * t))
        Yp = epsilon * (np.cos(t) - shape_beta * np.cos(2.0 * t))
        Xpp = epsilon * (-np.cos(t) + 2.0 * shape_beta * np.cos(2.0 * t))
        Ypp = epsilon * (-np.sin(t) + 2.0 * shape_beta * np.sin(2.0 * t))
        return X, Y, Xp, Yp, Xpp, Ypp

    if case == "bic_y":
        X = epsilon * (np.sin(t) - 0.5 * shape_beta * np.sin(2.0 * t))
        Y = epsilon * (-np.cos(t) + 0.5 * shape_beta * np.cos(2.0 * t))
        Xp = epsilon * (np.cos(t) - shape_beta * np.cos(2.0 * t))
        Yp = epsilon * (np.sin(t) - shape_beta * np.sin(2.0 * t))
        Xpp = epsilon * (-np.sin(t) + 2.0 * shape_beta * np.sin(2.0 * t))
        Ypp = epsilon * (np.cos(t) - 2.0 * shape_beta * np.cos(2.0 * t))
        return X, Y, Xp, Yp, Xpp, Ypp

    raise ValueError(f"Unsupported case: {case}")


def obstacle_position(case: CaseName, epsilon: float, shape_beta: float) -> float:
    """Vertical obstacle position used in each asymptotic case."""
    if case == "discrete_circle":
        return 0.6
    if case == "bic_x":
        return 0.0
    if case == "bic_y":
        # Leading-order approximation a = epsilon a_1,
        # with a_1 = -beta/12 + O(beta^2).
        return epsilon * (-shape_beta / 12.0)
    raise ValueError(case)


def boundary_nodes(M: int) -> np.ndarray:
    """Midpoint collocation nodes on the full interval [0, 2*pi]."""
    return (np.arange(M) + 0.5) * (2.0 * PI / M)


@lru_cache(maxsize=None)
def dipole_moment(case: CaseName, shape_beta: float) -> tuple[float, float]:
    if case == "discrete_circle":
        return 1.0, 0.0
    symmetry = "x" if case == "bic_x" else "y"
    mu, nu = dip.dipole(shape_beta, 0.0, 0.0, symmetry)
    return float(mu), float(nu)


def asymptotic_prediction(
    case: CaseName,
    epsilon: float,
    config: NumericalConfig,
) -> tuple[float, float, float, float]:
    """Return (kb_asymptotic, sigma_asymptotic, threshold, a)."""
    b = config.b
    a = obstacle_position(case, epsilon, config.shape_beta)
    lambda_1 = (PI / (2.0 * b)) ** 2
    lambda_2 = (PI / b) ** 2

    if case == "discrete_circle":
        mu = 1.0
        area = PI
        alpha = PI * a / b
        bracket = (
            PI * mu * np.sin(alpha / 2.0) ** 2
            - 0.5 * area * np.cos(alpha / 2.0) ** 2
        )
        sigma = epsilon**2 * PI**2 / (4.0 * b**3) * bracket
        if sigma <= 0:
            raise ValueError(
                "The selected discrete-mode geometry does not produce a "
                "positive leading-order decay rate."
            )
        threshold = lambda_1
    else:
        mu, _ = dipole_moment(case, config.shape_beta)
        sigma = epsilon**2 * PI**3 / b**3 * mu
        threshold = lambda_2

    k_squared = threshold - sigma**2
    if k_squared <= 0:
        raise ValueError(
            f"Non-positive asymptotic k^2 for case={case}, epsilon={epsilon}."
        )

    return float(np.sqrt(k_squared) * b), float(sigma), float(threshold), float(a)


def make_green_functions(
    kb: float,
    a: float,
    config: NumericalConfig,
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
    """Centred finite differences with respect to source coordinates."""
    dG_dxi = (G(x, y, xi + h, eta) - G(x, y, xi - h, eta)) / (2.0 * h)
    dG_deta = (G(x, y, xi, eta + h) - G(x, y, xi, eta - h)) / (2.0 * h)
    return dG_dxi, dG_deta


def weighted_normal_kernel(
    psi: float,
    theta: float,
    epsilon: float,
    case: CaseName,
    config: NumericalConfig,
    G: Callable[..., complex],
    G_regularized: Callable[..., complex],
) -> complex:
    """Compute w(theta) * dG/dn_q, including diagonal regularization."""
    x, y, _, _, _, _ = obstacle_geometry(psi, epsilon, config.shape_beta, case)
    xi, eta, xi_p, eta_p, xi_pp, eta_pp = obstacle_geometry(
        theta, epsilon, config.shape_beta, case
    )
    w = float(np.hypot(xi_p, eta_p))
    h = config.finite_difference_step

    if abs(psi - theta) > 1.0e-10:
        G_xi, G_eta = source_derivatives(G, x, y, xi, eta, h)
        return xi_p * G_eta - eta_p * G_xi

    G_xi_reg, G_eta_reg = source_derivatives(G_regularized, x, y, xi, eta, h)
    geometric_term = (xi_pp * eta_p - eta_pp * xi_p) / (4.0 * PI * w**2)
    regularized_term = xi_p * G_eta_reg - eta_p * G_xi_reg
    return geometric_term + regularized_term


def assemble_matrix(
    kb: float,
    epsilon: float,
    case: CaseName,
    M: int,
    config: NumericalConfig,
) -> tuple[np.ndarray, np.ndarray, Callable[..., complex]]:
    # Matrix assembly only needs the geometric position.  Do not recompute
    # the auxiliary dipole BEM inside every spectral evaluation.
    a = obstacle_position(case, epsilon, config.shape_beta)
    G, G_regularized = make_green_functions(kb, a, config)
    theta = boundary_nodes(M)

    K_weighted = np.empty((M, M), dtype=np.complex128)
    for i, psi in enumerate(theta):
        for j, source_theta in enumerate(theta):
            K_weighted[i, j] = weighted_normal_kernel(
                psi,
                source_theta,
                epsilon,
                case,
                config,
                G,
                G_regularized,
            )

    # Full contour [0, 2*pi]:
    # 1/2 u_i = (2*pi/M) sum_j K^w_ij u_j
    # => A = I - (4*pi/M) K^w.
    A = np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted
    return A, theta, G


def smallest_singular_pair(A: np.ndarray) -> tuple[float, np.ndarray]:
    _, singular_values, Vh = np.linalg.svd(A, full_matrices=False)
    vector = Vh.conj().T[:, -1]
    vector /= np.linalg.norm(vector)

    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))

    return float(singular_values[-1]), vector


def sigma_minimum(
    kb: float,
    epsilon: float,
    case: CaseName,
    M: int,
    config: NumericalConfig,
) -> float:
    A, _, _ = assemble_matrix(kb, epsilon, case, M, config)
    sigma_min, _ = smallest_singular_pair(A)
    return sigma_min


def allowed_kb_interval(
    case: CaseName,
    config: NumericalConfig,
) -> tuple[float, float]:
    b = config.b
    kb_1 = np.sqrt((PI / (2.0 * b)) ** 2) * b
    kb_2 = np.sqrt((PI / b) ** 2) * b
    margin = 1.0e-9

    if case == "discrete_circle":
        return 1.0e-6, kb_1 - margin
    return kb_1 + margin, kb_2 - margin


def find_spectral_candidate(
    case: CaseName,
    epsilon: float,
    M: int,
    config: NumericalConfig,
) -> tuple[float, float, np.ndarray, np.ndarray, Callable[..., complex]]:
    """Locate a local minimum of sigma_min(A(k)) near the asymptotic value."""
    kb_asymptotic, _, _, _ = asymptotic_prediction(case, epsilon, config)
    band_left, band_right = allowed_kb_interval(case, config)

    left = max(band_left, kb_asymptotic - config.search_half_width)
    right = min(band_right, kb_asymptotic + config.search_half_width)
    if not left < right:
        raise RuntimeError("Invalid search interval.")

    scan = np.linspace(left, right, config.scan_points)
    scan_values = np.array(
        [sigma_minimum(kb, epsilon, case, M, config) for kb in scan]
    )
    index = int(np.argmin(scan_values))

    local_left = scan[max(0, index - 1)]
    local_right = scan[min(len(scan) - 1, index + 1)]
    if local_left == local_right:
        local_left, local_right = left, right

    tiny = np.finfo(float).tiny
    result = minimize_scalar(
        lambda kb: np.log10(max(sigma_minimum(kb, epsilon, case, M, config), tiny)),
        bounds=(local_left, local_right),
        method="bounded",
        options={"xatol": config.minimizer_xatol},
    )
    if not result.success:
        raise RuntimeError(f"Spectral minimization failed: {result.message}")

    kb_candidate = float(result.x)
    A, theta, G = assemble_matrix(kb_candidate, epsilon, case, M, config)
    sigma_min, boundary_vector = smallest_singular_pair(A)
    return kb_candidate, sigma_min, boundary_vector, theta, G


def weighted_kernel_at_field_point(
    x: float,
    y: float,
    theta: float,
    epsilon: float,
    case: CaseName,
    config: NumericalConfig,
    G: Callable[..., complex],
) -> complex:
    xi, eta, xi_p, eta_p, _, _ = obstacle_geometry(
        theta, epsilon, config.shape_beta, case
    )
    G_xi, G_eta = source_derivatives(
        G, x, y, xi, eta, config.finite_difference_step
    )
    return xi_p * G_eta - eta_p * G_xi


def reconstruct_field(
    x: float,
    y_values: np.ndarray,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    epsilon: float,
    case: CaseName,
    config: NumericalConfig,
    G: Callable[..., complex],
) -> np.ndarray:
    delta_theta = 2.0 * PI / len(theta)
    field = np.empty(len(y_values), dtype=np.complex128)

    for i, y in enumerate(y_values):
        kernel_values = np.array(
            [
                weighted_kernel_at_field_point(
                    x, float(y), float(t), epsilon, case, config, G
                )
                for t in theta
            ],
            dtype=np.complex128,
        )
        field[i] = delta_theta * np.dot(kernel_values, boundary_vector)

    return field


def first_propagating_component(
    case: CaseName,
    epsilon: float,
    boundary_vector: np.ndarray,
    theta: np.ndarray,
    G: Callable[..., complex],
    config: NumericalConfig,
) -> dict[str, float]:
    """Project the reconstructed BIC field onto the first open Dirichlet mode."""
    if case == "discrete_circle":
        return {
            "ratio_left": math.nan,
            "ratio_right": math.nan,
            "coefficient_left": math.nan,
            "coefficient_right": math.nan,
            "field_norm_left": math.nan,
            "field_norm_right": math.nan,
        }

    b = config.b
    a = obstacle_position(case, epsilon, config.shape_beta)
    x_far = config.projection_x_over_b * b

    nodes, weights = leggauss(config.projection_quadrature_points)
    y_lower = -b - a
    y_upper = b - a
    y_values = 0.5 * (y_upper - y_lower) * nodes + 0.5 * (y_upper + y_lower)
    physical_weights = 0.5 * (y_upper - y_lower) * weights

    # In obstacle-centred coordinates, the guide-centred coordinate is y+a.
    phi_1 = np.cos(PI * (y_values + a) / (2.0 * b))
    norm_phi_squared = float(np.sum(physical_weights * np.abs(phi_1) ** 2))

    diagnostics: dict[str, float] = {}
    for label, x in (("left", -x_far), ("right", x_far)):
        field = reconstruct_field(
            x,
            y_values,
            boundary_vector,
            theta,
            epsilon,
            case,
            config,
            G,
        )

        coefficient = np.sum(
            physical_weights * field * np.conjugate(phi_1)
        ) / norm_phi_squared
        projected_norm = abs(coefficient) * np.sqrt(norm_phi_squared)
        field_norm = float(np.sqrt(np.sum(physical_weights * np.abs(field) ** 2)))
        ratio = float(projected_norm / max(field_norm, 1.0e-30))

        diagnostics[f"ratio_{label}"] = ratio
        diagnostics[f"coefficient_{label}"] = float(abs(coefficient))
        diagnostics[f"field_norm_{label}"] = field_norm

    return diagnostics


def run_case(
    case: CaseName,
    epsilon: float,
    M: int,
    config: NumericalConfig,
) -> SpectralResult:
    kb_asymptotic, sigma_asymptotic, threshold, a = asymptotic_prediction(
        case, epsilon, config
    )

    kb_numerical, sigma_min, boundary_vector, theta, G = find_spectral_candidate(
        case, epsilon, M, config
    )

    k_numerical = kb_numerical / config.b
    distance_squared = threshold - k_numerical**2
    sigma_numerical = (
        float(np.sqrt(distance_squared)) if distance_squared >= 0.0 else math.nan
    )

    relative_error_k = abs(kb_numerical - kb_asymptotic) / abs(kb_asymptotic)
    relative_error_sigma = (
        abs(sigma_numerical - sigma_asymptotic) / abs(sigma_asymptotic)
        if np.isfinite(sigma_numerical) and sigma_asymptotic != 0.0
        else math.nan
    )

    propagation = first_propagating_component(
        case, epsilon, boundary_vector, theta, G, config
    )

    passes_singular_value = sigma_min <= config.singular_value_tolerance
    if case == "discrete_circle":
        passes_propagating_test = True
    else:
        passes_propagating_test = (
            propagation["ratio_left"] <= config.propagating_ratio_tolerance
            and propagation["ratio_right"] <= config.propagating_ratio_tolerance
        )

    return SpectralResult(
        case=case,
        epsilon=float(epsilon),
        a=a,
        M=M,
        kb_asymptotic=kb_asymptotic,
        kb_numerical=kb_numerical,
        sigma_asymptotic=sigma_asymptotic,
        sigma_numerical=sigma_numerical,
        relative_error_k=float(relative_error_k),
        relative_error_sigma=float(relative_error_sigma),
        sigma_min=sigma_min,
        svd_score=float(-2.0 * np.log10(max(sigma_min, np.finfo(float).tiny))),
        propagating_ratio_left=propagation["ratio_left"],
        propagating_ratio_right=propagation["ratio_right"],
        first_mode_coefficient_left=propagation["coefficient_left"],
        first_mode_coefficient_right=propagation["coefficient_right"],
        field_norm_left=propagation["field_norm_left"],
        field_norm_right=propagation["field_norm_right"],
        passes_singular_value=passes_singular_value,
        passes_propagating_test=passes_propagating_test,
        accepted_local_candidate=(passes_singular_value and passes_propagating_test),
    )


def run_refinement(
    case: CaseName,
    epsilon: float,
    config: NumericalConfig,
) -> list[RefinementResult]:
    results: list[RefinementResult] = []
    previous_kb = math.nan

    for M in config.refinement_M:
        result = run_case(case, epsilon, M, config)
        change = (
            abs(result.kb_numerical - previous_kb) / abs(result.kb_numerical)
            if np.isfinite(previous_kb)
            else math.nan
        )
        results.append(
            RefinementResult(
                case=case,
                epsilon=epsilon,
                M=M,
                kb_numerical=result.kb_numerical,
                sigma_min=result.sigma_min,
                relative_change_from_previous=change,
                propagating_ratio_left=result.propagating_ratio_left,
                propagating_ratio_right=result.propagating_ratio_right,
            )
        )
        previous_kb = result.kb_numerical

    return results


def write_csv(path: Path, rows: Iterable[object]) -> None:
    rows = list(rows)
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def latex_quantitative_table(results: list[SpectralResult]) -> str:
    lines = [
        r"\begin{tabular}{lrrrrrr}",
        r"\hline",
        r"Case & $\varepsilon$ & $kb_{\mathrm{BEM}}$ & $kb_{\mathrm{asym}}$ & $E_k$ & $\sigma_{\min}$ & $R_{\mathrm{prop}}$ \\",
        r"\hline",
    ]
    for row in results:
        propagation = (
            max(row.propagating_ratio_left, row.propagating_ratio_right)
            if row.case != "discrete_circle"
            else math.nan
        )
        propagation_text = f"{propagation:.3e}" if np.isfinite(propagation) else "--"
        lines.append(
            f"{row.case} & {row.epsilon:.3f} & {row.kb_numerical:.8f} & "
            f"{row.kb_asymptotic:.8f} & {row.relative_error_k:.3e} & "
            f"{row.sigma_min:.3e} & {propagation_text} \\\\" 
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return "\n".join(lines)


def latex_refinement_table(results: list[RefinementResult]) -> str:
    lines = [
        r"\begin{tabular}{lrrrrr}",
        r"\hline",
        r"Case & $M$ & $kb_M$ & $\Delta_M$ & $\sigma_{\min}$ & $R_{\mathrm{prop}}$ \\",
        r"\hline",
    ]
    for row in results:
        change = (
            f"{row.relative_change_from_previous:.3e}"
            if np.isfinite(row.relative_change_from_previous)
            else "--"
        )
        propagation = (
            max(row.propagating_ratio_left, row.propagating_ratio_right)
            if np.isfinite(row.propagating_ratio_left)
            else math.nan
        )
        propagation_text = f"{propagation:.3e}" if np.isfinite(propagation) else "--"
        lines.append(
            f"{row.case} & {row.M:d} & {row.kb_numerical:.8f} & {change} & "
            f"{row.sigma_min:.3e} & {propagation_text} \\\\" 
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return "\n".join(lines)


def plot_case_results(
    case: CaseName,
    results: list[SpectralResult],
    output_directory: Path,
    config: NumericalConfig,
) -> None:
    eps = np.array([row.epsilon for row in results])
    kb_num = np.array([row.kb_numerical for row in results])
    kb_asym = np.array([row.kb_asymptotic for row in results])
    error = np.array([row.relative_error_k for row in results])
    sigma_min = np.array([row.sigma_min for row in results])

    plt.figure(figsize=(8, 5))
    plt.plot(eps, kb_num, "o-", label=r"$kb_{\mathrm{BEM}}$")
    plt.plot(eps, kb_asym, "s--", label=r"$kb_{\mathrm{asym}}$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$kb$")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"{case}_kb_comparison.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.semilogy(eps, error, "o-")
    plt.axhline(
        y=config.relative_error_tolerance,
        linestyle="--",
        label="Reporting tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"Relative error in $kb$")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"{case}_relative_error.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.semilogy(eps, sigma_min, "o-")
    plt.axhline(
        y=config.singular_value_tolerance,
        linestyle="--",
        label=r"$\sigma_{\min}$ tolerance",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$\sigma_{\min}(A(k_*))$")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"{case}_sigma_min.png", dpi=200)
    plt.close()

    if case != "discrete_circle":
        ratio = np.array(
            [max(row.propagating_ratio_left, row.propagating_ratio_right) for row in results]
        )
        plt.figure(figsize=(8, 5))
        plt.semilogy(eps, ratio, "o-")
        plt.axhline(
            y=config.propagating_ratio_tolerance,
            linestyle="--",
            label="Propagating-component tolerance",
        )
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"First propagating-mode fraction")
        plt.grid(True, alpha=0.4)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_directory / f"{case}_propagating_component.png", dpi=200)
        plt.close()


def sampled_accuracy_limit(
    results: list[SpectralResult],
    config: NumericalConfig,
) -> float:
    """Largest sampled epsilon satisfying the declared accuracy criterion."""
    valid = [
        row.epsilon
        for row in results
        if row.accepted_local_candidate
        and np.isfinite(row.relative_error_sigma)
        and row.relative_error_sigma <= config.relative_error_tolerance
    ]
    return max(valid) if valid else math.nan


def print_summary(
    case: CaseName,
    results: list[SpectralResult],
    config: NumericalConfig,
) -> None:
    accepted = [row for row in results if row.accepted_local_candidate]
    print(f"\n=== {case} ===")
    print(
        f"Accepted local candidates: {len(accepted)}/{len(results)} "
        f"using sigma_min <= {config.singular_value_tolerance:g}"
    )

    if accepted:
        max_k_error = max(row.relative_error_k for row in accepted)
        sigma_errors = [
            row.relative_error_sigma
            for row in accepted
            if np.isfinite(row.relative_error_sigma)
        ]
        print(f"Maximum relative error in kb: {max_k_error:.6%}")
        if sigma_errors:
            print(
                "Maximum relative error in spectral distance: "
                f"{max(sigma_errors):.6%}"
            )

        if case != "discrete_circle":
            max_propagating = max(
                max(row.propagating_ratio_left, row.propagating_ratio_right)
                for row in accepted
            )
            print(
                "Maximum first-propagating-mode fraction: "
                f"{max_propagating:.3e}"
            )

    if case == "discrete_circle":
        epsilon_limit = sampled_accuracy_limit(results, config)
        print(
            "Largest sampled epsilon satisfying the declared "
            f"{config.relative_error_tolerance:.1%} error tolerance: "
            f"{epsilon_limit}"
        )


CONFIG = NumericalConfig()


def main() -> None:
    output_directory = Path(CONFIG.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    cases: tuple[CaseName, ...] = (
        "discrete_circle",
        "bic_x",
        "bic_y",
    )
    epsilon_values = np.linspace(0.01, 0.20, 10)
    all_results: list[SpectralResult] = []

    for case in cases:
        case_results: list[SpectralResult] = []
        for epsilon in epsilon_values:
            print(
                f"Running case={case}, epsilon={epsilon:.4f}, "
                f"M={CONFIG.boundary_points}"
            )
            try:
                row = run_case(case, float(epsilon), CONFIG.boundary_points, CONFIG)
            except Exception as exc:
                print(f"  FAILED: {exc}")
                continue

            case_results.append(row)
            all_results.append(row)

        write_csv(output_directory / f"{case}_quantitative_comparison.csv", case_results)
        (output_directory / f"{case}_quantitative_table.tex").write_text(
            latex_quantitative_table(case_results), encoding="utf-8"
        )
        plot_case_results(case, case_results, output_directory, CONFIG)
        print_summary(case, case_results, CONFIG)

    refinement_results: list[RefinementResult] = []
    for case in ("discrete_circle", "bic_x"):
        print(
            f"\nMesh refinement for case={case}, "
            f"epsilon={CONFIG.refinement_epsilon}"
        )
        refinement_results.extend(
            run_refinement(case, CONFIG.refinement_epsilon, CONFIG)
        )

    write_csv(output_directory / "mesh_refinement.csv", refinement_results)
    (output_directory / "mesh_refinement_table.tex").write_text(
        latex_refinement_table(refinement_results), encoding="utf-8"
    )

    write_csv(output_directory / "all_quantitative_results.csv", all_results)
    (output_directory / "all_quantitative_table.tex").write_text(
        latex_quantitative_table(all_results), encoding="utf-8"
    )

    print("\nLaTeX quantitative table:\n")
    print(latex_quantitative_table(all_results))
    print("\nLaTeX mesh-refinement table:\n")
    print(latex_refinement_table(refinement_results))
    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
