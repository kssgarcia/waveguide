"""
Numerical validation of Theorem 2.1 for the discrete trapped mode.

Problem
-------
Infinite 2-D waveguide of half-width b with Dirichlet conditions on the
waveguide walls and a Neumann condition on a circular obstacle of radius
r = epsilon centered at (0, a).

For the circular obstacle (R0 = 1), Theorem 2.1 predicts, for sufficiently
small epsilon and a > a*, exactly one discrete trapped eigenvalue in

    0 < k^2 < Lambda_1,

with

    k^2 = Lambda_1 - sigma^2,

and the leading asymptotic approximation

    sigma = epsilon^2 * pi^2/(4 b^3)
            * [pi sin^2(alpha/2) - (pi/2) cos^2(alpha/2)],
    alpha = pi a / b.

What this script checks numerically
-----------------------------------
1. The selected geometry satisfies the leading-order condition a > a0*.
2. A multiscale scan covers the complete discrete spectral band away from
   the numerical cutoff singularity.
3. Every resolved local minimum of sigma_min(A(k)) is refined with a finer
   BEM discretization.
4. A candidate mode is accepted only if sigma_min(A(k)) <= tolerance.
5. The uniqueness statement passes only when EXACTLY ONE candidate is
   accepted in the scanned band.
6. The unique numerical mode is compared with the asymptotic prediction.

This is a numerical validation, not a mathematical proof of uniqueness.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

# Keep the same import convention used by the existing project scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice


PI = np.pi


@dataclass(frozen=True)
class Config:
    # Geometry
    b: float = 1.0
    a: float = 0.6

    # A deliberately small/moderate asymptotic range.  epsilon=0.01 is avoided
    # here because the predicted eigenvalue is only O(1e-8) in kb from the
    # first cutoff, making the numerical localization unnecessarily delicate.
    epsilon_values: tuple[float, ...] = (0.05, 0.07, 0.09, 0.11)

    # BEM resolution: cheap discovery first, accurate refinement second.
    coarse_M: int = 16
    fine_M: int = 32
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # Whole-band multiscale discovery grid.
    linear_scan_points: int = 23
    logarithmic_scan_points: int = 28
    cutoff_delta_floor: float = 1.0e-7  # minimum Lambda_1-k^2 inspected
    low_k_margin: float = 1.0e-5

    # Local refinement / acceptance.
    minimizer_xatol: float = 2.0e-9
    candidate_merge_tolerance: float = 2.0e-7
    singular_value_tolerance: float = 1.0e-6
    relative_sigma_error_tolerance: float = 0.05

    # Output
    output_directory: str = "theorem_2_1_discrete_validation"


@dataclass
class Candidate:
    kb: float
    sigma_min: float
    accepted: bool


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
    relative_error_kb: float
    relative_error_sigma: float
    sigma_min: float
    resolved_local_minima: int
    accepted_mode_count: int
    accepted_kb_values: str
    unique_mode_verified: bool
    asymptotic_agreement_verified: bool
    theorem_numerically_verified: bool


CONFIG = Config()


def lambda_1(config: Config) -> float:
    return (PI / (2.0 * config.b)) ** 2


def kb_cutoff(config: Config) -> float:
    # sqrt(Lambda_1) * b = pi/2.
    return math.sqrt(lambda_1(config)) * config.b


def critical_height_leading_order(config: Config) -> float:
    """Leading-order a_0* for a circle with R0=1, mu=1 and S=pi."""
    mu = 1.0
    area = PI
    return (2.0 * config.b / PI) * math.atan(math.sqrt(area / (2.0 * PI * mu)))


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    """Return (kb_asymptotic, sigma_asymptotic)."""
    alpha = PI * config.a / config.b
    mu = 1.0
    area = PI

    bracket = (
        PI * mu * math.sin(alpha / 2.0) ** 2
        - 0.5 * area * math.cos(alpha / 2.0) ** 2
    )
    sigma = (
        epsilon**2
        * PI**2
        / (4.0 * config.b**3)
        * bracket
    )

    if sigma <= 0.0:
        raise ValueError(
            "The leading-order sigma is not positive for this geometry; "
            "Theorem 2.1 does not predict the discrete mode in this setup."
        )

    k_squared = lambda_1(config) - sigma**2
    if k_squared <= 0.0:
        raise ValueError("The asymptotic prediction produced non-positive k^2.")

    return math.sqrt(k_squared) * config.b, sigma


def boundary_nodes(M: int) -> np.ndarray:
    return (np.arange(M) + 0.5) * (2.0 * PI / M)


def circle_geometry(
    t: np.ndarray | float,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Circle of radius epsilon in obstacle-centred coordinates."""
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

    # Nodes are midpoint-collocation nodes. Equality identifies the diagonal.
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
    # hence
    #   A = I - (4pi/M) K^w.
    return np.eye(M, dtype=np.complex128) - (4.0 * PI / M) * K_weighted


def smallest_singular_value(kb: float, epsilon: float, M: int, config: Config) -> float:
    A = assemble_matrix(kb, epsilon, M, config)
    singular_values = np.linalg.svd(A, compute_uv=False)
    return float(singular_values[-1])


def build_multiscale_band_grid(epsilon: float, config: Config) -> np.ndarray:
    """
    Cover 0 < k^2 < Lambda_1 with two complementary grids.

    - A linear kb grid covers the lower/intermediate part of the band.
    - A logarithmic grid in delta = Lambda_1-k^2 gives high resolution near
      the first cutoff, where Theorem 2.1 places the discrete eigenvalue.
    - Extra points are clustered around the asymptotic spectral distance.
    """
    lam1 = lambda_1(config)
    kb1 = kb_cutoff(config)
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2

    linear = np.linspace(config.low_k_margin, kb1 - config.low_k_margin, config.linear_scan_points)

    delta_max = lam1 * (1.0 - 1.0e-8)
    deltas = np.geomspace(
        config.cutoff_delta_floor,
        delta_max,
        config.logarithmic_scan_points,
    )
    logarithmic = config.b * np.sqrt(np.maximum(lam1 - deltas, 0.0))

    # Explicitly resolve the scale predicted by the theorem.
    factors = np.array(
        [0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00, 1.15, 1.35, 1.70, 2.5, 5.0, 10.0]
    )
    local_deltas = delta_asym * factors
    local_deltas = local_deltas[(local_deltas > config.cutoff_delta_floor) & (local_deltas < lam1)]
    local = config.b * np.sqrt(lam1 - local_deltas)

    grid = np.concatenate([linear, logarithmic, local, [kb_asym]])
    grid = grid[(grid > 0.0) & (grid < kb1)]
    return np.unique(np.sort(grid))


def local_minimum_brackets(scan: np.ndarray, values: np.ndarray) -> list[tuple[float, float]]:
    """Bracket every resolved sampled minimum, including numerical band edges."""
    brackets: list[tuple[float, float]] = []
    if len(scan) < 2:
        return brackets

    if values[0] <= values[1]:
        brackets.append((float(scan[0]), float(scan[1])))

    for i in range(1, len(scan) - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            brackets.append((float(scan[i - 1]), float(scan[i + 1])))

    if values[-1] <= values[-2]:
        brackets.append((float(scan[-2]), float(scan[-1])))

    return brackets


def merge_candidates(candidates: list[Candidate], config: Config) -> list[Candidate]:
    if not candidates:
        return []

    candidates = sorted(candidates, key=lambda c: c.kb)
    merged = [candidates[0]]
    for candidate in candidates[1:]:
        previous = merged[-1]
        if abs(candidate.kb - previous.kb) <= config.candidate_merge_tolerance:
            if candidate.sigma_min < previous.sigma_min:
                merged[-1] = candidate
        else:
            merged.append(candidate)
    return merged


def discover_and_refine_all_modes(
    epsilon: float,
    config: Config,
) -> tuple[list[Candidate], np.ndarray, np.ndarray, int]:
    """Discover candidates coarsely across the whole band, then refine with fine_M."""
    scan = build_multiscale_band_grid(epsilon, config)

    print(f"  coarse whole-band scan: {len(scan)} spectral points (M={config.coarse_M})")
    values = np.empty(len(scan), dtype=float)
    for i, kb in enumerate(scan):
        values[i] = smallest_singular_value(float(kb), epsilon, config.coarse_M, config)
        if (i + 1) % 15 == 0 or i + 1 == len(scan):
            print(f"    {i + 1:>3}/{len(scan)}")

    brackets = local_minimum_brackets(scan, values)

    # The asymptotic mode is expected extremely close to Lambda_1. Add one
    # explicit bracket around k_asym so the theorem-predicted candidate is not
    # lost merely because the coarse discovery grid is sparse at that point.
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
    delta_asym = sigma_asym**2
    lam1 = lambda_1(config)
    delta_left = min(lam1 * 0.999, delta_asym * 3.0)
    delta_right = max(config.cutoff_delta_floor, delta_asym * 0.30)
    k_left = config.b * math.sqrt(max(lam1 - delta_left, 0.0))
    k_right = config.b * math.sqrt(max(lam1 - delta_right, 0.0))
    if k_left < kb_asym < k_right:
        brackets.append((k_left, k_right))

    print(f"  resolved coarse local-minimum brackets: {len(brackets)}")

    tiny = np.finfo(float).tiny
    refined: list[Candidate] = []

    # Fine BEM is used only inside intervals that coarse discovery identified.
    for j, (left, right) in enumerate(brackets, start=1):
        print(f"  refining candidate {j}/{len(brackets)} with M={config.fine_M}")
        result = minimize_scalar(
            lambda kb: math.log10(
                max(smallest_singular_value(float(kb), epsilon, config.fine_M, config), tiny)
            ),
            bounds=(left, right),
            method="bounded",
            options={"xatol": config.minimizer_xatol},
        )
        if not result.success:
            continue

        kb_candidate = float(result.x)
        sigma_min = smallest_singular_value(kb_candidate, epsilon, config.fine_M, config)
        refined.append(
            Candidate(
                kb=kb_candidate,
                sigma_min=sigma_min,
                accepted=sigma_min <= config.singular_value_tolerance,
            )
        )

    return merge_candidates(refined, config), scan, values, len(brackets)


def validate_epsilon(epsilon: float, config: Config) -> tuple[ValidationResult, np.ndarray, np.ndarray]:
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
    candidates, scan, scan_values, minima_count = discover_and_refine_all_modes(epsilon, config)
    accepted = [candidate for candidate in candidates if candidate.accepted]

    unique = len(accepted) == 1
    accepted_text = ", ".join(f"{c.kb:.12f}" for c in accepted) if accepted else "none"

    if unique:
        mode = accepted[0]
        kb_num = mode.kb
        sigma_min = mode.sigma_min
        k_num = kb_num / config.b
        spectral_distance = lambda_1(config) - k_num**2
        sigma_num = math.sqrt(max(spectral_distance, 0.0))
        error_kb = abs(kb_num - kb_asym) / abs(kb_asym)
        error_sigma = abs(sigma_num - sigma_asym) / abs(sigma_asym)
        agreement = error_sigma <= config.relative_sigma_error_tolerance
    else:
        kb_num = math.nan
        sigma_min = math.nan
        sigma_num = math.nan
        error_kb = math.nan
        error_sigma = math.nan
        agreement = False

    result = ValidationResult(
        epsilon=epsilon,
        a=config.a,
        lambda_1=lambda_1(config),
        kb_cutoff=kb_cutoff(config),
        kb_asymptotic=kb_asym,
        kb_numerical=kb_num,
        sigma_asymptotic=sigma_asym,
        sigma_numerical=sigma_num,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        sigma_min=sigma_min,
        resolved_local_minima=minima_count,
        accepted_mode_count=len(accepted),
        accepted_kb_values=accepted_text,
        unique_mode_verified=unique,
        asymptotic_agreement_verified=agreement,
        theorem_numerically_verified=(unique and agreement),
    )
    return result, scan, scan_values


def write_csv(results: list[ValidationResult], output_directory: Path) -> None:
    path = output_directory / "theorem_2_1_validation.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0])))
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)


def plot_scan(
    epsilon: float,
    scan: np.ndarray,
    values: np.ndarray,
    result: ValidationResult,
    output_directory: Path,
) -> None:
    plt.figure(figsize=(8, 5))
    plt.semilogy(scan, values, "o-", markersize=3, label=r"coarse $\sigma_{\min}(A)$")
    plt.axvline(result.kb_asymptotic, linestyle="--", label=r"$kb_{\mathrm{asym}}$")
    if np.isfinite(result.kb_numerical):
        plt.axvline(result.kb_numerical, linestyle=":", label=r"$kb_{\mathrm{BEM}}$")
    plt.xlabel(r"$kb$")
    plt.ylabel(r"$\sigma_{\min}(A(k))$")
    plt.title(rf"Theorem 2.1, $\varepsilon={epsilon:.2f}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / f"scan_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()


def plot_summary(results: list[ValidationResult], output_directory: Path) -> None:
    eps = np.array([r.epsilon for r in results])
    kb_asym = np.array([r.kb_asymptotic for r in results])
    kb_num = np.array([r.kb_numerical for r in results])

    plt.figure(figsize=(8, 5))
    plt.plot(eps, kb_asym, "s--", label=r"$kb_{\mathrm{asym}}$")
    plt.plot(eps, kb_num, "o-", label=r"$kb_{\mathrm{BEM}}$")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$kb$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "kb_asymptotic_vs_bem.png", dpi=180)
    plt.close()

    error_sigma = np.array([r.relative_error_sigma for r in results])
    plt.figure(figsize=(8, 5))
    plt.plot(eps, error_sigma, "o-")
    plt.axhline(CONFIG.relative_sigma_error_tolerance, linestyle="--", label="5% tolerance")
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"relative error in $\sigma$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "relative_error_sigma.png", dpi=180)
    plt.close()


def print_result(result: ValidationResult, config: Config) -> None:
    status = "PASS" if result.theorem_numerically_verified else "FAIL"
    print(f"\n  epsilon = {result.epsilon:.3f}")
    print(f"  accepted modes in 0 < k^2 < Lambda_1: {result.accepted_mode_count}")
    print(f"  accepted kb: {result.accepted_kb_values}")
    print(f"  uniqueness: {'PASS' if result.unique_mode_verified else 'FAIL'}")

    if result.unique_mode_verified:
        print(f"  kb asymptotic = {result.kb_asymptotic:.12f}")
        print(f"  kb BEM        = {result.kb_numerical:.12f}")
        print(f"  sigma asym    = {result.sigma_asymptotic:.8e}")
        print(f"  sigma BEM     = {result.sigma_numerical:.8e}")
        print(f"  sigma_min(A)  = {result.sigma_min:.3e}")
        print(f"  relative sigma error = {result.relative_error_sigma:.3%}")
        print(
            f"  asymptotic agreement (<={config.relative_sigma_error_tolerance:.1%}): "
            f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
        )

    print(f"  THEOREM 2.1 NUMERICAL CHECK: {status}")


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    a0_star = critical_height_leading_order(config)
    print("=== Theorem 2.1: discrete trapped mode ===")
    print(f"b = {config.b}")
    print(f"a = {config.a}")
    print(f"leading-order a0* = {a0_star:.12f}")
    print(f"geometric existence condition a > a0*: {'PASS' if config.a > a0_star else 'FAIL'}")
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"sqrt(Lambda_1) b = {kb_cutoff(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")

    if config.a <= a0_star:
        raise RuntimeError("Selected a does not satisfy the leading-order existence condition a > a0*.")

    results: list[ValidationResult] = []
    for epsilon in config.epsilon_values:
        print(f"\n--- epsilon={epsilon:.3f} ---")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
        gap = kb_cutoff(config) - kb_asym
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted cutoff gap = {gap:.3e}")
        print(f"  predicted sigma = {sigma_asym:.3e}")

        result, scan, scan_values = validate_epsilon(epsilon, config)
        results.append(result)
        print_result(result, config)
        plot_scan(epsilon, scan, scan_values, result, output_directory)

    write_csv(results, output_directory)
    plot_summary(results, output_directory)

    passed = sum(result.theorem_numerically_verified for result in results)
    unique_passed = sum(result.unique_mode_verified for result in results)

    print("\n=== FINAL SUMMARY ===")
    print("Numerical interval checked: 0 < k^2 < Lambda_1 (excluding only a tiny numerical cutoff layer).")
    print(f"Exactly-one-mode checks passed: {unique_passed}/{len(results)}")
    print(
        f"Uniqueness + <= {config.relative_sigma_error_tolerance:.1%} asymptotic sigma error: "
        f"{passed}/{len(results)}"
    )

    for result in results:
        print(
            f"  epsilon={result.epsilon:.3f}: "
            f"modes={result.accepted_mode_count}, "
            f"unique={'PASS' if result.unique_mode_verified else 'FAIL'}, "
            f"asymptotic={'PASS' if result.asymptotic_agreement_verified else 'FAIL'}, "
            f"overall={'PASS' if result.theorem_numerically_verified else 'FAIL'}"
        )

    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
