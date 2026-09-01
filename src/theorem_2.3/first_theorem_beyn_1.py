"""
Theorem 2.1 — discrete trapped mode discovered with Beyn's contour method (v2).

This script keeps the same BEM discretization used by the refined Theorem 2.1
validator, but replaces the exhaustive whole-band sigma_min scan with a
nonlinear contour eigensolver of Beyn type.

Problem
-------
Find all real trapped-mode wavenumbers in

    0 < k^2 < Lambda_1,

for the circular obstacle centered at height a.  The nonlinear BEM eigenproblem
is

    A(kb) u = 0,

where kb = k*b and A(kb) is the same boundary-integral matrix used previously.

Beyn discovery
--------------
For a closed contour Gamma in the complex kb-plane and a probing matrix V,
compute

    S0 = (1/(2*pi*i)) int_Gamma A(z)^(-1) V dz,
    S1 = (1/(2*pi*i)) int_Gamma z A(z)^(-1) V dz.

If S0 = U Sigma W*, the reduced matrix

    B = U_r^* S1 W_r Sigma_r^(-1)

has eigenvalues approximating all nonlinear eigenvalues enclosed by Gamma.

Workflow
--------
1. Run a Beyn quadrature-convergence study on one contour enclosing the
   discrete spectral band.
2. Estimate rank(S0) using both a relative threshold and the dominant spectral
   gap, instead of a single hard threshold.
3. Keep only near-real eigenvalues that are actually inside the contour and
   the physical discrete band.
4. Locally refine each final Beyn candidate by minimizing sigma_min(A(kb)) on
   the real axis. This is verification, not global discovery.
5. Repeat the local verification for several BEM meshes M.
6. Count the persistent resolved modes and compare the branch nearest the
   theorem prediction with the leading asymptotic formula.

The script therefore separates:

    GLOBAL DISCOVERY : Beyn contour method
    LOCAL VALIDATION : SVD / sigma_min(A)
    ASYMPTOTIC TEST  : sigma_BEM versus sigma_asym

Important
---------
* The contour must not cross k=0 or the threshold sqrt(Lambda_1)*b, because the
  Green function changes analytic structure at thresholds.
* Beyn requires A(z) for complex z.  A startup smoke test checks whether the
  project's lattice_sums / Green-function implementation accepts complex kb.
* The default contour is kept farther from the first cutoff than in v1.
* The code reports Beyn convergence explicitly; failure of Beyn discovery is
  not interpreted as failure of the theorem.
* This is a numerical validation, not a mathematical proof of uniqueness.
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

# Same project import convention as the existing scripts.
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import lattice_sums as lattice


PI = np.pi


@dataclass(frozen=True)
class Config:
    # Geometry.
    b: float = 1.0
    a: float = 0.6
    epsilon_values: tuple[float, ...] = (0.05, 0.07, 0.09, 0.11)

    # Green function / BEM.
    lattice_terms: int = 200
    harmonic_order: int = 20
    finite_difference_step: float = 1.0e-6

    # ------------------------------------------------------------------
    # Beyn global discovery.
    # ------------------------------------------------------------------
    # BEM order used only for contour discovery. Local candidates are then
    # validated at the higher refinement_M values below.
    beyn_M: int = 24

    # Quadrature-convergence study. The last level is used for final global
    # discovery; the previous levels diagnose whether the contour moments and
    # eigenvalue count are converging.
    beyn_quadrature_levels: tuple[int, ...] = (96, 192, 384)

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

    # The contour is an ellipse enclosing the real discrete band. It stays a
    # small positive distance away from k=0 and from the first cutoff.
    beyn_low_k_margin: float = 1.0e-4
    beyn_cutoff_margin: float = 1.0e-5
    beyn_imag_half_height: float = 3.0e-2

    # A true discrete trapped mode is real. Beyn may return a small imaginary
    # part from quadrature/discretization error. Candidates below this threshold
    # are projected to the real axis and locally verified with SVD.
    beyn_real_axis_tolerance: float = 5.0e-5

    # Cluster duplicate near-real Beyn eigenvalues if numerical multiplicity
    # produces nearly identical values.
    beyn_cluster_tolerance: float = 5.0e-5

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
    near_singular_tolerance: float = 1.0e-4
    mesh_sigma_relative_tolerance: float = 0.03

    # Asymptotic accuracy is intentionally separate from mode existence.
    relative_sigma_error_tolerance: float = 0.05

    output_directory: str = "theorem_2_1_beyn_v2_validation"


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
    accepted_for_real_refinement: bool


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
    primary_candidate_real: float
    primary_candidate_imag: float
    max_candidate_shift_from_previous: float
    rank_stable_from_previous: bool


@dataclass
class ModeRefinementRow:
    epsilon: float
    candidate_index: int
    beyn_kb_real: float
    beyn_kb_imag: float
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
class ModeResult:
    epsilon: float
    candidate_index: int
    beyn_kb_real: float
    beyn_kb_imag: float
    kb_numerical: float
    sigma_numerical: float
    sigma_min_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    resolved: bool


@dataclass
class ValidationResult:
    epsilon: float
    a: float
    lambda_1: float
    kb_cutoff: float
    kb_asymptotic: float
    sigma_asymptotic: float
    beyn_final_quadrature_points: int
    beyn_estimated_rank: int
    beyn_rank_stable: bool
    beyn_near_real_candidates: int
    resolved_mode_count: int
    kb_numerical: float
    sigma_numerical: float
    sigma_min_final: float
    final_drop_factor: float
    final_relative_mesh_change: float
    relative_error_kb: float
    relative_error_sigma: float
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
    return (2.0 * config.b / PI) * math.atan(
        math.sqrt(area / (2.0 * PI * mu))
    )


def asymptotic_prediction(epsilon: float, config: Config) -> tuple[float, float]:
    """Return (kb_asymptotic, sigma_asymptotic)."""
    alpha = PI * config.a / config.b
    mu = 1.0
    area = PI

    bracket = (
        PI * mu * math.sin(alpha / 2.0) ** 2
        - 0.5 * area * math.cos(alpha / 2.0) ** 2
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


def smallest_singular_value(
    kb: float,
    epsilon: float,
    M: int,
    config: Config,
) -> float:
    A = assemble_matrix(complex(kb, 0.0), epsilon, M, config)
    singular_values = np.linalg.svd(A, compute_uv=False)
    return float(singular_values[-1])


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
    Estimate rank(S0) using a threshold estimate plus a dominant-gap estimate.

    Returns
    -------
    estimated_rank, threshold_rank, gap_rank, selected_gap_ratio

    The threshold rank is the conventional count s_j > tol*s_1.  We then look
    only at gaps up to that provisional rank and choose the largest ratio
    s_j/s_{j+1}.  If that gap is at least beyn_rank_gap_threshold, the gap rank
    is preferred.  This is intentionally transparent and is reported alongside
    quadrature convergence; it is not treated as a mathematical proof.
    """
    s = np.asarray(singular_values, dtype=float)
    if len(s) == 0 or not np.isfinite(s[0]) or s[0] <= 0.0:
        return 0, 0, 0, math.nan

    threshold = max(
        config.beyn_rank_absolute_tolerance,
        config.beyn_rank_relative_tolerance * s[0],
    )
    threshold_rank = int(np.sum(s > threshold))
    threshold_rank = max(1, min(threshold_rank, len(s)))

    if len(s) == 1:
        return 1, threshold_rank, 1, math.inf

    # Search the gap immediately after any singular value provisionally kept by
    # the threshold test. This avoids choosing arbitrary ratios deep in the
    # numerical tail.
    n_gaps = min(threshold_rank, len(s) - 1)
    floor = max(config.beyn_rank_absolute_tolerance, np.finfo(float).tiny)
    ratios = s[:n_gaps] / np.maximum(s[1 : n_gaps + 1], floor)
    best_index = int(np.argmax(ratios))
    selected_gap = float(ratios[best_index])
    gap_rank = best_index + 1

    if selected_gap >= config.beyn_rank_gap_threshold:
        rank = gap_rank
    else:
        rank = threshold_rank

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

    for idx, eig in enumerate(raw_eigenvalues):
        eig = complex(eig)
        inside = point_inside_ellipse(eig, config)
        inside_real_band = bool(0.0 < eig.real < cutoff)
        inside_contour_real_span = bool(left < eig.real < right)
        near_real = bool(abs(eig.imag) <= config.beyn_real_axis_tolerance)
        accepted = bool(
            inside and inside_real_band and inside_contour_real_span and near_real
        )
        accepted_count += int(accepted)
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
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
]:
    """
    Run all configured quadrature levels and return the highest-resolution run.

    The same probing matrix is used at every level because probing_matrix is
    seeded deterministically. Therefore S0 matrices can be compared directly.
    """
    all_diagnostics: list[BeynDiagnostics] = []
    all_eigen_rows: list[BeynEigenvalueRow] = []
    convergence_rows: list[BeynConvergenceRow] = []

    previous_S0: np.ndarray | None = None
    previous_candidates: list[complex] = []
    previous_rank: int | None = None

    final_raw = np.array([], dtype=np.complex128)
    final_s = np.array([], dtype=float)
    final_diag: BeynDiagnostics | None = None
    final_rows: list[BeynEigenvalueRow] = []
    final_candidates: list[complex] = []

    print("  Beyn quadrature-convergence study:")

    for Nq in config.beyn_quadrature_levels:
        raw, s, diag, rows, S0 = beyn_discover(epsilon, Nq, config)
        candidates = cluster_real_candidates(raw, rows, config)

        if previous_S0 is None:
            s0_change = math.nan
        else:
            s0_change = float(
                np.linalg.norm(S0 - previous_S0)
                / max(np.linalg.norm(S0), 1.0e-30)
            )

        shift = candidate_set_shift(previous_candidates, candidates)
        rank_stable = bool(previous_rank is not None and previous_rank == diag.estimated_rank)

        primary = min(candidates, key=lambda z: abs(z.imag)) if candidates else complex(math.nan, math.nan)
        convergence_rows.append(
            BeynConvergenceRow(
                epsilon=epsilon,
                quadrature_points=Nq,
                estimated_rank=diag.estimated_rank,
                threshold_rank=diag.threshold_rank,
                gap_rank=diag.gap_rank,
                selected_gap_ratio=diag.selected_gap_ratio,
                s0_sv1=float(s[0]) if len(s) > 0 else math.nan,
                s0_sv2=float(s[1]) if len(s) > 1 else math.nan,
                s0_sv3=float(s[2]) if len(s) > 2 else math.nan,
                s0_relative_change_from_previous=s0_change,
                accepted_candidates=len(candidates),
                primary_candidate_real=float(primary.real),
                primary_candidate_imag=float(primary.imag),
                max_candidate_shift_from_previous=shift,
                rank_stable_from_previous=rank_stable,
            )
        )

        change_text = f"{s0_change:.3e}" if np.isfinite(s0_change) else "--"
        shift_text = f"{shift:.3e}" if np.isfinite(shift) else "--"
        print(
            f"    Nq={Nq}: rank={diag.estimated_rank}, candidates={len(candidates)}, "
            f"S0-change={change_text}, candidate-shift={shift_text}"
        )

        all_diagnostics.append(diag)
        all_eigen_rows.extend(rows)

        previous_S0 = S0
        previous_candidates = candidates
        previous_rank = diag.estimated_rank

        final_raw = raw
        final_s = s
        final_diag = diag
        final_rows = rows
        final_candidates = candidates

    assert final_diag is not None
    return (
        final_raw,
        final_s,
        final_diag,
        final_rows,
        final_candidates,
        all_diagnostics,
        all_eigen_rows,
        convergence_rows,
    )


# ---------------------------------------------------------------------------
# Local real-axis verification of Beyn candidates
# ---------------------------------------------------------------------------


def candidate_brackets(candidates: list[complex], config: Config) -> list[tuple[float, float]]:
    """
    Build non-overlapping local real-axis brackets around Beyn candidates.

    Midpoints between neighboring Beyn candidates prevent two very close modes
    from collapsing into the same local minimization.
    """
    if not candidates:
        return []

    xs = np.array(sorted(z.real for z in candidates), dtype=float)
    left_band = config.beyn_low_k_margin
    right_band = kb_cutoff(config) - config.beyn_cutoff_margin
    brackets: list[tuple[float, float]] = []

    for i, x in enumerate(xs):
        neighbor_left = 0.5 * (xs[i - 1] + x) if i > 0 else left_band
        neighbor_right = 0.5 * (x + xs[i + 1]) if i + 1 < len(xs) else right_band

        left = max(left_band, neighbor_left, x - config.local_refine_half_width)
        right = min(right_band, neighbor_right, x + config.local_refine_half_width)

        if not left < x < right:
            # Near a threshold, keep a very small but valid bracket.
            tiny_width = max(1.0e-8, 0.1 * config.local_refine_half_width)
            left = max(left_band, x - tiny_width)
            right = min(right_band, x + tiny_width)

        if not left < right:
            raise RuntimeError(f"Could not build local bracket around Beyn candidate kb={x}.")
        brackets.append((float(left), float(right)))

    return brackets


def refine_candidate_for_M(
    epsilon: float,
    M: int,
    left: float,
    right: float,
    config: Config,
) -> tuple[float, float, float, float, float, bool]:
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
    sv = smallest_singular_value(kb, epsilon, M, config)
    drop = min(left_value, right_value) / max(sv, tiny)

    width = right - left
    edge_margin = 0.005 * width
    interior = (kb > left + edge_margin) and (kb < right - edge_margin)
    return kb, sv, left_value, right_value, drop, interior


def run_candidate_refinement(
    epsilon: float,
    candidate_index: int,
    beyn_candidate: complex,
    bracket: tuple[float, float],
    config: Config,
) -> tuple[list[ModeRefinementRow], ModeResult]:
    rows: list[ModeRefinementRow] = []
    previous_sigma = math.nan
    left, right = bracket

    print(
        f"  candidate {candidate_index}: Beyn kb={beyn_candidate.real:.12f} "
        f"{beyn_candidate.imag:+.3e}i, local bracket=[{left:.12f}, {right:.12f}]"
    )

    for M in config.refinement_M:
        kb, sv, left_value, right_value, drop, interior = refine_candidate_for_M(
            epsilon, M, left, right, config
        )
        sigma_bem = sigma_from_kb(kb, config)
        change = (
            abs(sigma_bem - previous_sigma) / max(abs(sigma_bem), 1.0e-30)
            if np.isfinite(previous_sigma)
            else math.nan
        )

        row = ModeRefinementRow(
            epsilon=epsilon,
            candidate_index=candidate_index,
            beyn_kb_real=float(beyn_candidate.real),
            beyn_kb_imag=float(beyn_candidate.imag),
            M=M,
            kb=kb,
            sigma_bem=sigma_bem,
            sigma_min=sv,
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
            f"sv_min={sv:.3e}, drop={drop:.2e}, mesh_change={change_text}, "
            f"interior={'yes' if interior else 'no'}"
        )

    final = rows[-1]
    previous = rows[-2] if len(rows) >= 2 else rows[-1]
    final_change = abs(final.sigma_bem - previous.sigma_bem) / max(
        abs(final.sigma_bem), 1.0e-30
    )

    resolved = bool(
        final.minimum_is_interior
        and final.sigma_min <= config.near_singular_tolerance
        and final.drop_factor >= config.minimum_drop_factor
        and final_change <= config.mesh_sigma_relative_tolerance
    )

    result = ModeResult(
        epsilon=epsilon,
        candidate_index=candidate_index,
        beyn_kb_real=float(beyn_candidate.real),
        beyn_kb_imag=float(beyn_candidate.imag),
        kb_numerical=final.kb,
        sigma_numerical=final.sigma_bem,
        sigma_min_final=final.sigma_min,
        final_drop_factor=final.drop_factor,
        final_relative_mesh_change=final_change,
        resolved=resolved,
    )
    print(f"    resolved={'YES' if resolved else 'no'}")
    return rows, result


# ---------------------------------------------------------------------------
# Validation / reporting
# ---------------------------------------------------------------------------


def validate_epsilon(
    epsilon: float,
    config: Config,
) -> tuple[
    ValidationResult,
    list[BeynDiagnostics],
    list[BeynEigenvalueRow],
    list[BeynConvergenceRow],
    list[ModeRefinementRow],
    list[ModeResult],
    np.ndarray,
    list[BeynEigenvalueRow],
]:
    kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)

    (
        final_raw,
        final_s,
        final_diag,
        final_eigen_rows,
        candidates,
        all_beyn_diagnostics,
        all_eigen_rows,
        convergence_rows,
    ) = run_beyn_convergence_study(epsilon, config)

    brackets = candidate_brackets(candidates, config)
    print(f"  final clustered near-real Beyn candidates = {len(candidates)}")

    all_refinement: list[ModeRefinementRow] = []
    mode_results: list[ModeResult] = []

    for i, (candidate, bracket) in enumerate(zip(candidates, brackets, strict=True), start=1):
        rows, mode = run_candidate_refinement(epsilon, i, candidate, bracket, config)
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

    # This is deliberately stricter than simply finding one local SVD minimum:
    # global one-mode support requires a stable final Beyn rank of one plus one
    # locally resolved mode. Candidate positional stability is reported but is
    # not required because local SVD is the position-certification stage.
    one_mode_count_supported = bool(
        final_diag.estimated_rank == 1
        and rank_stable
        and len(candidates) == 1
        and len(resolved) == 1
    )

    if resolved:
        expected = min(resolved, key=lambda mode: abs(mode.kb_numerical - kb_asym))
        kb_num = expected.kb_numerical
        sigma_num = expected.sigma_numerical
        sv_final = expected.sigma_min_final
        drop_final = expected.final_drop_factor
        mesh_change = expected.final_relative_mesh_change
        error_kb = abs(kb_num - kb_asym) / abs(kb_asym)
        error_sigma = abs(sigma_num - sigma_asym) / abs(sigma_asym)
        asymptotic_ok = error_sigma <= config.relative_sigma_error_tolerance
    else:
        kb_num = math.nan
        sigma_num = math.nan
        sv_final = math.nan
        drop_final = math.nan
        mesh_change = math.nan
        error_kb = math.nan
        error_sigma = math.nan
        asymptotic_ok = False

    if not candidate_stable and candidates:
        print(
            "  NOTE: final Beyn candidate position is still moving by more than "
            f"{config.beyn_candidate_convergence_tolerance:.1e} between the last "
            "two Nq levels. Local SVD refinement remains the trusted position."
        )

    result = ValidationResult(
        epsilon=epsilon,
        a=config.a,
        lambda_1=lambda_1(config),
        kb_cutoff=kb_cutoff(config),
        kb_asymptotic=kb_asym,
        sigma_asymptotic=sigma_asym,
        beyn_final_quadrature_points=final_diag.quadrature_points,
        beyn_estimated_rank=final_diag.estimated_rank,
        beyn_rank_stable=rank_stable,
        beyn_near_real_candidates=len(candidates),
        resolved_mode_count=len(resolved),
        kb_numerical=kb_num,
        sigma_numerical=sigma_num,
        sigma_min_final=sv_final,
        final_drop_factor=drop_final,
        final_relative_mesh_change=mesh_change,
        relative_error_kb=error_kb,
        relative_error_sigma=error_sigma,
        unique_mode_verified=one_mode_count_supported,
        asymptotic_agreement_verified=asymptotic_ok,
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
    )


def write_dataclass_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


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
    plt.savefig(output_directory / f"beyn_s0_convergence_epsilon_{epsilon:.2f}.png", dpi=180)
    plt.close()

    real_candidates = np.array([row.primary_candidate_real for row in rows], dtype=float)
    finite = np.isfinite(real_candidates)
    if np.any(finite):
        plt.figure(figsize=(7, 4.5))
        plt.plot(nq[finite], real_candidates[finite], "o-")
        plt.axhline(kb_cutoff(CONFIG), linestyle="--", label=r"$\sqrt{\Lambda_1}b$")
        plt.xlabel(r"Beyn quadrature points $N_q$")
        plt.ylabel(r"primary Beyn candidate $kb$")
        plt.title(rf"Beyn candidate convergence, $\varepsilon={epsilon:.2f}$")
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
        kb = np.array([row.kb for row in subset])

        plt.figure(figsize=(7, 4.5))
        plt.semilogy(M, sv, "o-")
        plt.xlabel(r"$M$")
        plt.ylabel(r"$\sigma_{\min}(A(k_*))$")
        plt.title(
            rf"Beyn candidate {candidate_index}, $\varepsilon={epsilon:.2f}$: SVD refinement"
        )
        plt.tight_layout()
        plt.savefig(
            output_directory / f"candidate_{candidate_index}_svd_epsilon_{epsilon:.2f}.png",
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
            output_directory / f"candidate_{candidate_index}_kb_epsilon_{epsilon:.2f}.png",
            dpi=180,
        )
        plt.close()


def plot_summary(results: list[ValidationResult], output_directory: Path, config: Config) -> None:
    eps = np.array([row.epsilon for row in results])
    counts = np.array([row.resolved_mode_count for row in results])

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
        plt.axhline(config.relative_sigma_error_tolerance, linestyle="--", label="5% tolerance")
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"relative error in $\sigma$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_directory / "relative_error_sigma.png", dpi=180)
        plt.close()


def print_result(result: ValidationResult, config: Config) -> None:
    print("\n  --- Beyn-v2 + local-SVD validation result ---")
    print(f"  final Beyn quadrature Nq       = {result.beyn_final_quadrature_points}")
    print(f"  final Beyn estimated rank      = {result.beyn_estimated_rank}")
    print(f"  Beyn rank stable (last 2 Nq)   = {'YES' if result.beyn_rank_stable else 'no'}")
    print(f"  final near-real candidates     = {result.beyn_near_real_candidates}")
    print(f"  locally resolved modes         = {result.resolved_mode_count}")
    print(
        "  one-mode count supported      = "
        f"{'PASS' if result.unique_mode_verified else 'FAIL'}"
    )
    print(f"  kb asymptotic                  = {result.kb_asymptotic:.12f}")
    print(f"  kb BEM                         = {result.kb_numerical:.12f}")
    print(f"  sigma asym                     = {result.sigma_asymptotic:.8e}")
    print(f"  sigma BEM                      = {result.sigma_numerical:.8e}")
    print(f"  final sigma_min(A)             = {result.sigma_min_final:.3e}")
    print(f"  final minimum drop             = {result.final_drop_factor:.2e}")
    print(f"  final mesh change in sigma     = {result.final_relative_mesh_change:.3%}")
    print(f"  relative sigma error           = {result.relative_error_sigma:.3%}")
    print(
        f"  asymptotic accuracy <= {config.relative_sigma_error_tolerance:.1%}: "
        f"{'PASS' if result.asymptotic_agreement_verified else 'FAIL'}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    config = CONFIG
    output_directory = Path(config.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    a0_star = critical_height_leading_order(config)
    left, right, center, rx, ry = contour_geometry(config)

    print("=== Theorem 2.1 validation with Beyn contour discovery (v2) ===")
    print(f"b = {config.b}")
    print(f"a = {config.a}")
    print(f"leading-order a0* = {a0_star:.12f}")
    print(
        f"geometric condition a > a0*: "
        f"{'PASS' if config.a > a0_star else 'FAIL'}"
    )
    print(f"Lambda_1 = {lambda_1(config):.12f}")
    print(f"sqrt(Lambda_1)b = {kb_cutoff(config):.12f}")
    print(f"epsilons = {config.epsilon_values}")
    print(
        "Beyn ellipse: "
        f"real endpoints=[{left:.8f}, {right:.8f}], "
        f"center={center:.8f}, rx={rx:.8f}, ry={ry:.3e}"
    )
    print(
        f"Beyn settings: M={config.beyn_M}, "
        f"Nq-levels={config.beyn_quadrature_levels}, "
        f"probe_dim={config.beyn_probe_dim}"
    )
    print(
        f"rank diagnostics: rel_tol={config.beyn_rank_relative_tolerance:.1e}, "
        f"gap_threshold={config.beyn_rank_gap_threshold:.1e}"
    )
    print(f"local verification M = {config.refinement_M}")

    if config.a <= a0_star:
        raise RuntimeError(
            "The selected geometry does not satisfy the leading-order condition a > a0*."
        )

    # Warn if the fixed cutoff margin is so large that the asymptotic branch
    # would be excluded. This is only a safety diagnostic; Beyn itself does not
    # use the asymptotic value for discovery.
    for epsilon in config.epsilon_values:
        kb_asym, _ = asymptotic_prediction(epsilon, config)
        if kb_asym >= right:
            print(
                f"WARNING: epsilon={epsilon:.3f} has kb_asym={kb_asym:.12f} "
                f"to the right of the contour endpoint {right:.12f}. "
                "Reduce beyn_cutoff_margin before interpreting discovery."
            )

    print("\nChecking complex-kb support required by Beyn...")
    complex_support_smoke_test(config.epsilon_values[0], config)
    print("  complex-kb assembly: PASS")

    summaries: list[ValidationResult] = []
    all_beyn_diagnostics: list[BeynDiagnostics] = []
    all_eigen_rows: list[BeynEigenvalueRow] = []
    all_convergence_rows: list[BeynConvergenceRow] = []
    all_refinement_rows: list[ModeRefinementRow] = []
    all_mode_results: list[ModeResult] = []

    for epsilon in config.epsilon_values:
        print(f"\n=== epsilon={epsilon:.3f} ===")
        kb_asym, sigma_asym = asymptotic_prediction(epsilon, config)
        print(f"  predicted kb = {kb_asym:.12f}")
        print(f"  predicted sigma = {sigma_asym:.8e}")
        print(f"  predicted cutoff gap = {kb_cutoff(config) - kb_asym:.3e}")

        (
            result,
            beyn_diags,
            eigen_rows,
            convergence_rows,
            refinement_rows,
            mode_results,
            final_s0_singular_values,
            final_eigen_rows,
        ) = validate_epsilon(epsilon, config)

        summaries.append(result)
        all_beyn_diagnostics.extend(beyn_diags)
        all_eigen_rows.extend(eigen_rows)
        all_convergence_rows.extend(convergence_rows)
        all_refinement_rows.extend(refinement_rows)
        all_mode_results.extend(mode_results)

        print_result(result, config)
        plot_beyn_contour(epsilon, final_eigen_rows, output_directory, config)
        plot_s0_singular_values(epsilon, final_s0_singular_values, output_directory)
        plot_beyn_convergence(epsilon, convergence_rows, output_directory)
        plot_candidate_refinement(epsilon, refinement_rows, output_directory)

    write_dataclass_csv(output_directory / "summary.csv", summaries)
    write_dataclass_csv(output_directory / "beyn_diagnostics.csv", all_beyn_diagnostics)
    write_dataclass_csv(output_directory / "beyn_convergence.csv", all_convergence_rows)
    write_dataclass_csv(output_directory / "beyn_raw_eigenvalues.csv", all_eigen_rows)
    write_dataclass_csv(output_directory / "mode_refinement.csv", all_refinement_rows)
    write_dataclass_csv(output_directory / "mode_results.csv", all_mode_results)
    plot_summary(summaries, output_directory, config)

    one_mode_supported = sum(row.unique_mode_verified for row in summaries)
    asymptotic_passed = sum(row.asymptotic_agreement_verified for row in summaries)

    print("\n=== FINAL SUMMARY ===")
    print("Global discovery: Beyn contour method with Nq convergence study")
    print("Local certification: SVD minima + BEM mesh refinement")
    print("Interval enclosed: 0 < k^2 < Lambda_1, excluding endpoint margins")
    print(
        f"Stable Beyn rank=1 + one locally resolved mode: "
        f"{one_mode_supported}/{len(summaries)}"
    )
    print(
        f"Leading asymptotic sigma within {config.relative_sigma_error_tolerance:.1%}: "
        f"{asymptotic_passed}/{len(summaries)}"
    )

    for row in summaries:
        sigma_error_text = (
            f"{row.relative_error_sigma:.3%}"
            if np.isfinite(row.relative_error_sigma)
            else "--"
        )
        print(
            f"  epsilon={row.epsilon:.3f}: "
            f"Nq={row.beyn_final_quadrature_points}, "
            f"rank={row.beyn_estimated_rank}, "
            f"rank_stable={'yes' if row.beyn_rank_stable else 'no'}, "
            f"near-real={row.beyn_near_real_candidates}, "
            f"resolved={row.resolved_mode_count}, "
            f"one-mode={'PASS' if row.unique_mode_verified else 'FAIL'}, "
            f"sigma_error={sigma_error_text}"
        )

    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
