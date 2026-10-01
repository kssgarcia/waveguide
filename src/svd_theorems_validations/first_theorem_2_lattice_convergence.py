"""Lattice-sum convergence study for Theorem 2.1.

This script reuses the numerical routines from ``first_theorem_2.py`` and
keeps the BEM discretization fixed at M=16.  It varies the number of lattice
terms

    N_lat = 25, 50, 100, 200

for every epsilon in the standard epsilon sweep.  The expected-mode local
refinement is run independently for each (epsilon, N_lat) pair.  The script
reports the resulting kb, sigma, singular-value diagnostic, and the change
relative to the N_lat=200 result.

This is a discretization/convergence study only.  It intentionally does not
run the paper a-sweep or the whole-band uniqueness scan.  The underlying
``first_theorem_2.py`` local minimizer is left unchanged.
"""

from __future__ import annotations

import csv
import io
import math
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import first_theorem_2 as theorem


M_FIXED = 16
N_LAT_VALUES: tuple[int, ...] = (25, 50, 100, 200)
EPSILON_VALUES: tuple[float, ...] = theorem.Config().epsilon_values
OUTPUT_DIRECTORY = Path("theorem_2_1_lattice_convergence_M16")
PARALLEL_WORKERS = 4
LATTICE_RELATIVE_TOLERANCE = 0.03


@dataclass
class LatticeConvergenceRow:
    epsilon: float
    M: int
    N_lat: int
    kb: float
    cutoff_gap: float
    sigma_bem: float
    sigma_asymptotic: float
    relative_sigma_error_to_asymptotic: float
    sigma_min: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_Nlat200: float = math.nan
    relative_kb_change_from_Nlat200: float = math.nan
    lattice_converged_vs_Nlat200: bool = False


def config_for_lattice_terms(base_config: theorem.Config, n_lat: int) -> theorem.Config:
    """Fix every BEM order used by this study at M=16."""
    return replace(
        base_config,
        lattice_terms=int(n_lat),
        refinement_M=(M_FIXED,),
        uniqueness_scan_M=M_FIXED,
        uniqueness_refine_M=M_FIXED,
        parallel_workers=1,
    )


def evaluate_point(
    epsilon: float,
    n_lat: int,
    base_config: theorem.Config,
) -> tuple[LatticeConvergenceRow, str]:
    """Run the expected-mode refinement for one (epsilon, N_lat) pair."""
    config = config_for_lattice_terms(base_config, n_lat)
    refinement = theorem.run_expected_branch_refinement(epsilon, config)
    if len(refinement) != 1 or refinement[0].M != M_FIXED:
        raise RuntimeError(
            f"Expected exactly one refinement row at M={M_FIXED}, "
            f"got {[row.M for row in refinement]}"
        )

    row = refinement[0]
    _, sigma_asym = theorem.asymptotic_prediction(epsilon, config)
    relative_sigma_error = abs(row.sigma_bem - sigma_asym) / max(
        abs(sigma_asym), 1.0e-30
    )
    result = LatticeConvergenceRow(
        epsilon=float(epsilon),
        M=M_FIXED,
        N_lat=int(n_lat),
        kb=float(row.kb),
        cutoff_gap=float(theorem.kb_cutoff(config) - row.kb),
        sigma_bem=float(row.sigma_bem),
        sigma_asymptotic=float(sigma_asym),
        relative_sigma_error_to_asymptotic=float(relative_sigma_error),
        sigma_min=float(row.sigma_min),
        drop_factor=float(row.drop_factor),
        minimum_is_interior=bool(row.minimum_is_interior),
    )
    return result, ""


def _worker(payload: tuple[float, int, theorem.Config]) -> tuple[LatticeConvergenceRow, str]:
    epsilon, n_lat, config = payload
    # Keep worker logs out of the terminal; the parent prints one concise line
    # per completed point.
    buffer = io.StringIO()
    with redirect_stdout(buffer), redirect_stderr(buffer):
        row, _ = evaluate_point(epsilon, n_lat, config)
    return row, buffer.getvalue()


def complete_reference_comparisons(
    rows: list[LatticeConvergenceRow],
) -> None:
    """Add changes relative to N_lat=200 for each epsilon."""
    reference = {
        row.epsilon: row
        for row in rows
        if row.N_lat == max(N_LAT_VALUES)
    }
    for row in rows:
        ref = reference.get(row.epsilon)
        if ref is None:
            continue
        row.relative_sigma_change_from_Nlat200 = abs(row.sigma_bem - ref.sigma_bem) / max(
            abs(ref.sigma_bem), 1.0e-30
        )
        row.relative_kb_change_from_Nlat200 = abs(row.kb - ref.kb) / max(
            abs(ref.kb), 1.0e-30
        )
        row.lattice_converged_vs_Nlat200 = bool(
            row.relative_sigma_change_from_Nlat200 <= LATTICE_RELATIVE_TOLERANCE
        )


def write_csv(path: Path, rows: list[LatticeConvergenceRow]) -> None:
    if not rows:
        return
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def plot_sigma_vs_epsilon(
    rows: list[LatticeConvergenceRow],
    output_directory: Path,
) -> None:
    plt = theorem.plt
    plt.figure(figsize=(8, 5))
    for n_lat in N_LAT_VALUES:
        subset = sorted((row for row in rows if row.N_lat == n_lat), key=lambda r: r.epsilon)
        plt.plot(
            [row.epsilon for row in subset],
            [row.sigma_bem for row in subset],
            "o-",
            markersize=4,
            label=rf"$N_{{\mathrm{{lat}}}}={n_lat}$, $M={M_FIXED}$",
        )
    asymptotic = sorted(
        {row.epsilon: row.sigma_asymptotic for row in rows}.items()
    )
    plt.plot(
        [epsilon for epsilon, _ in asymptotic],
        [sigma for _, sigma in asymptotic],
        "k--",
        linewidth=1.8,
        label=r"$\sigma_{\mathrm{asym}}$",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$\sigma_{\mathrm{BEM}}$")
    plt.title(rf"Lattice-sum convergence at fixed $M={M_FIXED}$")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "sigma_vs_epsilon_by_Nlat_M16.png", dpi=220)
    plt.close()


def plot_relative_change(
    rows: list[LatticeConvergenceRow],
    output_directory: Path,
) -> None:
    plt = theorem.plt
    plt.figure(figsize=(8, 5))
    for n_lat in N_LAT_VALUES[:-1]:
        subset = sorted(
            (
                row
                for row in rows
                if row.N_lat == n_lat and math.isfinite(row.relative_sigma_change_from_Nlat200)
            ),
            key=lambda r: r.epsilon,
        )
        plt.semilogy(
            [row.epsilon for row in subset],
            [max(row.relative_sigma_change_from_Nlat200, 1.0e-16) for row in subset],
            "o-",
            markersize=4,
            label=rf"$N_{{\mathrm{{lat}}}}={n_lat}$ vs. 200",
        )
    plt.axhline(
        LATTICE_RELATIVE_TOLERANCE,
        color="k",
        linestyle="--",
        linewidth=1.2,
        label=rf"tolerance ({LATTICE_RELATIVE_TOLERANCE:.0%})",
    )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$|\sigma_N-\sigma_{200}|/|\sigma_{200}|$")
    plt.title(rf"Lattice convergence in $\sigma$ at fixed $M={M_FIXED}$")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "relative_sigma_change_vs_Nlat_M16.png", dpi=220)
    plt.close()


def plot_cutoff_gap(
    rows: list[LatticeConvergenceRow],
    output_directory: Path,
) -> None:
    plt = theorem.plt
    plt.figure(figsize=(8, 5))
    for n_lat in N_LAT_VALUES:
        subset = sorted((row for row in rows if row.N_lat == n_lat), key=lambda r: r.epsilon)
        plt.semilogy(
            [row.epsilon for row in subset],
            [max(row.cutoff_gap, 1.0e-16) for row in subset],
            "o-",
            markersize=4,
            label=rf"$N_{{\mathrm{{lat}}}}={n_lat}$",
        )
    plt.xlabel(r"$\varepsilon$")
    plt.ylabel(r"$kb_{\mathrm{cutoff}}-kb_{\mathrm{BEM}}$")
    plt.title(rf"Cutoff gap at fixed $M={M_FIXED}$")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "cutoff_gap_vs_Nlat_M16.png", dpi=220)
    plt.close()


def main() -> None:
    base_config = replace(
        theorem.Config(),
        refinement_M=(M_FIXED,),
        uniqueness_scan_M=M_FIXED,
        uniqueness_refine_M=M_FIXED,
    )
    output_directory = OUTPUT_DIRECTORY
    output_directory.mkdir(parents=True, exist_ok=True)
    theorem.configure_parallel_environment()

    jobs = [
        (float(epsilon), int(n_lat), base_config)
        for n_lat in N_LAT_VALUES
        for epsilon in EPSILON_VALUES
    ]
    workers = max(1, min(PARALLEL_WORKERS, len(jobs)))

    print("=== Theorem 2.1 lattice-sum convergence study ===")
    print(f"fixed BEM order M = {M_FIXED}")
    print(f"N_lat values = {N_LAT_VALUES}")
    print(f"epsilon values = {EPSILON_VALUES}")
    print(f"workers = {workers}")

    rows: list[LatticeConvergenceRow] = []
    if workers == 1:
        for index, (epsilon, n_lat, config) in enumerate(jobs, start=1):
            row, _ = evaluate_point(epsilon, n_lat, config)
            rows.append(row)
            print(
                f"[{index:>3}/{len(jobs)}] N_lat={n_lat:>3}, "
                f"epsilon={epsilon:.3f}, sigma={row.sigma_bem:.8e}"
            )
    else:
        context = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
            futures = {executor.submit(_worker, job): job[:2] for job in jobs}
            for index, future in enumerate(as_completed(futures), start=1):
                epsilon, n_lat = futures[future]
                row, _ = future.result()
                rows.append(row)
                print(
                    f"[{index:>3}/{len(jobs)}] N_lat={n_lat:>3}, "
                    f"epsilon={epsilon:.3f}, sigma={row.sigma_bem:.8e}"
                )

    rows.sort(key=lambda row: (row.epsilon, row.N_lat))
    complete_reference_comparisons(rows)
    write_csv(output_directory / "lattice_convergence_M16.csv", rows)
    plot_sigma_vs_epsilon(rows, output_directory)
    plot_relative_change(rows, output_directory)
    plot_cutoff_gap(rows, output_directory)

    print("\n=== convergence summary relative to N_lat=200 ===")
    for epsilon in EPSILON_VALUES:
        subset = [row for row in rows if math.isclose(row.epsilon, epsilon)]
        print(f"epsilon={epsilon:.3f}")
        for row in subset:
            print(
                f"  N_lat={row.N_lat:>3}: "
                f"relative sigma change={row.relative_sigma_change_from_Nlat200:.3e}, "
                f"within 3%={'YES' if row.lattice_converged_vs_Nlat200 else 'no'}"
            )
    print(f"\nFiles written to: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
