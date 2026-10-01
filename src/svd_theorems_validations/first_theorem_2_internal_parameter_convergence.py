"""Internal-parameter convergence study for Theorem 2.1.

This compact study keeps the physical problem and the BEM discretization fixed
at

    epsilon = 0.10, M = 32, N_lat = 200,

and varies only the two remaining numerical parameters used by the lattice
Green function and the source finite differences:

    harmonic_order       = 5, 10, 15, 20, 30, 40
    finite_difference_step = 1e-4, 1e-5, 1e-6, 1e-7, 1e-8.

The expected branch is refined in the same sigma-based search used by
``first_theorem_2.py``.  The output records sigma_BEM, the singular-value
diagnostic, and changes relative to a stated reference value.  The singular
value is a diagnostic of a deep minimum; sigma_BEM is the quantity used for
numerical convergence.
"""

from __future__ import annotations

import csv
import io
import math
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt  # pyright: ignore[reportMissingImports]

import first_theorem_2 as theorem


EPSILON = 0.10
M_FIXED = 32
N_LAT_FIXED = 200
HARMONIC_ORDER_VALUES: tuple[int, ...] = (5, 10, 15, 20, 30, 40)
FINITE_DIFFERENCE_STEPS: tuple[float, ...] = (1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8)
BASE_HARMONIC_ORDER = 20
BASE_FINITE_DIFFERENCE_STEP = 1.0e-6
HARMONIC_REFERENCE = 40
FINITE_DIFFERENCE_REFERENCE = 1.0e-6
OUTPUT_DIRECTORY = Path("theorem_2_1_internal_parameter_convergence_eps010")


@dataclass
class InternalParameterRow:
    parameter: str
    value: float
    epsilon: float
    M: int
    N_lat: int
    harmonic_order: int
    finite_difference_step: float
    kb: float
    sigma_bem: float
    sigma_asymptotic: float
    relative_sigma_error_to_asymptotic: float
    sigma_min: float
    drop_factor: float
    minimum_is_interior: bool
    relative_sigma_change_from_reference: float = math.nan
    relative_sigma_change_from_baseline: float = math.nan
    diagnostic_passed: bool = False


def study_config(base: theorem.Config, *, harmonic_order: int, step: float) -> theorem.Config:
    """Keep the physical and BEM settings fixed while changing one parameter."""
    return replace(
        base,
        lattice_terms=N_LAT_FIXED,
        harmonic_order=int(harmonic_order),
        finite_difference_step=float(step),
        refinement_M=(M_FIXED,),
        uniqueness_scan_M=M_FIXED,
        uniqueness_refine_M=M_FIXED,
        parallel_workers=1,
    )


def evaluate(
    base: theorem.Config,
    parameter: str,
    value: float,
    *,
    harmonic_order: int,
    step: float,
) -> InternalParameterRow:
    config = study_config(base, harmonic_order=harmonic_order, step=step)
    # The refinement routine prints progress; keep the CSV-producing study quiet.
    buffer = io.StringIO()
    with redirect_stdout(buffer), redirect_stderr(buffer):
        refinement = theorem.run_expected_branch_refinement(EPSILON, config)
    if len(refinement) != 1 or refinement[0].M != M_FIXED:
        raise RuntimeError(
            f"Expected one refinement row at M={M_FIXED}; "
            f"received {[row.M for row in refinement]}"
        )

    refined = refinement[0]
    _, sigma_asymptotic = theorem.asymptotic_prediction(EPSILON, config)
    relative_error = abs(refined.sigma_bem - sigma_asymptotic) / max(
        abs(sigma_asymptotic), 1.0e-30
    )
    return InternalParameterRow(
        parameter=parameter,
        value=float(value),
        epsilon=EPSILON,
        M=M_FIXED,
        N_lat=N_LAT_FIXED,
        harmonic_order=harmonic_order,
        finite_difference_step=step,
        kb=float(refined.kb),
        sigma_bem=float(refined.sigma_bem),
        sigma_asymptotic=float(sigma_asymptotic),
        relative_sigma_error_to_asymptotic=float(relative_error),
        sigma_min=float(refined.sigma_min),
        drop_factor=float(refined.drop_factor),
        minimum_is_interior=bool(refined.minimum_is_interior),
    )


def relative_change(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-30)


def complete_diagnostics(rows: list[InternalParameterRow]) -> None:
    by_parameter = {
        parameter: [row for row in rows if row.parameter == parameter]
        for parameter in ("harmonic_order", "finite_difference_step")
    }
    references = {
        "harmonic_order": next(
            row for row in by_parameter["harmonic_order"] if row.value == HARMONIC_REFERENCE
        ),
        "finite_difference_step": next(
            row
            for row in by_parameter["finite_difference_step"]
            if math.isclose(row.value, FINITE_DIFFERENCE_REFERENCE, rel_tol=0.0, abs_tol=1.0e-20)
        ),
    }
    baselines = {
        "harmonic_order": next(
            row for row in by_parameter["harmonic_order"] if row.value == BASE_HARMONIC_ORDER
        ),
        "finite_difference_step": next(
            row
            for row in by_parameter["finite_difference_step"]
            if math.isclose(row.value, BASE_FINITE_DIFFERENCE_STEP, rel_tol=0.0, abs_tol=1.0e-20)
        ),
    }
    for parameter, parameter_rows in by_parameter.items():
        reference = references[parameter]
        baseline = baselines[parameter]
        for row in parameter_rows:
            row.relative_sigma_change_from_reference = relative_change(
                row.sigma_bem, reference.sigma_bem
            )
            row.relative_sigma_change_from_baseline = relative_change(
                row.sigma_bem, baseline.sigma_bem
            )
            row.diagnostic_passed = bool(
                row.sigma_min <= 1.0e-4
                and row.drop_factor >= 100.0
                and row.minimum_is_interior
            )


def write_csv(path: Path, rows: list[InternalParameterRow]) -> None:
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def plot_sweep(
    rows: list[InternalParameterRow],
    parameter: str,
    output_directory: Path,
) -> None:
    subset = sorted(
        (row for row in rows if row.parameter == parameter),
        key=lambda row: row.value,
    )
    x = [row.value for row in subset]
    sigma = [row.sigma_bem for row in subset]
    changes = [max(row.relative_sigma_change_from_reference, 1.0e-18) for row in subset]
    if parameter == "harmonic_order":
        xlabel = r"harmonic order $L_h$"
        stem = "harmonic_order"
        reference_label = rf"reference $L_h={HARMONIC_REFERENCE}$"
    else:
        xlabel = r"finite-difference step $h$"
        stem = "finite_difference_step"
        reference_label = rf"reference $h={FINITE_DIFFERENCE_REFERENCE:.0e}$"

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(x, sigma, "o-", color="tab:blue")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(r"$\sigma_{\mathrm{BEM}}$")
    axes[0].set_title(r"Spectral distance")
    axes[0].grid(True, linestyle="--", alpha=0.35)

    axes[1].semilogy(x, changes, "o-", color="tab:orange", label=reference_label)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(r"$|\sigma-\sigma_{\mathrm{ref}}|/|\sigma_{\mathrm{ref}}|$")
    axes[1].set_title(r"Relative change in $\sigma$")
    axes[1].grid(True, which="both", linestyle="--", alpha=0.35)
    axes[1].legend()

    figure.suptitle(
        rf"Internal-parameter convergence: $\varepsilon={EPSILON}$, "
        rf"$M={M_FIXED}$, $N_{{\mathrm{{lat}}}}={N_LAT_FIXED}$"
    )
    figure.tight_layout()
    figure.savefig(output_directory / f"{stem}_convergence.png", dpi=220)
    plt.close(figure)


def main() -> None:
    output_directory = OUTPUT_DIRECTORY
    output_directory.mkdir(parents=True, exist_ok=True)
    theorem.configure_parallel_environment()
    base = theorem.Config()

    rows: list[InternalParameterRow] = []
    for harmonic_order in HARMONIC_ORDER_VALUES:
        rows.append(
            evaluate(
                base,
                "harmonic_order",
                float(harmonic_order),
                harmonic_order=harmonic_order,
                step=BASE_FINITE_DIFFERENCE_STEP,
            )
        )
    for step in FINITE_DIFFERENCE_STEPS:
        rows.append(
            evaluate(
                base,
                "finite_difference_step",
                step,
                harmonic_order=BASE_HARMONIC_ORDER,
                step=step,
            )
        )

    complete_diagnostics(rows)
    rows.sort(key=lambda row: (row.parameter, row.value))
    write_csv(output_directory / "internal_parameter_convergence.csv", rows)
    plot_sweep(rows, "harmonic_order", output_directory)
    plot_sweep(rows, "finite_difference_step", output_directory)

    print("=== Theorem 2.1 internal-parameter convergence study ===")
    print(f"epsilon={EPSILON}, M={M_FIXED}, N_lat={N_LAT_FIXED}")
    print(f"CSV: {output_directory / 'internal_parameter_convergence.csv'}")
    for row in rows:
        print(
            f"{row.parameter:>23s}={row.value:.8g}: "
            f"sigma={row.sigma_bem:.12e}, "
            f"relative change={row.relative_sigma_change_from_reference:.3e}, "
            f"diagnostic={'Sí' if row.diagnostic_passed else 'No'}"
        )


if __name__ == "__main__":
    main()
