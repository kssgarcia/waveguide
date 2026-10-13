# Corrected numerical validation — Theorem 2.1
Generation: corrected; assembly: theorem-2.1-periodic-green-v11.
Legacy/pre-fix evidence (not reused): /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_1/legacy/pre_fix.
Existence is local certification plus configured physical diagnostics; uniqueness additionally requires converged S0/S1 moments, agreement of threshold/gap rank, stable rank/candidate and tighter contour.

Scope: circular Neumann obstacle with S=pi and mu=1. The computation is a numerical validation/support study, not a mathematical proof.

- b = 1
- baseline a = 0.6
- fixed boundary order: M = 24 for Beyn discovery, SVD certification, and physical diagnostics
- leading a0* = 0.391826552031
- first cutoff sqrt(Lambda_1)b = 1.570796326795
- formal small-epsilon reporting window: epsilon <= 0.1

## Claim audit

| Claim | Status | Numerical evidence | Limitation |
|---|---|---|---|
| Existence/uniqueness and leading asymptotic | not-run | No finite baseline results are available. | Run the baseline validation. |
| Non-existence on sampled subcritical small-epsilon geometries | partial | 0/4 sampled subcritical cases have no resolved interior mode; zero-supported cases use an independent whole-band absolute-SVD scan when Beyn rank is contaminated by the cutoff. | Numerical non-detection has a finite resolution floor and is not a mathematical proof of absence arbitrarily close to the threshold. |
| Critical height a*(epsilon) = a0* + O(epsilon) | not-established | No numerical zero/one transition was successfully bracketed. | Critical-height study must produce finite brackets. |
| Numerical implementation sensitivity at fixed M | not-run | Internal convergence study was not run. | Fixed-M one-at-a-time sensitivity is representative, not exhaustive over every geometry; no mesh-convergence claim is made. |
| Analyticity in epsilon and epsilon log epsilon | theoretical-only | The code reports numerical consistency with the leading expansion but does not test analyticity. | Analyticity is a theorem-level functional property and cannot be established from finitely many numerical samples. |

## Critical-height brackets
Upper bounds require local certified existence; they do not assert uniqueness. The legacy column name a_upper_one_mode denotes this existence bound.

| epsilon | zero-side a | one-side a | a_c estimate | width | (a_c-a0*)/epsilon | status |
|---:|---:|---:|---:|---:|---:|---|
| 0.04 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |
| 0.06 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |
| 0.08 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |
| 0.1 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |

## Subcritical non-existence search

Across 20 admissible subcritical paper-sweep points: clean zero=0, zero-supported=0, ambiguous=20.
`zero-supported` means the Beyn moment rank was not cleanly zero, but an independent whole-band absolute-SVD search found no resolved interior eigenvalue above the numerical resolution floor.

## Interpretation boundary

Points outside the configured small-epsilon reporting window are retained as finite-size exploration. Turning points or other non-monotone behavior there are not used as evidence for or against the small-obstacle asymptotic statement.
