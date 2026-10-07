# Corrected numerical validation — Theorem 2.1
Generation: corrected; assembly: theorem-2.1-periodic-green-v8.
Legacy/pre-fix evidence (not reused): /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/theorem_2_1_publication_ready.
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
| Uniqueness for sampled supercritical small-epsilon cases | partial | 0/5 sampled cases have stable Beyn rank one and converged S0/S1 moments, stable candidate, tighter-contour check and one certified fixed-M mode; 5/5 also pass off-grid BIE + decay checks. | Finite numerical sampling supports, but does not prove, theorem-wide uniqueness. |
| Local mode existence for sampled supercritical small-epsilon cases | supported | 5/5 have explicit singular-vector and equilibrated certificates plus configured physical diagnostics. | Local existence support does not imply a converged global count. |
| Leading asymptotic sigma = C(a) epsilon^2 + O(epsilon^3 |log epsilon|) | numerical-consistency | On epsilon <= 0.1, relative error of sigma/epsilon^2 from C(a) changes from 0.097% at epsilon=0.02 to 5.560% at epsilon=0.1; Q range [0.0461, 0.902]. | Bounded Q on a finite sample is numerical consistency with the Big-O remainder, not a proof of the asymptotic estimate. |
| Non-existence on sampled subcritical small-epsilon geometries | partial | 0/4 sampled subcritical cases have no resolved interior mode; zero-supported cases use an independent whole-band absolute-SVD scan when Beyn rank is contaminated by the cutoff. | Numerical non-detection has a finite resolution floor and is not a mathematical proof of absence arbitrarily close to the threshold. |
| Critical height a*(epsilon) = a0* + O(epsilon) | not-established | No numerical zero/one transition was successfully bracketed. | Critical-height study must produce finite brackets. |
| Numerical implementation sensitivity at fixed M | documented | finite_difference_step at fixed M=24: max delta_sigma=9.826e-08; harmonic_order at fixed M=24: max delta_sigma=5.111e-09; lattice_terms at fixed M=24: max delta_sigma=4.315e-09 | Fixed-M one-at-a-time sensitivity is representative, not exhaustive over every geometry; no mesh-convergence claim is made. |
| Analyticity in epsilon and epsilon log epsilon | theoretical-only | The code reports numerical consistency with the leading expansion but does not test analyticity. | Analyticity is a theorem-level functional property and cannot be established from finitely many numerical samples. |

## Small-epsilon asymptotic data

| epsilon | sigma_BEM | sigma/epsilon^2 | C(a) | relative sigma error | scaled remainder Q |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 1.492323556e-03 | 3.730808889e+00 | 3.734417247e+00 | 0.0966% | 0.0461188 |
| 0.04 | 5.938091873e-03 | 3.711307421e+00 | 3.734417247e+00 | 0.6188% | 0.179487 |
| 0.06 | 1.321736501e-02 | 3.671490281e+00 | 3.734417247e+00 | 1.6851% | 0.37278 |
| 0.08 | 2.310363831e-02 | 3.609943485e+00 | 3.734417247e+00 | 3.3332% | 0.616029 |
| 0.1 | 3.526797438e-02 | 3.526797438e+00 | 3.734417247e+00 | 5.5596% | 0.901681 |

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
