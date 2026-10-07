# Numerical validation audit — Theorem 2.1

Scope: circular Neumann obstacle with S=pi and mu=1. The computation is a numerical validation/support study, not a mathematical proof.

- b = 1
- baseline a = 0.6
- leading a0* = 0.391826552031
- first cutoff sqrt(Lambda_1)b = 1.570796326795
- formal small-epsilon reporting window: epsilon <= 0.1

## Claim audit

| Claim | Status | Numerical evidence | Limitation |
|---|---|---|---|
| Existence and uniqueness for sampled supercritical small-epsilon cases | supported | 5/5 sampled cases have stable Beyn rank one and exactly one locally resolved multi-M mode; 5/5 also pass off-grid BIE + decay checks. | Finite numerical sampling supports, but does not prove, theorem-wide uniqueness. |
| Leading asymptotic sigma = C(a) epsilon^2 + O(epsilon^3 |log epsilon|) | supported | On epsilon <= 0.1, relative error of sigma/epsilon^2 from C(a) changes from 0.096% at epsilon=0.02 to 5.557% at epsilon=0.1; Q range [0.0461, 0.901]. | Bounded Q on a finite sample is numerical consistency with the Big-O remainder, not a proof of the asymptotic estimate. |
| Non-existence on sampled subcritical small-epsilon geometries | partial | 1/4 sampled subcritical cases have no resolved interior mode; zero-supported cases use an independent whole-band absolute-SVD scan when Beyn rank is contaminated by the cutoff. | Numerical non-detection has a finite resolution floor and is not a mathematical proof of absence arbitrarily close to the threshold. |
| Critical height a*(epsilon) = a0* + O(epsilon) | partial | 1 epsilon values bracket the zero/one transition; (a_c-a0*)/epsilon lies in [-1.5, -1.5] with maximum final a-bracket width 1.200e-01. | Finite epsilon values demonstrate bounded sampled scaling, not a limit proof. |
| Numerical discretization/implementation convergence | documented | M: max delta_sigma=8.141e-06; finite_difference_step: max delta_sigma=3.880e-08; harmonic_order: max delta_sigma=1.240e-02; lattice_terms: max delta_sigma=3.003e-09 | One-at-a-time convergence study is representative rather than exhaustive over every geometry. |
| Analyticity in epsilon and epsilon log epsilon | theoretical-only | The code reports numerical consistency with the leading expansion but does not test analyticity. | Analyticity is a theorem-level functional property and cannot be established from finitely many numerical samples. |

## Small-epsilon asymptotic data

| epsilon | sigma_BEM | sigma/epsilon^2 | C(a) | relative sigma error | scaled remainder Q |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 1.492325501e-03 | 3.730813753e+00 | 3.734417247e+00 | 0.0965% | 0.0460567 |
| 0.04 | 5.938120860e-03 | 3.711325537e+00 | 3.734417247e+00 | 0.6183% | 0.179346 |
| 0.06 | 1.321751118e-02 | 3.671530882e+00 | 3.734417247e+00 | 1.6840% | 0.372539 |
| 0.08 | 2.310409510e-02 | 3.610014860e+00 | 3.734417247e+00 | 3.3312% | 0.615676 |
| 0.1 | 3.526906992e-02 | 3.526906992e+00 | 3.734417247e+00 | 5.5567% | 0.901206 |

## Critical-height brackets

| epsilon | zero-side a | one-side a | a_c estimate | width | (a_c-a0*)/epsilon | status |
|---:|---:|---:|---:|---:|---:|---|
| 0.04 | 0.271826552 | 0.391826552 | 0.331826552 | 1.200e-01 | -1.5 | partial-bracket-ambiguous |
| 0.06 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |
| 0.08 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |
| 0.1 | nan | 0.511826552 | nan | nan | nan | unbracketed-or-ambiguous |

## Subcritical non-existence search

Across 20 admissible subcritical paper-sweep points: clean zero=0, zero-supported=1, ambiguous=19.
`zero-supported` means the Beyn moment rank was not cleanly zero, but an independent whole-band absolute-SVD search found no resolved interior eigenvalue above the numerical resolution floor.

## Interpretation boundary

Points outside the configured small-epsilon reporting window are retained as finite-size exploration. Turning points or other non-monotone behavior there are not used as evidence for or against the small-obstacle asymptotic statement.
