# Numerical validation report — Theorem 2.3(iii)

This run uses the exact x-axis symmetry reduction (odd in y) described in Remark 2.5.

## Geometry and dipole strength

- beta = 0.5
- a = 0
- fixed physical boundary order M = 24; Beyn/SVD use the odd reduced operator
- X-even defect = 0.000e+00
- Y-odd defect = 0.000e+00
- minimum reference speed = 5.000000e-01
- circle calibration mu = 1, nu = 1.560e-14
- final shape mu = 1.07977772066, nu = 1.281e-15
- C = pi^3 mu/b^3 = 33.4798867602

## Small-epsilon branch

| epsilon | odd Beyn rank | resolved modes | kb | sigma/eps^2 | Q(eps) | physical |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.0200 | 1 | 1 | 3.14156425747 | 33.393235 | 1.1075 | PASS |
| 0.0400 | 1 | 1 | 3.14115053022 | 32.940213 | 4.19148 | PASS |
| 0.0600 | 1 | 1 | 3.13947473679 | 32.038242 | 8.54032 | PASS |
| 0.0800 | 1 | 1 | 3.13544403809 | 30.6963 | 13.7762 | PASS |
| 0.1000 | 1 | 1 | 3.12819756571 | 28.980061 | 19.5425 | PASS |

## Claim audit

- **supported** — Theorem 2.3(iii) symmetry assumptions (a=0, X even, Y odd, nu=0): a=0.000e+00, X-even defect=0.000e+00, Y-odd defect=0.000e+00, MFS nu=1.281e-15, mu=1.0797777. Limitation: nu=0 is checked numerically in addition to the exact parametrization symmetry.
- **supported** — One odd embedded trapped mode in [Lambda_1,Lambda_2) for sampled small epsilon: 5/5 small-epsilon cases have stable odd-sector Beyn rank one and exactly one resolved fixed-M SVD mode; 5/5 pass the embedded-field diagnostics. Limitation: Finite sampling supports but does not prove theorem-wide uniqueness.
- **supported** — sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|): At the smallest sampled epsilon, sigma/epsilon^2=33.393235 versus C=pi^3 mu/b^3=33.479887; scaled remainders are finite over the reporting window. Limitation: Bounded finite-sample scaled remainder is numerical evidence, not a proof of Big-O.
- **documented** — Numerical implementation sensitivity at fixed M: Representative one-at-a-time sensitivity at fixed M=24 in finite-difference step, lattice truncation, and harmonic order is saved to internal_convergence.csv. Limitation: Representative convergence is not repeated at every epsilon.
- **theoretical-only** — Analyticity in epsilon and epsilon log epsilon: The paper proves analyticity; finite numerical samples cannot establish analyticity. Limitation: Not a numerically provable statement from a finite sweep.

## Interpretation boundary

The odd-parity reduction is essential: it removes the first even propagating channel, turning the embedded full-space eigenvalue into a discrete eigenvalue of the odd sector. The full reconstructed field is then checked to be odd, to suppress the first open channel, and to decay exponentially.
Points above the configured small-epsilon reporting window are finite-size exploration and are not used to define the formal asymptotic claim.
