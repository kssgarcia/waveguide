# Numerical validation of Theorem 2.3(iv)

## Configuration
- b = 1.0
- beta = 0.5
- fixed BEM order M = 24
- epsilon values = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)

## Independent reference-obstacle quantities
- area S = 3.53429173529
- mu = 1.17693827901
- nu = -2.858e-17
- a1 from theorem formula (2.7) = -0.229718010365
- a1 from equivalent formula (5.47) = -0.229718010365
- relative (2.7)/(5.47) mismatch = 0.000e+00
- Example 5.4 leading value -beta/12 = -0.0416666666667
- circle calibration: mu=1, a1=-3.701e-17

## Methodological note
Unlike Theorem 2.3(iii), no parity reduction removes the continuous spectrum here. The reported BIC is therefore obtained from the full BEM matrix, with the O(epsilon^2) placement correction and sigma tuned numerically at fixed M.

## Claim audit
### Theorem 2.3(iv) symmetry assumptions: X odd, Y even, nu=0
- status: **supported**
- evidence: X-odd defect=0.000e+00, Y-even defect=0.000e+00, nu=-2.858e-17.
- limitation: Numerical symmetry/dipole checks do not replace the analytic symmetry argument.

### Placement law a = epsilon a1 + O(epsilon^2)
- status: **partial**
- evidence: a1 from formula (2.7)=-0.229718010365; scaled correction range [-2.8, 0.0148].
- limitation: Bounded scaled corrections on a finite epsilon sample are numerical consistency with O(epsilon^2), not a proof of the Big-O statement.

### One non-radiating embedded trapped mode for sampled tuned geometries
- status: **partial**
- evidence: 1/5 small-epsilon cases pass full-matrix singularity, x-even field, open-channel suppression, off-grid BIE and decay checks.
- limitation: Finite sampling supports the theorem numerically; it does not prove existence for every sufficiently small epsilon.

### sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)
- status: **supported**
- evidence: relative error of sigma/epsilon^2 from C changes from 0.204% at epsilon=0.02 to 12.895% at epsilon=0.1; scaled-remainder range [0.95,20.4].
- limitation: Finite asymptotic data support the stated scaling but do not prove the remainder bound.

### Uniqueness on the sampled real embedded band
- status: **supported**
- evidence: 1/1 scanned small-epsilon cases found no additional resolved non-radiating real-axis candidate.
- limitation: A finite real-axis SVD screen is numerical evidence, not a theorem-wide proof and not a validated interval-arithmetic exclusion.

### Analyticity in epsilon and epsilon log epsilon
- status: **theoretical-only**
- evidence: This is established analytically in the paper, not by the finite computation.
- limitation: No finite numerical experiment can prove analyticity.

## Per-epsilon summary

| eps | a_num | (a-eps a1)/eps^2 | kb | sigma | rel sv | open/BIC physical |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.020 | -0.0057125497 | -2.79547 | 3.1415588800 | 0.0145672 | 1.425e-08 | FAIL/NA |
| 0.040 | -0.0091817831 | 0.00433579 | 3.1410655166 | 0.0575484 | 1.529e-11 | PASS |
| 0.060 | -0.013729941 | 0.0147611 | 3.1390599772 | 0.126122 | 4.384e-10 | FAIL/NA |
| 0.080 | -0.018619239 | -0.037781 | 3.1342164740 | 0.215155 | 7.742e-09 | FAIL/NA |
| 0.100 | -0.022848614 | 0.0123187 | 3.1254701972 | 0.317869 | 7.008e-09 | FAIL/NA |
| 0.120 | -0.03115597 | -0.249292 | 3.1124633762 | 0.426821 | 4.582e-06 | FAIL/NA |
| 0.140 | -0.035600886 | -0.175529 | 3.0956072723 | 0.535556 | 5.762e-06 | FAIL/NA |
| 0.160 | -0.037871117 | -0.0436029 | 3.0759572280 | 0.63882 | 7.878e-07 | FAIL/NA |
| 0.180 | -0.041465883 | -0.00360004 | 3.0549397902 | 0.732767 | 1.238e-08 | FAIL/NA |
| 0.200 | -0.045898044 | 0.00113894 | 3.0339590472 | 0.815289 | 1.921e-10 | FAIL/NA |