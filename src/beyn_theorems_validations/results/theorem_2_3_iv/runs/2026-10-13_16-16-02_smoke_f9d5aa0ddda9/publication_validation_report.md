# Numerical validation of Theorem 2.3(iv)

## Configuration
- assembly version = `2026-10-13-Y-signed-channel-certified-v3`
- Result provenance: corrected/v3. Historical embedded_validation and v2 artifacts are legacy/pre-fix.
- Existence and full-band uniqueness are separate. The finite screen never certifies full-band uniqueness.
- Formula (2.7) is retained; the printed Example 5.4 coefficient disagrees with its own (5.47) substitution.
- b = 1.0
- beta = 0.5
- fixed BEM order M = 24
- epsilon values = (0.02, 0.1, 0.2)
- Green assembly uses exact vertical periodic reduction before evaluating the local lattice expansion.

## Independent reference-obstacle quantities
- area S = 3.53429173529
- mu = 1.17693827901
- nu = -2.858e-17
- a1 from theorem formula (2.7) = -0.229718010365
- a1 from algebraically equivalent formula (5.47) = -0.229718010365
- relative (2.7)/(5.47) mismatch = 0.000e+00
- a1 MFS source-scale/quadrature worst relative variation = 3.625e-16
- Example 5.4 leading value -beta/12 = -0.0416666666667
- circle calibration: mu=1, a1=-3.701e-17

The value from (5.47) is an algebraic cross-check using the same computed Psi. The independent numerical check of a1 is the epsilon->0 placement-branch intercept below.

## Independent a1 estimate from the BIC placement branch
Insufficient resolved small-epsilon placements for an independent a1 intercept fit.

## Small-beta Example 5.4 diagnostic
- fitted slope a1/beta = -0.498626301208; paper slope -1/12 = -0.0833333333333; relative difference = 498.352%.
- This discrepancy is reported as unresolved; the script does not replace formula (2.7) by -beta/12.

## Methodological note
Unlike Theorem 2.3(iii), no parity reduction removes the continuous spectrum here. The full BEM matrix is used. Formula (2.7) supplies only the leading placement coefficient: a=epsilon*a1+c*epsilon^2, with c determined numerically. The final acceptance requires both relative and row-scaled singularity criteria plus physical non-radiation diagnostics.

## Claim audit
### Theorem 2.3(iv) symmetry assumptions: X odd, Y even, nu=0
- status: **supported**
- evidence: X-odd defect=0.000e+00, Y-even defect=0.000e+00, nu=-2.858e-17.
- limitation: Numerical symmetry/dipole checks do not replace the analytic symmetry argument.

### Direct numerical evaluation of a1 from theorem formula (2.7)
- status: **supported**
- evidence: a1=-0.229718010365 from the exterior-Laplace MFS; max source-scale/quadrature relative variation=3.625e-16.
- limitation: The (5.47) value is an algebraic cross-check using the same Psi, not a fully independent PDE discretization. The BIC-branch intercept below is the independent check.

### Finite-sample consistency with placement law a = epsilon a1 + O(epsilon^2)
- status: **partial**
- evidence: a1 from formula (2.7)=-0.229718010365; scaled correction range [-0.00454, -0.00274]; insufficient verified points for an independent epsilon->0 intercept.
- limitation: The full-BEM placement intercept is independent of the MFS a1 calculation, but finite epsilon sampling still cannot prove the O(epsilon^2) remainder.

### One non-radiating embedded trapped mode for sampled tuned geometries
- status: **partial**
- evidence: 2/2 small-epsilon cases pass full-matrix singularity (relative and row-scaled smin), x-even field, signed open-channel suppression/phase consistency, off-grid BIE and decay checks.
- limitation: Finite sampling supports the theorem numerically; it does not prove existence for every sufficiently small epsilon.

### Finite-sample consistency with sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)
- status: **partial**
- evidence: Using physically verified BICs only, relative error of sigma/epsilon^2 from C changes from 0.201% at epsilon=0.02 to 12.895% at epsilon=0.1; scaled-remainder range [0.938,20.4].
- limitation: Finite asymptotic data support the stated scaling but do not prove the remainder bound.

### No additional resolved non-radiating candidate on the sampled real band
- status: **supported**
- evidence: 2/2 scanned small-epsilon cases found no additional candidate passing the strengthened singularity and open-channel filters.
- limitation: This is a finite real-axis screen and must not be interpreted as a proof of theorem-wide uniqueness.

### Uniqueness of the embedded mode on the full open band
- status: **ambiguous/not-established**
- evidence: A finite real-axis screen is not a certified spectral count.
- limitation: Unscreened cutoff intervals and finite-grid exclusion are not covered by a global certificate.

### Example 2.4 / 5.4 small-beta coefficient a1 ~ -beta/12
- status: **unresolved-inconsistency**
- evidence: Direct (2.7) computations give fitted a1/beta slope -0.4986263 versus paper -1/12=-0.083333333.
- limitation: Substituting the printed f0, f1 into (5.47) gives -beta/2, not -beta/12; the MFS agrees with that substitution. The displayed example is internally inconsistent; (2.7) remains the numerical definition, without forcing a printed reference value.

### Analyticity in epsilon and epsilon log epsilon
- status: **theoretical-only**
- evidence: This is established analytically in the paper, not by the finite computation.
- limitation: No finite numerical experiment can prove analyticity.

## Per-epsilon summary

| eps | a_num | (a-eps a1)/eps^2 | kb | sigma | rel sv | row-scaled smin | physical |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.020 | -0.0045954577 | -0.00274379 | 3.1415588782 | 0.0145676 | 9.545e-12 | 1.899e-11 | PASS |
| 0.100 | -0.023017221 | -0.00454196 | 3.1254703840 | 0.317867 | 3.052e-12 | 6.267e-12 | PASS |
| 0.200 | -0.04589081 | 0.00131981 | 3.0339589480 | 0.81529 | 4.018e-12 | 8.134e-12 | PASS |