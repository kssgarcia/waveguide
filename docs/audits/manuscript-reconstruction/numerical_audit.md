# Independent numerical evidence audit for the reconstructed article

Audited 2026-10-13. This audit preserves the existing computed outputs. No production sweep or numerical implementation was changed. The code was read directly after the supervising agent's CodeGraph availability check failed; the current CSV values were checked independently of their acceptance labels.

## Data admitted

All three datasets use full boundary order M=24, b=1 and ten epsilon values 0.02,0.04,...,0.20. The X reduced operator has dimension 12; M remains the full boundary dimension, not a physical parameter.

| Case | Exact run under `src/beyn_theorems_validations/results/<case>/runs/` | Scope |
|---|---|---|
| theorem_2_1 | 2026-10-13_17-52-48_full_5605ecd253a8 | 10 baseline heights a=.6; 99 geometry points; four critical-height brackets; 15 sensitivity variants |
| theorem_2_3_iii | 2026-10-13_21-23-19_full_fe445a6af342 | 10 x-axis-symmetric shapes beta=.5; 15 sensitivity variants; odd-y discovery |
| theorem_2_3_iv | 2026-10-13_20-16-45_full_63388d81f58d | 10 y-axis-symmetric tuned placements beta=.5; 15 sensitivity variants; full embedded-mode certification |

Current `latest_run.json` pointers agree with these directories. Earlier `embedded_validation/` folders, SVD folders and pre-fix outputs must not be mixed into principal results. The Oct13 completed-run provenance audit records implementation fingerprints; `recomputed_acceptance.json` additionally records fresh SHA256 values for each admitted summary, physical and sensitivity table.

## Recomputed acceptance

`recompute_acceptance.py`, run from the repository root with the project Python, independently asserts dispersion identities, scaled sigma, asymptotic relative error, finite residuals, thresholds, spectral drop, wall/BIE/decay diagnostics, X full-vector leakage, Y root residual and signed-placement identity, and the intentional limits on uniqueness. All 30 principal cases satisfy the relevant quantitative thresholds. This is an audit of fixed discretizations, not a proof of the continuous eigenproblem.

| Maximum over ten principal cases | 2.1 | 2.3(iii) X | 2.3(iv) Y |
|---|---:|---:|---:|
| normalized spectral residual | 6.708e-10 | 9.294e-11 | 1.293e-11 |
| row-equilibrated vector residual | 6.880e-10 | 9.567e-11 | 1.319e-11 |
| wall relative residual | 7.653e-10 | 2.452e-11 | 4.753e-11 |
| off-grid BIE relative residual | 7.112e-6 | 4.033e-3 | 2.117e-3 |
| relative measured-decay error | 4.379e-7 | 2.897e-8 | 3.574e-8 |
| first open-channel relative amplitude | closed | 1.303e-10 | 3.537e-11 |

The BIE tolerance is 5e-3, substantially weaker than the spectral and wall tests. In particular the smallest-epsilon X BIE value 4.033e-3 is close to this tolerance. A paper must state this quantitative bound rather than claiming uniform machine-precision physical verification. Zero wall residual is also partly built into the Dirichlet Green representation, so it is an implementation consistency check. Near-identical decay rates use the same spectral root and Green formulation; they are not independent operator convergence proofs.

Spectral certification thresholds are normalized residual <=1e-8, row-equilibrated residual <=1e-8, relative singular value <=1e-4, spectral drop >=100 and an interior minimum. The row-scaled minimum singular value thresholds are 1e-6 (2.1/X), 1e-8 (Y). Y also uses smin normalized by a reference row scale <=1e-4. These are implementation acceptance criteria, not rigorously established error bounds.

The normalized residual of an exact singular vector equals smin/smax. It is explicitly recorded but cannot be counted as statistically or mathematically independent evidence. Row equilibration and physical checks are the additional protections against misleading scale separation (`first_theorem_beyn.py:633–672`, `second_theorem_X.py:1490–1513`, `second_theorem_Y.py:1290–1307`).

## Discovery, counting and scope

2.1 uses the full BEM matrix, contour moments, local SVD refinement and field diagnostics. All ten baseline roots have one-mode contour evidence. Maximum final relative changes of S0 and S1 are 5.0092e-4, below 1e-3, with stable rank, stable candidate, controlled solves and tighter-contour agreement. Code locators: `first_theorem_beyn.py:687–699` and geometry classification at `3867–3924`. The statement allowed is one numerical mode in the inspected finite contour at the specified discretization, not theorem-wide uniqueness.

The 99 geometry rows contain 79 `one` and 20 `zero-supported`. Each absence row has complete diagnostics, zero agreement between threshold/gap/selected ranks and controlled linear solve residual. All twenty absence rows have complete two-level full-determinant exclusion files (listed in the recomputed JSON), with zero winding at 260 and516 contour points, controlled phase increments, singularity checks and solved contour discretizations. Seventeen of these twenty also contain suspicious unresolved local SVD minima; it is the independently resolved exclusion, not their failed local refinement, that supports absence in those seventeen cases. This is finite-resolution argument-principle evidence for the discretized matrix inside a recorded ellipse. It stops at kb=1e-4 and a cutoff gap of 1e-12 where the independent exclusion is used. Endpoint intervals and continuous-problem non-existence are not established.

**Rank-zero nuance:** relative moment changes can be O(1) because moments are near zero. Code uses an absolute-smallness criterion for all four moment snapshots when successive ranks are both zero (`first_theorem_beyn.py:708–714`), rather than pretending those relative changes are below 1e-3. Describe this as absolute rank-zero moment control, and preserve the independent exclusion distinction. NaN kb/sigma in absence rows means no eigenvalue is reported, not that absence was inferred from NaNs.

Critical heights are four numerical brackets (eps=.04,.06,.08,.10) of width 9.375e-4 in `critical_height.csv`, centered near a0*=0.3918265520306073. The midpoint is an estimate with half-bracket uncertainty 4.6875e-4. Curvature of the scaled midpoint shift should not be interpreted as a resolved higher-order coefficient.

X uses Q* A Q with an odd-y lifting matrix Q. Final threshold/gap/selected ranks agree at one in all ten cases, with final S0/S1 changes <=7.1231e-4. Uniqueness is supported **only in the odd sector inside the finite contour**. The lower threshold margin is 1e-4 and the upper margin is case dependent. The even sector and omitted margins are not excluded (`second_theorem_X.py:1930–1961`, flags at `2074–2076`). Theorem-wide/full-space uniqueness must not be claimed.

Y solves an embedded-mode tuning problem using the full matrix, not a justified global Beyn count. Its finite real-band screen has no additional resolved BICs but deliberately sets `uniqueness_supported=False`. Omitted intervals are present at both thresholds (`second_theorem_Y.py:2248` and each summary row's `screen_diagnostics_json`). No full-space uniqueness claim is admissible.

## X parity: imposed versus independently observed

The standard X physical reconstruction uses the vector lifted from the odd sector (`second_theorem_X.py:1718`). Odd parity and open-channel suppression in that field are partly imposed; they cannot independently prove that a full-space eigenvector has the same parity.

To address this issue, this audit recomputed a full M=24 SVD **without imposing parity** at eps=.02,.10,.20. The smallest full singular vectors have ||Rv+v||/||v|| <=4.049e-10, overlap with the odd-discovered vectors >=0.9999999999999989, normalized full residual <=4.811e-11 and singular-value separation >=9.166e9. `X_independent_full_vectors.json` records these new results. This independently supports compatibility of the odd discovery roots with isolated singular vectors of the full discretized operator. It does not create a full-spectrum eigenvalue count or a separate independently reconstructed physical field.

## Asymptotic claims and placement

The asymptotic formulas are inputs for initial guesses and comparison curves; the reported root values are refined numerical quantities. Relative sigma errors at epsilon=.02 are 0.096623% (2.1), 0.258817% (X), 0.201065% (Y). At epsilon=.20 they are 23.5282%,45.0419%,44.1467%, respectively. Under a descriptive 5% accuracy cutoff, sampled agreement holds through .08 (2.1), .06 (X/Y). The full sampled range must not be described as quantitatively agreeing with a leading-order asymptotic formula.

The X/Y MFS quantities in the principal shape are mu=1.0797777206648984 / 1.1769382790137215 and nu=1.2807e-15 / -2.8575e-17. X-axis reflection and Y-axis reflection both give numerically vanishing nu for these selected shapes. This does not test Theorem2.3(ii)'s contrasting nonzero-nu case.

Y directly computes a1=-0.22971801036478778 using formula(2.7) (`second_theorem_Y.py:760–786`); the BIC root tuning uses a=epsilon*a1+c*epsilon^2. The fitted intercept from the five physical BICs eps<=.10 is -0.22967826186578144, relative difference 0.017303%. The quadratic-fit intercept has relative difference 0.037224%; these are finite-window fits, not exact limits. `placement_a1_fit.csv` is the source. Leading spectral coefficient is 36.49247391473535, versus33.47988676015933 for X and3.7344172474882895 for baseline2.1.

The Example5.4 discrepancy remains substantive. For beta=.02,.05,.10 the direct (2.7)/(5.47) values agree internally and have a fitted slope -0.4986263012, whereas the printed -beta/12 coefficient has slope -0.0833333333. These formula cross-checks share the same Psi and are not independent PDE solves (`second_theorem_Y.py:783–789`). Preserve the discrepancy and direct a1 provenance. Do not identify -beta/12 as the main a1 or silently replace it with a favorable root-fitted coefficient.

## Sensitivity and independent checks

All studies keep M=24. Fifteen variants examine derivative step, harmonic order and lattice terms one at a time near epsilon=.10; this is parameter sensitivity, not a boundary mesh convergence study. Y retains twelve certified and three ambiguous variants; h=1e-5,h=1e-4 and harmonic_order=5 must remain visible as rejected low-accuracy settings. Reference settings pass; ambiguity must not be portrayed as a physical contradiction or hidden.

Previously stored independent modal checks use 12000/24000 terms for representative roots, and Y also varies field quadrature/source interpolation and verifies placement starting at a=0. These checks remain available in `independent_modal_operator.json` (X), `independent_validation.json`, `independent_placement.json` and `internal_physical_certificates.json` (Y). They support the Green/operator implementation at sampled points; they do not establish all parameter regimes. Exact periodic reduction is implemented for ordinary and regularized assembly evaluations; current outputs carry the corrected assembly versions.

## Claims not tested or still limited

Theorem2.3(i) absence away from Lambda2 and2.3(ii) contrasting nonzero-nu geometries are not established by these three ten-point datasets. No analyticity theorem, exact Big-O assertion, continuum eigenvalue proof, full-space embedded-mode uniqueness, endpoint exclusion or boundary-mesh convergence follows from these outputs. Principal evidence can support existence/localization and strong sampled physical non-radiation diagnostics; the theorem provides the analytical statements.

## Norm definitions and audit execution

Spectral backward residuals use the matrix spectral norm ||A||2=smax and Euclidean vector norm. Row equilibration divides each nonzero row by its Euclidean norm before computing the additional residual. Green evaluation and the physical operator are unchanged by this diagnostic scaling.

The off-grid quantity is max |(1/2)u_interpolated - K u| divided by the largest of the target half-density, integral term and boundary-density peak (plus a numerical floor). It uses shifted targets and trigonometric density interpolation (`first_theorem_beyn.py:2680–2740`; X:1671; Y:1866). This checks off-grid consistency of the boundary integral equation. It is **not a direct hypersingular evaluation of the obstacle normal derivative**, not a continuum error estimate, and not an independently proved Neumann boundary error bound.

Wall/centerline/parity values divide pointwise sampled field residuals by max(boundary-density peak, sampled cross-section peak). First-channel amplitudes use a transverse projection normalized by the first-mode norm and transverse field norm. Decay uses a linear fit of log transverse L2 norms versus physical longitudinal distance. Such normalized quantities are scale invariant but use different norms and sample sets; tables must not conflate them.

Executed audit experiments: full-space X SVD at three representative existing roots; raw-data recomputation for all thirty principal cases and twenty geometry absences; `python -m unittest discover -s src/beyn_theorems_validations/tests -q`, 90 tests passed in1.950s. Pytest is not installed in the project Python; the repository tests are unittest tests and were run using that supported runner. No expensive sweep was rerun. No manuscript source or numerical solver code was edited by this audit.
