# Author-directed adversarial review — working draft

This is an internal scientific assessment authorized by the manuscript owner's pasted reconstruction instructions, not an assigned journal referee report or an editorial decision. It uses the authorized local manuscript, repository theory/code and computed outputs; no manuscript text was sent to a citation search or external review service. No target journal has been specified, so journal-specific confidentiality, submission, declarations and AI-assistance requirements remain for the accountable authors. No independent human specialist review is asserted. Scope is numerical/mathematical claim alignment; bibliography and rendered-figure audits are being handled separately.

## Materials and reproduction limits

Read the rewritten `waveguide_my_report/main.tex` in full and compared the principal quantitative statements with current summary/physical CSVs, the numerical audit and direct code. Main plots and table generation were still under construction when this assessment was made, so their final renderings are outside this review. Prior audit recomputed thirty principal acceptances, twenty geometry exclusions, three unconstrained full-X singular vectors, and ran90 regressions; this review did not launch a new sweep.

## Required corrections

### R1 — Row-equilibrated residual refers to the wrong vector

**Location:** `main.tex`, discretization/certification paragraph containing “mediante el residual de su vector singular” (line158 at this review).

**Observation:** This describes the residual as computed using a singular vector of B. In all three certification functions, v is the smallest right singular vector of A (of A_odd for X); the same v is evaluated under B. The singular values of B are computed separately, without extracting its right singular vector.

**Evidence:** `first_theorem_beyn.py:648–661`; `second_theorem_X.py:1493–1500`; `second_theorem_Y.py:1293–1298`.

**Why it matters:** With the actual diagnostic, r_B=||Bv||2/(||B||2||v||2) need not equal smin(B)/smax(B). Replacing v with a singular vector of B changes the diagnostic and falsely implies its redundancy with B's singular-value ratio.

**Requested action:** State explicitly that B tests the same A singular vector and define r_B accordingly. Retain the independent absolute threshold on smin(B). No code or numerical rerun is needed.

### R2 — Make the signed Y root coupling reproducible

**Location:** `main.tex`, Y joint tuning paragraph following equation(tuning).

**Observation:** The prose correctly separates real crossing/coupling from final complex certification, but omits that the signed source coupling uses Re(v) after phase alignment and removes the outgoing Green prefactor. “Una proyección firmada” is insufficient to recover the two-component root equation from the text.

**Evidence:** `second_theorem_Y.py:1347–1363` evaluates outgoing coefficients using `vec.real.astype(complex)` and forms

`F1 = Re(lambda)/scale`,
`F2 = Re{i 2b q1 [C_left(Re v)+C_right(Re v)]/2}/epsilon^2`,

where q1=sqrt(k^2-Lambda1), and the vector/eigenvalue follows an overlap-selected branch. Final acceptance evaluates the full complex matrix and complex outgoing coefficients. The projection is supported by the x-even real-density branch; it is not by itself a guarantee that a general complex eigenvector does not radiate.

**Why it matters:** The root residual and physical complex cancellation are distinct quantities. A reproducible method must identify which density and normalization enter the root solve.

**Requested action:** Add the two residual equations or an equivalent compact precise description including Re(v), branch phase alignment and the subsequent complex checks. No additional experiment is needed for the stated bounded evidence.

## Clarifications that would improve appraisal

### R3 — Gauge statement for MFS

**Location:** MFS paragraph saying the constraint sum c_j=0 removes the logarithmic term “y fija el gauge de la expansión”.

**Observation/evidence:** Zero total source removes the logarithmic far-field term; the chosen logarithmic-source representation without an added constant selects the vanishing far-field constant. An additive harmonic constant is not excluded by the zero-charge condition alone.

**Requested action:** Say that the zero-charge constraint removes the logarithmic term and the no-constant representation enforces the stated far-field gauge.

### R4 — Operational tolerances need a bounded interpretation

**Location:** Acceptance methods and principal results.

**Observation:** The manuscript appropriately labels tolerances as operational and acknowledges no continuum a posteriori estimate. X's worst off-grid BIE residual is4.033e-3 compared with acceptance5e-3. These checks are evidence of consistency at fixed M, not high-accuracy Neumann error bounds.

**Requested action:** Preserve the existing explicit qualifications and quantitative maxima during polishing; do not replace them with a generic robust/high-precision claim. A continuous-error claim would require new mathematics/mesh evidence and is unnecessary for the present bounded article.

## Claim checks that pass

- Thirty principal modes and the reported maximum residuals agree with raw output, including the comparatively large off-grid integral residuals.
- Ten-point epsilon range, b=1,M=24, lattice200,h=1e-6,H=20 and X dimension12 match saved configurations.
- Twenty geometry absences all have complete zero full-determinant counts at260/516points; seventeen is the subset with suspicious unresolved SVD minima. This distinction is now correct.
- Four critical brackets and deterministic width9.375e-4 are stated with appropriate uncertainty interpretation.
- Sigma asymptotic errors and the sampled 5%windows are accurate. The text does not treat rejected asymptotic agreement as absence of a mode.
- X's partial circularity from odd-lifted field diagnostics is acknowledged; the additional full-A SVD results accurately establish independently observed odd singular vectors at three roots without claiming a full count.
- Full-space embedded-mode uniqueness is not claimed. X's count is explicitly restricted to the odd sector and finite contours; Y's finite screen is exploratory.
- Theorem2.3(i)/(ii) are separated from the experiments; the symmetric configurations do not test nonzero-nu exclusion.
- The normalized A residual is correctly acknowledged as equal to smin/smax for an exact singular vector in the same spectral norm.
- Off-grid BIE consistency is distinguished from direct hypersingular obstacle Neumann verification and continuum error bounds.
- Direct Y a1, finite-window intercept fits and the Example5.4 discrepancy are not conflated.
- Fixed M, finite cutoff exclusions, lack of exact Big-O/analyticity proof, shared decay reconstruction and built-in walls are stated clearly.

## Completion conditions for the supervising author

Correct R1 and provide the precise Y coupling as R2, then verify generated tables and all plotted source transformations against the pinned runs. Complete citation/source-version checks, particularly that the discrepancy is attributed to the supplied preprint rather than asserted of an unexamined published revision. Compile and inspect the final PDF. Human authors must review the mathematical signs, claimed physical interpretation and target-journal requirements. This assessment cannot constitute submission approval.
