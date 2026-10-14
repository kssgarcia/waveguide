# Scientific evidence audit

This is an author-directed, AI-assisted working manuscript. Computational/source checks below are not human sign-off, an editorial decision, or a mathematical proof. Target journal and human approvals are pending. No raw solver output was overwritten.

## Frozen datasets and provenance

| ID | Current dataset | Purpose |
|---|---|---|
| D21 | `../src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-13_17-52-48_full_5605ecd253a8` | 10 baseline modes,99 placements,4 critical brackets,15 sensitivities |
| DX | `../src/beyn_theorems_validations/results/theorem_2_3_iii/runs/2026-10-13_21-23-19_full_fe445a6af342` | 10 odd-sector roots and full field diagnostics,15 sensitivities |
| DY | `../src/beyn_theorems_validations/results/theorem_2_3_iv/runs/2026-10-13_20-16-45_full_63388d81f58d` | 10 tuned full-space BICs,MFS/a1,15 sensitivities |

The generator checks pinned hashes in `scripts/source_lock.json`, copies admitted CSVs into `data/source/`, and records each figure/table transformation in `data/publication_manifest.json`. Numerical source code, configurations and raw runs remain outside the manuscript folder. `docs/references/WG_03.03.25.pdf` is the directly audited theoretical preprint, not assumed byte-equivalent to the final journal article.

## Claim–evidence matrix

Status means computational support at the stated scope. `Supported` never means a continuous-operator proof. D21/DX/DY expand to the paths above; auxiliary audit filenames below are in `../docs/audits/manuscript-reconstruction/`.

| ID | Claim/source | Numerical method | Exact evidence | Status | Limitations / allowed statement |
|---|---|---|---|---|---|
| T21-E | existence side | full BEM, Beyn discovery, sigma SVD, physical checks | D21 summary.csv; mode_results.csv; physical_diagnostics.csv | supported | 10 baseline cases, fixed M24; continuum error not bounded |
| T21-N | non-existence side | rank zero plus independently resolved determinant winding | D21 paper_a_sweep_M24.csv; geometry_diagnostics/*/contour_exclusion.json | partially supported | 20 finite-ellipse exclusions; lower kb1e-4 and upper gap1e-12; no threshold-wide proof |
| T21-U | uniqueness | converged S0/S1, rank agreement, candidate and tighter contour | D21 summary.csv; beyn_diagnostics.csv | partially supported | uniqueness inside finite contours for discretized matrix; not theorem-wide |
| T21-C | critical-height behavior | classification and bisection | D21 critical_height.csv; critical_height_probes.csv | partially supported | four brackets width9.375e-4; no resolved higher coefficient |
| T21-A | eigenvalue asymptotics | independent refined sigma vs leading formula | D21 summary.csv | supported | finite epsilon .02-.20; 5percent agreement only .02-.08; no BigO proof |
| T21-P | physical trapped-mode diagnostics | walls, offgrid BIE, measured transverse decay | D21 physical_diagnostics.csv | supported | integral consistency not independent hypersingular Neumann derivative |
| T23-I | absence away from Lambda2 | finite exploratory screening only | DY screen diagnostics | not tested | does not test quantified absence away from Lambda2 |
| T23-II | role of nonzero nu | symmetric geometries have nu approximately zero | DX/DY MFS convergence.csv | not tested | no nonsymmetric nu nonzero experiment |
| X-S | x-axis symmetry and odd-y sector | geometry and reflection operator | DX shape_diagnostics.csv; source reference_obstacle_geometry | supported | beta .5 centered a0; physical parity mostly imposed in production |
| X-D | discrete spectrum after odd reduction | theoretical transverse odd-y restriction | WG Remark2.5 pp4-5; DX odd_lift_matrix | supported | analytical sector argument; no full-space count |
| X-B | Beyn rank/counting and moments | graded contour quadrature and moment/rank gates | DX beyn_diagnostics.csv; beyn_convergence.csv | supported | finite analytic matrix contours only |
| X-R | SVD refinement | sigma minimization and original/equilibrated residual | DX mode_results.csv | supported | fixed discretization |
| X-P | full field and open channel | full-operator check plus reconstructed field | DX physical_diagnostics.csv; X_independent_full_vectors.json | supported | independent unconstrained parity checked .02/.10/.20; BIE max.00403 |
| X-A | sigma asymptotics | mu MFS and sigma leading formula | DX mu_mfs_convergence.csv; summary.csv | supported | 5percent agreement .02-.06; error45percent at .20 |
| X-U | uniqueness scope | odd-sector converged Beyn count | DX summary.csv | partially supported | odd-sector contour uniqueness only; full-space unresolved |
| Y-S | y-axis symmetry and truly embedded spectrum | shape symmetry and theoretical channel spectrum | DY shape_diagnostics.csv; WG Remark2.5 | supported | x parity retains n1 open channel |
| Y-M | numerical a1 from equation2.7 | exterior Laplace MFS and direct quadrature2.7 | DY mfs_convergence.csv; a1_robustness.csv | supported | beta.5 a1=-.229718010365; no use of -beta/12 |
| Y-L | placement law | joint root tuning a=eps*a1+c*eps2; finite extrapolation | DY summary.csv; placement_a1_fit.csv; independent_placement.json | partially supported | finite epsilon consistency; asymptotic remainder not proven |
| Y-A | eigenvalue asymptotics | mu and independently tuned sigma | DY summary.csv | supported | 5percent agreement .02-.06; error44percent at.20 |
| Y-R | full BEM singularity | full complex eigenvalue, signed radiation root, SVD residuals | DY tuned_modes.csv; summary.csv | supported | track branch and condition Jacobian; fixed M24 |
| Y-O | complex open-channel cancellation | left/right complex outgoing modal coefficients and projection | DY physical_diagnostics.csv; independent_validation.json | supported | finite amplitude tolerances; not exact zero |
| Y-P | physical BIC | walls, BIE, parity, radiation, measured decay | DY physical_diagnostics.csv; independent_validation.json | supported | BIE max.00212; no continuum error bound |
| Y-U | uniqueness scope | finite sampled scan only | DY summary.csv; additional_candidates.csv | unresolved | global uniqueness never established |
| Y-E | minus beta/12 discrepancy | printed formula5.47 algebra vs printed example coefficient | WG p24; theory_bibliography.md; DY a1_small_beta_diagnostic.csv | contradicted | which source expression should be corrected remains unresolved |

## Definitions and quantitative decisions

M24 is a boundary discretization choice. b1,epsilon .02:.02:.20, lattice200, harmonic20, finite-difference1e-6; X/Y shape_beta.5. All30 principal spectral/physical decisions were recomputed from raw columns by `recompute_acceptance.py` before rewriting. Both S0 and S1 must converge for count support; rank-zero uses absolute-smallness control. All20 absent geometry cases have resolved full-determinant exclusions (260/516effective points);17 of these also have suspicious unresolved SVD minima. Failure, NaN and lack of a discovered root never supply absence evidence.

r_A=||Av||2/(||A||2||v||2), with spectral matrix norm. B scales each row of A by its Euclidean norm; r_B uses the SAME v, and smin(B) is a separate quantity. For an exact singular vector r_A equals smin/smax, so these are not independent evidence. Absolute thresholds operate on the controlled scale. Off-grid BIE is max|u_interp/2-Ku| normalized by boundary/target/integral amplitudes; it is not a directly evaluated hypersingular Neumann derivative or a continuum error bound. Wall/parity tests and modal decay partly inherit structure of the representation.

The X production field uses an odd lift. Fresh unconstrained full matrix SVD at .02/.10/.20 independently observed odd vectors; saved measurements in `X_independent_full_vectors.json`. Publication field reconstruction uses unrestricted full matrices at existing roots, retaining complex values and no root retuning. These checks support compatibility, not full-spectrum uniqueness.

At5percent descriptive asymptotic accuracy,2.1 sampled agreement extends through.08; X/Y through.06. This cutoff is not used to classify existence. Y direct a1=-.22971801036478778; finite-window linear intercept -.22967826186578144 differs.017303percent. Exact BigO and analyticity are theoretical. Y sensitivity h1e-5,h1e-4,H5 remains ambiguous/rejected.

## Figures and tables

Eleven newly generated scientific figures replace all legacy numerical figures. Geometry uses current parameterizations; spectral scaling uses all30roots; Beyn panel uses actual moment changes/rank/candidate drifts; critical figure uses actual99classifications and deterministic brackets; placement uses direct a1 and tuned roots; physical diagnostics use scalar measurements; sensitivity retains failed variants; MFS figure preserves the printed-example disagreement; field/decay uses archived measured complex arrays at savedroots. See per-artifact source lists, hashes, normalization and transformations in `data/publication_manifest.json`, and field reconstruction manifests in `data/field_reconstructions/`.

Generated tables list all10root values percase plus representative physical diagnostics; source CSVs supply fullprecision. Rounded display digits are not uncertainty estimates. Auxiliary configuration and numeric macros are reproducible. No statistical confidence intervals are claimed. Critical brackets are deterministic numerical classifications, not rigorous continuum enclosures.

## Superseded evidence and changes

The previous manuscript depended on legacy SVD outputs, M32/M40 mesh claims and Y placement -epsilon beta/12. Twelve old result images matched legacy images byte-for-byte. These displays/tables/claims are replaced, while old figures remain archived outside the manuscript folder and unused. The original manuscript/bibliography is preserved outside the paper in `../docs/audits/manuscript-reconstruction/original-manuscript/`.

No numerical solver code was changed in this manuscript task. New work:90unittest regressions,30raw-acceptance recomputations,20exclusion rechecks,3unconstrained X singularvectors,2saved-root field/decay reconstructions and publication regeneration. No expensive full sweep was rerun.

## Citation/source checks

Core references: source theory; Linton1992integral-method lineage; Linton1998Greenevaluation; Evans1994symmetryexistence context; Sargent/Mestel2019integralwaveguide context; Beyn2012contourmatrixmethod. Public publisher/Crossref metadata were checked in `bibliography_metadata.json`; actual scientific source support is recorded in `theory_bibliography.md`. Peripheral unverified/redundant references were removed. Linton1998publisher pages377-401 differ from local masthead/endpage377-402; publisher metadata used without claiming the mismatch resolved. Skills citation Kassis2026 verified against arXiv primaryrecord, is procedural assistance only.

## Unresolved issues and submission status

Printed Example5.4 gives-beta/12 while its displayed ingredients in5.47 give-beta/2; direct2.7is used. It is unresolved which source expression should be corrected and whether the final published version fixes it. FixedM24 has no boundary-mesh error estimate. Full-space embedded uniqueness, endpoint exclusions, theorem2.3(i)/(ii), uniform asymptotic remainders and analyticity are not established by the datasets. These are explicit scope limits, not concealed successful claims.

The reconstructed manuscript can be circulated for scientific review. It is not marked ready for submission: authors must verify content/source locators and author list, approve AI disclosure, funding/conflict/contribution/data-access declarations, select the journal and check its template/policies. No declarations or approvals have been invented. A journal requiring continuum error/convergence estimates would require additional numerical analysis; current article restricts claims accordingly. The preprint/example discrepancy should be checked with the final article or source authors before a definitive statement about the published text.

Added kb_vs_epsilon_by_a and kb_vs_epsilon_each_a from the pinned D21 sweep, with figure data in data/kb_vs_epsilon_by_a.csv. Certified finite roots only; exclusions remain absent values.

The kb-height figures now include dashed leading-order curves from the saved C1(a), with sigma=C1(a)*epsilon² and kb=sqrt((pi/2)²-sigma²), only for positive C1 and real kb. The threshold is dotted. Dense theoretical samples are exported separately in data/kb_vs_epsilon_theory.csv; they are predictions, not additional computed eigenvalues.
