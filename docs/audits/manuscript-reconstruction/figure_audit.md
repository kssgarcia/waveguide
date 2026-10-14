# Figure and table audit before manuscript reconstruction

Audit date: 2026-10-13. Applied repository-local `scientific-visualization/SKILL.md` and its publication-integrity guidance. No manuscript text, manuscript figures, or numerical raw outputs were modified by this audit. The destination journal is unspecified: the design below is a provisional general journal layout, not a claim of publisher compliance.

## Current draft provenance

`waveguide_my_report/main.tex` contains six table environments and sixteen figure environments (some contain two panels). It embeds `images/geo_problema.jpeg`, `images/bem_potential.png`, obstacle images, and thirteen result plots. Twelve of the fifteen PNG files in `images/results` are **byte-identical** to corresponding files under `src/svd_theorems_validations/images/results`: every current spectral/radiation PNG except `t21_kb_vs_epsilon_and_a.png`. The two obstacle PNGs and height-sweep PNG differ; difference does not establish current scientific provenance. All existing numerical plots should be regenerated from current Beyn data, rather than retained on the basis of appearance. Geometry may be redrawn directly from current parametrizations. The reference field `bem_potential.png` needs an identifiable generating experiment before reuse.

Draft lines 1026–1050 label M=32/40 data as mesh convergence. Current preferred datasets are fixed M=24, so these rows and conclusions cannot serve as current convergence evidence. Draft Y captions at lines 1370 and 1388 use `a=-epsilon beta/12`; this is superseded by direct equation-(2.7) evaluation and numerical tuning. Old propagating-channel fraction figures cannot substitute for current complex open-channel projections and physical certificates.

## Frozen source datasets

Use explicit run IDs, not a mutable `latest_run.json`, in publication scripts and provenance manifests:

- D21: `src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-13_17-52-48_full_5605ecd253a8`.
- DX: `src/beyn_theorems_validations/results/theorem_2_3_iii/runs/2026-10-13_21-23-19_full_fe445a6af342`.
- DY: `src/beyn_theorems_validations/results/theorem_2_3_iv/runs/2026-10-13_20-16-45_full_63388d81f58d`.

All three pointers currently identify these runs. The shared audit `docs/audits/latest-complete-runs-2026-10-13/verification.json` provides prior completeness and physical-check context. Recompute manuscript acceptance from raw columns and numerical tolerances; PASS labels alone are insufficient.

| Data | Scientific use | Constraints |
|---|---|---|
| All `summary.csv` | sigma, sigma/epsilon², asymptotic coefficient, normalized/equilibrated residual, existence | 10 epsilon values .02–.20; fixed M=24 |
| D21/DX `beyn_convergence.csv`, `beyn_diagnostics.csv` | rank agreement, S0/S1 changes, candidate drift and quadrature refinement | Missing first-level changes mean not applicable, not zero |
| D21 `paper_a_sweep_M24.csv` | 99-point placement/epsilon diagram, 79 one-mode and 20 supported exclusions | Check geometry validity and diagnostic completeness; excluded finite contours do not prove full spectrum absence |
| D21 `critical_height.csv`, `critical_height_probes.csv` | Four brackets around critical height | Plot actual brackets, not confidence intervals; bracket width .0009375 |
| All `physical_diagnostics.csv` | wall, BIE, parity, open-channel and decay comparisons | Scalar checks are not field arrays |
| DX `mu_mfs_convergence.csv`, `mu_circle_calibration.csv` | auxiliary Laplace/MFS sensitivity and circle benchmark | Internal residuals and independent identities have different meanings |
| DY `mfs_convergence.csv`, `a1_robustness.csv` | direct mu, nu, a1 and equation-(5.47) consistency | Preserve source-placement and quadrature variation |
| DY `placement_a1_fit.csv`, `tuned_modes.csv` | actual placement, epsilon*a1 prediction, scaled correction, root status | Fit selection/point count must accompany extrapolation; finite fit not proof |
| DY `a1_small_beta_diagnostic.csv`, `a1_small_beta_fit.csv` | unresolved printed Example 5.4 discrepancy | No substitution of -beta/12 as computed a1 |
| All `internal_convergence.csv`; DY `internal_physical_certificates.json` | H, lattice and finite-difference sensitivity | One-parameter sensitivity, not multi-M convergence; include failed variants |
| DY `additional_candidates.csv`, screen diagnostics | explain finite screening scope and unresolved candidates | Screen has no rigorous full-space uniqueness count |

## Inspection of current representative plots

Viewed a contact sheet made directly from seven latest-run PNGs: all three `sigma_over_epsilon_squared.png`, D21 `beyn_s0_convergence_epsilon_0.02.png`, DX/DY harmonic-order plots, and DY `placement_a1_fit.png`. Numerical curves and theoretical lines are consistent with the corresponding variables, but publication redesign is needed:

1. D21 plot entitled “Beyn moment convergence” shows singular values of S0 only. Singular-value stabilization does not establish convergence of both moment matrices. Replace or supplement with actual S0/S1 relative changes and candidate drift.
2. DX harmonic-order plot labelled “Numerical convergence” is an internal sensitivity plot at fixed M. Caption it accordingly. Retained zeros on symlog axes are exact numerical equality to the reference, not demonstrated arbitrarily small error.
3. DY correctly exposes rejected H=5 using a cross while stable higher orders pass. Preserve this distinction in regenerated figures; do not filter rows before plotting.
4. Leading-coefficient horizontal lines and finite-epsilon departures are physically informative. Limit claims to agreement approaching the small-epsilon end; no large-epsilon validity claim from smooth curves.
5. DY placement fit displays a narrow vertical range with a dotted direct a1 reference. This is legitimate for position encoding but must state the fit range and model; a linear extrapolation is not independent proof of an asymptotic expansion.

## Proposed publication panels

Use consistent b-scaled coordinates, mathematical labels and marker/line redundancy; PDF vector line plots plus PNG previews. Suggested main-paper set:

1. **Geometry and spectral setting**, three panels: circular obstacle at a, X-symmetric profile, Y-symmetric profile at tuned a. Draw walls y/b=±1 and axes; inset or aligned spectral axis marks Lambda1/Lambda2 and distinguishes full/odd-y spectrum. Source is actual geometry functions/config, not copied article artwork. Equal aspect ratio; display epsilon=.10 and disclose magnification if a reference profile is enlarged.
2. **Scaled eigenvalue comparison**, three panels using each summary: sigma/epsilon² versus epsilon, horizontal direct leading coefficient. One color per dataset with marker and distinct dashed asymptotic line. Either separate y scales with labels or divide by the coefficient and label the exact normalization. Show all ten accepted baseline points. Companion table reports representative .02/.10/.20 errors; scaled remainders in supplement cannot assert an asymptotic error exponent from ten points.
3. **Beyn certification**, D21 and DX columns at representative .02/.10/.20: (a) threshold/gap rank by quadrature, (b) S0 and S1 relative changes, (c) candidate drift. Source raw diagnostics/convergence CSV. Do not put DY in an eigenvalue-count panel. Mention X counts apply to odd-y restriction and finite contours. Show tighter-contour check numerically in table or caption.
4. **Critical height and exclusion**, D21 panels: placement/epsilon scatter with explicit one-mode, zero-supported, ambiguous categories; bracket plot a versus epsilon with exact lower/upper endpoints and leading a0* reference. Overlay no invented smooth boundary; bracket bars are deterministic numerical enclosures, not stochastic error bars. NaN sigma in exclusions stays absent rather than becoming sigma=0.
5. **BIC placement**, DY panels: a_BIC/epsilon versus epsilon with direct a1, and (a_BIC-epsilon*a1)/epsilon². The first panel may show finite-window fitted line as a separate dotted extrapolation and disclose fit selection. Second panel exposes departure of finite tuning from leading formula. Supplement shows small-beta disagreement with -beta/12 separately; never blend contradictory predictions into a fitted “theory” curve.
6. **Physical diagnostics**, X/Y panels: normalized spectral residual, off-grid BIE residual, first-channel relative amplitude, and decay relative discrepancy versus epsilon. Threshold lines come from actual configuration and legends explain distinct diagnostics. Do not compare these quantities as if one universal tolerance applied. Show wall/parity diagnostics in a compact table or supplement; tiny inherited odd parity is not an independent discovery.
7. **Representative full fields and decay**, one D21 and one DY example (.10) with DX parity example if space permits. Real phase-fixed fields use a diverging zero-centered colormap; absolute magnitude uses a sequential map. Mask obstacle interior and a disclosed near-boundary exclusion band; never interpolate over the obstacle. Align near-field spatial axes. Far-field decay panel should use scaled x sigma/b or explicitly extended physical x, since small sigma gives a long decay length. Plot measured cross-section L2 norms for both signs, normalized at the first sampled section, with exp(-sigma*distance/b) comparison. Include first-channel coefficient/projection panel for BIC where available.

## Field/decay source gap and feasible reconstruction

The latest three run directories contain **zero `.npz` and zero `.npy` files**; diagnostics CSVs contain fitted rates and relative errors, not transverse samples, field grids, boundary vectors, or norm arrays. Plotting exp(-sigma*x) from scalar sigma alone would be a theoretical curve, not measured decay evidence.

Saved roots permit a targeted postprocessing experiment without a full sweep. Reassemble the current operator at saved epsilon, a, kb, recover a singular vector, and reconstruct actual fields with the existing functions. In D21 use `assemble_matrix`, `boundary_nodes`, `make_extended_field_green`, `reconstruct_field_at_points`; in DX use `odd_singular_pair` for the production diagnostic and preferably independently assemble full matrix for publication parity check; in DY use `full_singular_pair`, `modal_open_channel_coefficients`, `reconstruct_field_at_points`. Save roots, configuration, vector phase/normalization, coordinates, complex field, quadrature weights, transverse samples, L2 norms and open-channel projections to a new versioned publication-data directory. Compute and compare scalar diagnostics against existing rows before use. Do not tune the eigenvalue during plotting.

Near-obstacle lattice reconstruction loops scale roughly with field points × M × harmonic expansion size; use a modest grid (~120×60) first and exclude points close to the boundary where off-surface quadrature is unreliable. Far-separated sections can use the existing modal representation efficiently. Check modal/lattice overlap and source interpolation sensitivity; do not claim that a denser image increases the underlying M=24 resolution. State vectors are determined only up to complex phase: choose phase by one nonzero reference component, retain the complex field, report any imaginary leakage. Source-oversampled diagnostics do not increase BEM M.

## Tables

Replace six old tables with compact generated tables, preferably 3 representative epsilon rows per theorem in the main text plus full ten-row data in supplement. Include epsilon, kb, sigma, sigma/epsilon², relative asymptotic error, normalized spectral residual, off-grid BIE residual and physical status. X column for count explicitly says odd-y count; DY table includes a_leading, a_BIC, scaled placement correction, open-channel amplitude. A separate configuration/sensitivity table records M=24, H, lattice terms, finite-difference step and tolerances. Include all 15 sensitivity variants and mark DY's three rejected variants (H=5, h=1e-5, h=1e-4); seven reference-quality rows pass physical certificates. Neither absence of extra screen candidates nor one discovered root creates a full uniqueness column labelled YES.

## Publication integrity requirements

Source manifest must include frozen run ID, CSV SHA256, generating script and exact transformations for every panel/table. Keep raw outputs intact. Record deterministic brackets and sensitivity differences honestly; no statistical CI exists in these calculations. Log zeros require visible zero/equality categories or documented symlog handling, never arbitrary floors interpreted as measured values. Never connect missing/rejected rows as accepted data. Use common diagnostic scales when comparing cases, clear units, accessible marker/line semantics, and final-size review. Target-journal dimensions, font and export requirements remain to be checked once a journal is selected. Export and inspect PDF metadata and compiled manuscript placement; automated DPI checks alone are not scientific acceptance.

## Post-audit implementation and verification

After the supervising audit gate, `waveguide_my_report/scripts/generate_publication.py` generated nine current-data figure PDFs and opaque 300dpi PNG previews, five compact tables plus numeric macros, CSV source copies, measured section data, and `data/publication_manifest.json`. The invocation and exact artifact names are documented in `scripts/GENERATION.md`; direct dependencies are pinned in `scripts/requirements.txt`. No raw solver output or legacy image was replaced, no eigenvalue was retuned, and no full sweep was launched.

The generator freezes 47 scientific source files via SHA256 and checks source changes before processing. An in-process deliberately mismatched digest was rejected; all 24 figure/table output hashes matched after the final default replay. All three baseline datasets retain ten rows; all sensitivity variants remain in the exported source tables. Final visual inspection covered all nine figures. It exposed and corrected two presentation defects during development: initial rank-axis bounds clipped first-level threshold rank 3, and a1 robustness initially used collocation order even though order 160 is fixed while boundary quadrature changes. The delivered rank plot shows rank 3 against gap rank 1 at the first D21 quadrature level, both actual moment-matrix changes and the correct 1e-3 tolerance. No visual axis or category hides an unsuccessful certification.

Versioned X/Y field reconstruction at epsilon=.10 uses unrestricted full M=24 matrices and saved roots. Its archive is `data/field_reconstructions/2026-10-14_03-09-49` (UTC timestamp; local date remains October13). X normalized full-matrix residual is 1.4866e-11, observed odd boundary parity defect 1.0510e-10 and production-section channel ratio 4.1176e-12. Y residual is 2.7570e-12, observed even boundary parity defect 4.9620e-11 and production-section channel ratio 5.1898e-12. Seven-point batched-versus-original scalar reconstruction discrepancies are 3.0559e-11 and 1.5353e-11 respectively. Measured production-section decay rates match prior diagnostic rates to roughly 1e-12. New extended measured sections, complex fields and modal projections are preserved; the exponential reference curve is not substituted for measured norms.

`figure_export_metadata.json` records delivered PNG dimensions/DPI/opacity and PDF page sizes. PDFs are single-page, 7.1 in wide, and contain embedded CID TrueType text; field maps are intentionally rasterized inside vector PDFs. PNG alpha is uniformly 255. The final field plot uses a zero-centered diverging scale, disclosed phase/amplitude normalization and a saved near-obstacle mask. Exact target-journal requirements remain unverified because no target journal was specified. Final placement inside the compiled article still requires the supervising manuscript QA.
