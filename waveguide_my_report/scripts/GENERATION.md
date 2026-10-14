# Publication artifacts

From the repository root (Python with NumPy, pandas, SciPy, Matplotlib):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python waveguide_my_report/scripts/generate_publication.py
```

This checks `source_lock.json` against SHA256 hashes, copies the frozen latest baseline CSVs and run configurations into `data/source/`, regenerates eleven figure PDFs/300-dpi PNG previews, five tables plus numeric macros, and writes `data/publication_manifest.json`. It reuses the hash-checked completed field reconstruction. Publication outputs are replaced by this command; solver raw outputs and legacy images are preserved. PDF metadata may vary by creation timestamp even when scientific data do not.

To create a **new versioned** field reconstruction at the same saved roots:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python waveguide_my_report/scripts/generate_publication.py --reconstruct-fields
```

This assembles unrestricted full M=24 matrices for X and Y at epsilon=.10, takes the smallest right singular vector, and reconstructs complex fields and transverse sections. It does not search for or retune eigenvalues. Complex grids, masks, vectors, sections, norms, weights and channel coefficients are saved in `data/field_reconstructions/<UTC timestamp>/`. The `latest.json` pointer is updated only after both examples pass spectral, scalar/batched field and decay checks. The generator also exports `data/field_sections.csv`, containing measured section norms and complex first-channel projections. A failed reconstruction preserves its in-progress folder and leaves the previous pointer intact.

The near-field kernel calls the existing Y batched Green implementation with explicit X or Y source coordinates and normals; seven points per case are compared with the current case-specific scalar reconstruction. Far-separated sections use the same double-layer modal formula as the solver. The full X matrix is used without imposing parity on its singular vector. Source geometry and b are unchanged, and M remains 24. The .025b near-boundary exclusion band protects the field image from unresolved near-surface quadrature; it is saved as a mask, not filled or interpolated. The field image is phase-fixed using a boundary-vector component and normalized by its own maximum amplitude. Decay curves use actual sampled cross-section norms, with each side normalized at its first section. The exponential reference is shown separately.

`source_lock.json` freezes scientific provenance, not merely paths. If a source differs, the command stops; updating that lock requires a new scientific audit, rather than silently adopting new results. Legacy images are preserved outside the manuscript folder under `docs/audits/manuscript-reconstruction/original-manuscript/images/`, unused by the reconstructed manuscript. Exact target-journal size/font/export requirements are pending journal selection. Default figures are 7.1 inches wide, use vector text/lines and rasterized field maps, and have matched PNG previews; these choices alone do not establish publisher compliance.

Generated figure filenames:

- `geometry`: actual physical obstacle scale at epsilon=.10, three cases;
- `spectral_scaling`: normalized leading coefficient comparison and signed scaled remainders;
- `beyn_convergence`: rank agreement, both moment changes and candidate drift, D21/X only;
- `critical_height`: all 99 geometry classifications and four deterministic brackets;
- `placement`: all ten tuned Y placements and scaled correction;
- `physical_diagnostics`: normalized spectral, BIE, open-channel and decay checks;
- `sensitivity`: all 15 parameter variants per case, with rejected Y variants explicit;
- `mfs_coefficients`: auxiliary mu sensitivity, direct a1 robustness and printed-example discrepancy;
- `fields_decay`: saved-root unrestricted X/Y fields and measured decay.

Tables `theorem_2_1.tex`, `theorem_X.tex`, `theorem_Y.tex` contain all ten baseline rows. `diagnostics.tex` contains .02/.10/.20 physical diagnostics in all cases; `configuration.tex` records implementation parameters. `numeric_macros.tex` contains data-derived coefficients, representative errors and diagnostic maxima; full raw columns remain in companion CSVs. The X rank refers only to the odd-y discovery contour, not full-space uniqueness. Y receives no Beyn eigenvalue-count claim.

Added kb_vs_epsilon_by_a and kb_vs_epsilon_each_a from the pinned D21 sweep, with figure data in data/kb_vs_epsilon_by_a.csv. Certified finite roots only; exclusions remain absent values.
