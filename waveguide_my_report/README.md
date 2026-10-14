# Numerical waveguide article

`main.tex` is the reconstructed Spanish journal-style manuscript; `references.bib` contains the seven retained verified references. `figures/` contains eleven current PDF/PNG displays, `tables/` generated LaTeX, `data/` frozen source snapshots and versioned measured field reconstructions. `evidence_audit.md` bounds claims and records limitations. The manuscript remains subject to author approval and venue-specific preparation.

## Compile

From this directory, with TeX Live/MacTeX and latexmk:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Output: `build/main.pdf`. Citation/label resolution is managed by latexmk/BibTeX. Missing figures fail compilation; there are no placeholder boxes.

## Regenerate current figures and tables

From repository root, using the numerical environment (NumPy,SciPy,pandas,Matplotlib):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python waveguide_my_report/scripts/generate_publication.py
```

The script uses the exact three runs and hashes in `scripts/source_lock.json`, not mutable latest pointers. Source changes fail rather than silently substituting new evidence. It replaces generated publication displays/tables, preserves solver runs, and reuses archived checked field arrays. To create a fresh versioned field reconstruction at the saved roots:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python waveguide_my_report/scripts/generate_publication.py --reconstruct-fields
```

This does not repeat eigenvalue searches or adjust physical roots. Arrays include complex field, phase-fixed boundary vector, coordinates, cross-section norms and modal projections. See data/publication_manifest.json for transformations and per-artifact provenance.

## Audit reproduction

```sh
python docs/audits/manuscript-reconstruction/recompute_acceptance.py
python -m unittest discover -s src/beyn_theorems_validations/tests -p 'test_*.py' -q
```

From repository root. Neither is a full numerical sweep. Detailed independent source/numerical/display/reviewer assessments are in `docs/audits/manuscript-reconstruction/`. Original manuscript is archived there. Legacy `images/` and `prism-uploads/` were moved intact to `docs/audits/manuscript-reconstruction/original-manuscript/`; they are unused and are not current evidence.

No paper submission or external disclosure is performed by these commands. Full-space uniqueness of X/Y and fixed-M continuum convergence remain outside the supported claims. The Example5.4coefficient discrepancy is retained explicitly.

Added kb_vs_epsilon_by_a and kb_vs_epsilon_each_a from the pinned D21 sweep, with figure data in data/kb_vs_epsilon_by_a.csv. Certified finite roots only; exclusions remain absent values.
