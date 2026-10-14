# Independent theory review of reconstructed manuscript

Read the complete current `waveguide_my_report/main.tex` and `references.bib` on 2026-10-14, checked the source preprint and previous primary-source bibliography verification. No manuscript changes made. The manuscript was being edited concurrently, so refer to equation labels rather than treating snapshot line numbers as immutable.

## Required notation corrections

1. **Define exterior Laplace coordinates and radius in `eq:laplace`.** Currently `X,Y` denote bounded contour functions `X(t),Y(t)` immediately before the far-field expression, while `r` is undefined. A far-field expansion cannot literally be in the bounded boundary coordinates. Introduce an independent exterior point `(zeta,eta)` with `r=sqrt(zeta^2+eta^2)`, write `Psi(zeta,eta)=-mu eta/r^2-nu zeta/r^2+O(r^-2)` as `r->infty`, and use `Psi(X(t),Y(t))` only for its boundary trace in `eq:dipoles` and `eq:a1`. The circle calibration must use the same exterior coordinates. Also specify the domain of `Delta Psi=0` as the plane exterior to Gamma. This is a real mathematical ambiguity, not merely prose style.

2. **Define the Green difference variables before `eq:lattice`.** The periodic function introduces `X,Y,r,vartheta` without defining their relation to its argument; those capitals already designate the obstacle contour. Explicitly define local separations `(Delta x,Delta y)` and `r=sqrt(Delta x^2+Delta y^2)`, with `vartheta=atan2(Delta y,Delta x)`. The reflected evaluation uses its own reduced vertical difference. Renaming these variables avoids accidentally identifying the Green distance with the reference obstacle radius. Specify `P=(x,y)`, `q=(xi,eta)` and the diagonal limit if those are not introduced elsewhere.

3. **Make source-version attribution explicit at first theorem presentation.** The bibliography cites the journal version, but the equations/Example5.4 audit was performed against local March4,2025 preprint only. The final Discussion correctly discloses no line-by-line published comparison. Add one early sentence identifying the local preprint as the version used for equation/theorem numbering and the discrepancy assessment. Do not imply all published equations were verified unchanged. The body already avoids calling the discrepancy a proved error of the published article.

## Checks passed

- Problem PDE and mixed conditions are correct. `\\partial_nu` in TeX parses as `\\partial` with subscript `n`, followed by `u`; it is **not** an undefined control sequence or a derivative of Greek nu. Writing `\\partial_n u` is clearer but mathematically identical.
- Normal `(-Y',X')/speed` points into an anticlockwise obstacle and outward from its exterior domain. The circle `mu=1,nu=0,S=pi` calibration matches the signed exterior problem.
- Thresholds, modal normalization `integral phi_n^2=b`, and modal coefficient `-1/(2b gamma_n)` are consistent with `(Delta+k^2)G=+delta`. Open-channel branch `gamma=-i sqrt(k^2-Lambda)` gives outgoing `exp(+i k_ax |x-xi|)` for time factor `exp(-i omega t)`.
- Period `4b`, reflected separation `y+eta+2b`, and periodic wrap formula are correct. Free-space `-iH0/4`, regularization removing `Y0/4`, double-layer representation/jump, and negative curvature term agree. The matrix factor `I-(4pi/M)K^w` follows from `u/2=(2pi/M)K^w u`.
- Theorem2.1 coefficient, critical height and logarithmic remainder match local source (2.4)-(2.5). Main numerical coefficient and circle critical height are consistent. Critical coefficient caveat and positive-displacement limitation are stated.
- Theorem2.3 coefficient and integral a1 match local (2.6)-(2.7). Shifting integration from `[-pi,pi]` to `[0,2pi]` is valid for the stated periodic contour.
- Odd-y channel removal is correct: phi1 even/phi2 odd. X sector-counting is clearly separated from global full-space uniqueness. The Y geometry preserves x parity and cannot eliminate all first-channel radiation; no global Beyn count is claimed.
- Numerical root/refinement, eigenvalue-real-part tuning with subsequent complex acceptance, physical open-channel cancellation and decay checks are distinguished. No absence is inferred from failed/ambiguous refinement.
- The text correctly observes that a right singular vector's normalized residual equals `smin/smax` in exact arithmetic; it does not falsely call those independent scientific evidence. Row-equilibrated checks add an explicit scale-control criterion.
- MFS signed log fundamental solution and source dipole moments agree with the implementation. Zero total sources eliminates log r; using the displayed kernel without an additive constant provides zero far-field constant. Source scale0.65, order sequence and quadrature description match the stated method.
- Source Example5.4 discrepancy is correctly restricted to the supplied preprint. Its printed f and (5.47) give numerator `-3pi beta/2`, denominator `3pi+O(beta)`, yielding `-beta/2`; direct a1 uses (2.7). No invented correction is asserted.
- FixedM24 limitations, incomplete coverage, finite asymptotic window, lack of uniform Big-O/analyticity proof and dependence of built-in wall/decay/parity diagnostics are acknowledged.

## Bibliography consistency

Citation-key scan found **zero undefined and zero unused entries**. Seven retained entries correspond to verified sources. Beyn2012 volume436(10),pages3839-3863,DOI10.1016/j.laa.2011.03.030 are correct. Linton1998 DOI10.1023/A:1004377501747 and current publisher pages377-401 are correct, with local print377-402 discrepancy retained in the audit. Sargent2019 issue year is correct despite2018 online date. Zhevandrov2025 `pages={8}` denotes article8 rather than a page span; an article-number annotation or `number={1},note={Article8}` is optional depending on journal style. Kassis2026 authors/title/arxiv/DOI are verified against arxiv2609.00065; it is appropriately cited only for procedural assistance, not numerical evidence.

No bibliography existence or DOI mismatch found among retained entries. This check does not imply that every scientific sentence citing a contextual paper was independently reproduced from that paper's full text. Existing local Linton1992 content supports integral-method lineage; journal metadata verifies existence only for other contextual literature.

## Remaining review scope

This review validates theorem representation, signs and bibliography, not all CSV-derived numbers or regenerated figure contents. The parent numerical and visual audits remain necessary. After the coordinate definitions and early source-version attribution are clarified, no material theory/sign correction is identified from this review. The paper remains evidence at a fixed finite matrix dimension, not a continuous-operator proof or validated numerical error bound.
