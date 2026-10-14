# Theory and bibliography audit

Audited 2026-10-14. This document audits the supplied March 4, 2025 preprint `docs/references/WG_03.03.25.pdf`, the local Linton PDFs, the current manuscript, and publisher-deposited Crossref metadata. It does not change the manuscript. Public citation metadata are saved in `bibliography_metadata.json`; that file is provenance, not a substitute for checking an article's scientific content.

## Exact problem and conventions

The source problem (1.1)-(1.2) is `-Delta u=k^2 u` on the infinite strip `-b<y<b` outside the physical obstacle `(epsilon X(t),a+epsilon Y(t))`, with Dirichlet walls and homogeneous Neumann obstacle condition. The contour has zero-mean coordinates. Reference quantities `S,mu,nu,Psi,a1` are calculated on the inflated contour, before epsilon scaling. Thresholds are `Lambda_n=pi^2 n^2/(4b^2)`.

For an anticlockwise parametrization the source normal `n=(-Y',X')/|r'|` points into the obstacle, outward from the fluid exterior. The source defines `Delta Psi=0`, `partial_n Psi=n2`, zero far-field constant, and `Psi=-mu y/r^2-nu x/r^2+O(r^-2)`. Its (2.3) is `mu=(S+integral n2 Psi ds)/(2pi)`, `nu=integral n1 Psi ds/(2pi)`. These signs are consistent: a unit circle gives `Psi=-y/r^2`, `mu=1`, `nu=0`, `S=pi`. Reversing the normal requires changing every derivative convention coherently. Homogeneous Neumann data alone cannot detect an incorrectly reversed orientation.

Linton 1998 (2.3)-(2.4) uses `(Delta+k^2)G=+delta`, `G=-i H0^(1)/4`, outgoing waves for time convention `exp(-i omega t)`. Its singular part is `Y0(kr)/4`, positive `(2pi)^-1 log r`. The draft's Green representation and double-layer jump signs are compatible with this convention: `u=D u`, `u/2=K u` on the obstacle with the fluid-outward normal. The diagonal curvature term is `-kappa/(4pi)` before multiplying by arclength, with positive curvature for anticlockwise contours. Do not mix this with the widespread `+i H0/4`, `-delta` convention.

The strip image formula is Linton (4.5), translated from width `2b`: `G_D=G_periodic^(4b)(x-xi,y-eta)-G_periodic^(4b)(x-xi,y+eta+2b)`. The auxiliary period is **4b**, not 2b. At beta=0 its arguments may be reduced exactly modulo 4b before evaluating the local harmonic expansion. This changes coordinates only; it is not a lattice truncation approximation. Ordinary and regularized calls must use the same reduction, preserving the physical source singularity. The expansion (2.46) has local domain `r<d`; reduction is needed especially for the reflected separation. Regularization removes `Y0/4`, leaving `-i J0/4` for the free-space term, while retaining the full image term.

## Theorem meanings and limits

| Source | Exact scientific claim | What a numerical paper may claim |
|---|---|---|
| 2.1, (2.4)-(2.5), p.4 | For sufficiently small epsilon, a critical height `a*=a0*+O(epsilon)` separates a unique discrete trapped mode for `a>a*` from absence for `a<=a*`; `a0*=(2b/pi)atan sqrt(S/(2pi mu))`. `sigma=epsilon^2 pi^2/(4b^3)[pi mu sin^2(alpha/2)-(S/2)cos^2(alpha/2)]+O(epsilon^3 log epsilon)`, `alpha=pi a/b`, `k^2=Lambda1-sigma^2`. | Fixed-contour existence/counting and finite-grid height transition evidence. Not an epsilon-uniform theorem, analyticity proof, or global absence from a failed search. The source statement uses the positive displacement convention; do not silently assert it for negative a without reflection/geometry analysis. |
| 2.3(i), p.4 | No embedded eigenvalues away from Lambda2 by a fixed squared-eigenvalue margin sigma0>0 for sufficiently small epsilon. | A finite screening experiment is only partial support and cannot prove the quantified theorem. |
| 2.3(ii), p.4 | If horizontal dipole nu!=0, no eigenvalues on `[Lambda1,Lambda2)` for sufficiently small epsilon. | X/Y symmetric examples both have nu=0; they do not test this absence criterion. |
| 2.3(iii), (2.6), pp.4-5 | X-axis symmetry requires a=0, X even/Y odd, nu=0. A unique embedded eigenvalue has `sigma=epsilon^2 pi^3 mu/b^3+O(epsilon^3 log epsilon)`. Odd-y restriction removes the first channel and makes this interval discrete in that sector. | Beyn counting in the odd-y sector is justified for the finite analytic matrix problem; full reconstructed fields add physical evidence. Sector uniqueness is not full-space uniqueness. |
| 2.3(iv), (2.7), pp.4-5 | Y-axis symmetry requires X odd/Y even, nu=0. The displacement is specifically tuned `a=epsilon a1+O(epsilon^2)`, with `a1= integral (Y X'-3X Y')(Y-Psi)dt/[2(2S+pi mu)]`. Sigma follows (2.6). Even/odd-x restriction retains the radiation channel. | Full-operator singularity plus complex open-channel cancellation and decay support a tuned BIC. A sampled scan or Beyn rank of the embedded full problem alone does not establish a full-space spectral count. |
| Remark 2.5 | Analyticity in epsilon and epsilon log epsilon; symmetry-sector distinction. | Analyticity is theoretical, not demonstrated by finitely many computed points. |
| Remark 2.6 | Uniqueness depends on obstacle smallness and may fail for larger obstacles. | Do not extend uniqueness indiscriminately to every simulated epsilon. |

For sigma with nonzero leading coefficient, a gap of order epsilon^4 follows algebraically. Near the critical height the leading coefficient vanishes, so that generic order statement should be qualified. The remainder is `O(epsilon^3 log epsilon)`, not an ordinary analytic-power expansion in epsilon alone.

## Example 5.4 inconsistency: verified algebra, unresolved source correction

Rendered and inspected source p.24. The printed (5.47) is `a1=[integral (Y X'-3X Y') f dt]/(2S+pi mu)`, where `f=L0 Y=(Y-Psi)/2`. The same page supplies `X=sin t-beta sin(2t)/2`, `Y=-cos t+beta cos(2t)/2`, `f=-cos t+(beta/2)(cos(2t)-1/2)+O(beta^2)` and denominator `3pi+O(beta)`.

Direct differentiation/integration of those **printed expressions** gives numerator `-(3pi/2) beta+O(beta^2)` and hence `a1=-beta/2+O(beta^2)`, not the printed `-beta/12+O(beta^2)`. The zero-order numerator vanishes, so the first-order denominator perturbation cannot restore the factor six. The symbolic calculation was independently executed with SymPy. This establishes an internal inconsistency of the supplied preprint; it does not establish which source equation the authors intended to correct. Numerical integration of (2.7), with controlled MFS sign/gauge and geometry, remains the appropriate principal value. Preserve the disagreement explicitly and never substitute the printed example to obtain a desired BIC. The peer-reviewed article text has not been compared line by line to this preprint; its DOI metadata alone does not settle whether the final version fixes the issue.

## Current draft changes required

The draft's reference-contour definitions and principal Green signs are usable. Its methods/results architecture, however, is legacy: it describes M32/M40 refinement and singular-value/band scans, lacks the actual current Beyn moment/rank gates, and describes Y placement as `-epsilon beta/12` at beta=.10. The pinned latest Y run `2026-10-13_20-16-45_full_63388d81f58d/run_manifest.json` sets `shape_beta=0.5`, `bem_M=24`, `harmonic_order=20`; it computes (2.7) and tunes a second-order correction. The pinned latest X run `2026-10-13_21-23-19_full_fe445a6af342/run_manifest.json` has the same beta/M/harmonic settings. The earlier audit text incorrectly said beta=.45; this correction follows the authoritative manifests. That entire experiment description and its tables/figures must be replaced, not cosmetically revised. Current full-space uniqueness is not established for X or Y. Do not transfer the theoretical uniqueness statement into a numerical certification claim. Explicitly state that 2.3(i)-(ii) are not fully tested by the symmetric examples.

## Bibliography audit

Verified public metadata at Crossref DOI records (URLs included in JSON), plus publisher Linton page [Springer](https://link.springer.com/article/10.1023/A:1004377501747). Relevant bibliography corrections/additions:

- Zhevandrov et al. 2025: title/authors/DOI match; Acta Applicandae Mathematicae **196, article 8**, issue 1. `pages={8}` represents an article number, not a page range. Local PDF is a March 4 preprint, public record March 7; make provenance explicit.
- Linton/Evans 1992: JFM **245,349-365**, DOI `10.1017/S0022112092000491`, verified from local title page and DOI. This is lineage for integral trapped-mode methods, not direct evidence that all its wall conditions equal this mixed quantum problem.
- Linton 1998: add DOI `10.1023/A:1004377501747`. Public Springer/Crossref citation says **33,377-401**; local PDF first-page masthead says377-402. Record the discrepancy and prefer publisher's current citation metadata in final bibliography. The local PDF also ends on printed page402 (references), independently confirmed by text extraction. Thus this is an actual mismatch between supplied print pagination and current publisher metadata; document it rather than claiming resolution.
- **Missing essential reference:** Wolf-Jürgen Beyn, *An integral method for solving nonlinear eigenvalue problems*, Linear Algebra and its Applications **436(10),3839-3863 (2012)**, DOI `10.1016/j.laa.2011.03.030`. Verified [publisher](https://www.sciencedirect.com/science/article/pii/S0024379511002540), [author preprint](https://arxiv.org/abs/1003.1580). Its analytic matrix problem and contour/invertibility assumptions must be stated, not just cited.
- Sargent/Mestel: **84(2),312-344,2019**, DOI `10.1093/imamat/hxy060`; online-first2018 does not change issue year2019.
- Evans/Levitin/Vassiliev 1994, Khallaf/Parnovski/Vassiliev2000, Harari/Patlashenko/Givoli1998: listed title/authors/volume/pages/DOIs match depositor records.
- Ihlenburg1998 book, Pozrikidis2002 book: title/year/DOI match DOI records. BibTeX author/editor/ISBN refinements are not essential if these peripheral FEM/software discussions are removed.
- Kang et al.2023, Cao et al.2021, Martí-Sabaté et al.2024: names/title/DOI/journal volume and pages/article numbers match records. Last two use article identifiers104502 and122, not ordinary page ranges. Photonic/acoustic/elastic BIC literature is context, not evidence for this mixed-boundary quantum operator.
- Koshelev et al.: title/authors/66(5),494-517 match DOI deposit; Crossref date records2021-12 (DOI/manuscript timing) despite listed issue2023. If retained, verify issue-year on publisher instead of treating deposit timing as issue year.
- Fahy/Walker: Crossref API request failed; do not declare metadata verified. Sargent thesis and Klaren MIT thesis have not been checked against repository landing pages. These three are optional context; remove them unless scientific relevance and authoritative metadata are established.

No reference should be described as confirming the present numerical values without reading the corresponding content. An unused reference may be removed rather than retained to simulate a broad literature review. The retained core bibliography can be short: source theory, Green evaluation, integral trapped-mode lineage, symmetry existence theorem, Beyn method, and precisely relevant BIC context.

## Optional scientific-skill citation

Independently verified against the [primary arXiv record](https://arxiv.org/abs/2609.00065), accessed 2026-10-14: Timothy Kassis, Vinayak Agarwal, Yuhuan He, Darshil Patel, Aubrey M. Brueckner; *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents* (2026), arXiv:2609.00065, version2 dated September2,2026, DOI10.48550/arXiv.2609.00065. This documents the procedural skill library; it provides no evidence for waveguide physics or the numerical theorem claims. Include only if acknowledging tools/procedural assistance is appropriate. No journal publication is established by this record.

```bibtex
@misc{Kassis2026ScientificSkills,
  author = {Kassis, Timothy and Agarwal, Vinayak and He, Yuhuan and Patel, Darshil and Brueckner, Aubrey M.},
  title = {Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents},
  year = {2026},
  eprint = {2609.00065},
  archivePrefix = {arXiv},
  primaryClass = {cs.CL},
  doi = {10.48550/arXiv.2609.00065},
  url = {https://arxiv.org/abs/2609.00065v2}
}
```
