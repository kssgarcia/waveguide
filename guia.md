# Updated English Speech

## 1

Good morning. My name is Kevin Sepúlveda García, from Universidad Militar Nueva Granada.

Today I will present our work, _Trapped Modes and BICs in Perturbed Waveguides: A Numerical Validation_.

The main objective is to compare asymptotic predictions for trapped modes with independent numerical computations obtained using the Boundary Element Method, or BEM.

We consider both discrete trapped modes and embedded trapped modes, also known as Bound States in the Continuum, or BICs.

**Transition:**
Let me first explain the difference between these two types of localized states.

---

## 2

We consider an infinite waveguide, where waves can in principle propagate indefinitely in the longitudinal direction.

However, a small localized perturbation — in our case, a small rigid obstacle — can generate a mode whose energy remains confined near the obstacle.

We study two cases.

A **discrete trapped mode** has an eigenvalue below the first propagation threshold. In this spectral region, all longitudinal components are evanescent, so localization is relatively natural.

A **BIC** is more subtle. Its eigenvalue lies inside the continuous spectrum, where a propagating transverse mode is already available, but the state nevertheless remains localized.

So our main question is: **how accurately do the numerical eigenvalues reproduce the asymptotic predictions?**

**Transition:**
We begin with the corresponding Helmholtz eigenvalue problem.

---

## 3

We consider a two-dimensional infinite waveguide with walls at

$$
y=\pm b,
$$

and a small obstacle with boundary \(\gamma\).

The field satisfies

$$
-\Delta u=k^2u.
$$

We impose homogeneous Dirichlet conditions on the waveguide walls,

$$
u=0,
$$

and a homogeneous Neumann condition on the rigid obstacle,

$$
\frac{\partial u}{\partial n}=0.
$$

Finally,

$$
u(x,y)\rightarrow0
\qquad\text{as}\qquad |x|\rightarrow\infty,
$$

which expresses spatial localization.

So we look for special values of \(k\) for which a non-trivial localized solution exists.

**Transition:**
The position of \(k^2\) relative to the waveguide cut-offs determines the spectral regime.

---

## 4

The transverse cut-off values are

$$
\Lambda_n=
\left(\frac{n\pi}{2b}\right)^2.
$$

The first one, \(\Lambda_1\), marks the beginning of the continuous spectrum.

For

$$
0<k^2<\Lambda_1,
$$

all transverse components are evanescent. This is the spectral region of the **discrete trapped mode**.

For

$$
\Lambda_1<k^2<\Lambda_2,
$$

the first transverse mode becomes propagating. This is the region where we study the **embedded trapped mode**.

This distinction is important numerically. Below \(\Lambda_1\), finding the eigenvalue is enough to obtain longitudinal decay. Inside the continuous spectrum, a nearly singular system is not enough: we must also verify that the propagating component is cancelled.

**Transition:**
The asymptotic theory gives explicit predictions for where these eigenvalues should occur.

---

## 5

For the discrete mode, the theory predicts

$$
k^2=\Lambda_1-\sigma^2.
$$

For the embedded mode,

$$
k^2=\Lambda_2-\sigma^2,
\qquad
\Lambda_1<k^2<\Lambda_2.
$$

The corresponding expressions for \(\sigma\) are asymptotic formulas derived for

$$
\varepsilon\ll1.
$$

For the embedded case, additional geometric symmetry conditions are required. We therefore consider one configuration symmetric with respect to the \(x\)-axis and another symmetric with respect to the \(y\)-axis.

Our purpose is not to prove these formulas numerically. We want to determine how accurately the independently computed eigenvalues follow them as \(\varepsilon\) changes.

**Transition:**
For that independent numerical computation, we use BEM.

---

## 6

BEM is especially convenient because the waveguide is infinite.

With a domain method, we would have to truncate the guide and introduce an additional treatment such as a PML or a Dirichlet-to-Neumann condition.

Instead, we use a Green function specifically adapted to the waveguide.

It already incorporates the wall conditions and the unbounded longitudinal direction.

Therefore, we only need to discretize the boundary of the obstacle.

**Transition:**
The next two slides show first how this Green function is constructed, and then how it enters the BEM system.

---

# 7 — Green function slide

This slide focuses on the Green function used in the computation.

The Dirichlet waveguide Green function is constructed as the difference of two periodic Green functions:

$$
G_s(P,q)
=
G_0^{4b}(X,\bar y-\bar\eta)
-
G_0^{4b}(X,\bar y+\bar\eta).
$$

This construction enforces the homogeneous Dirichlet condition on the waveguide walls.

The periodic Green function itself is evaluated through the lattice-sum representation shown here.

There is one additional numerical issue: when the source and observation points coincide, the Green function contains the usual logarithmic singularity.

For diagonal evaluation we therefore use the regularized function

$$
G_{s,\mathrm{reg}}(P,q)
=
G_s(P,q)
-
\frac14Y_0\!\left(k|P-q|\right).
$$

So the important point from this slide is that we have both the physical waveguide Green function and a regularized version used when evaluating the singular diagonal terms.

**Transition:**
With this Green function available, we can now reduce the original PDE to an equation only on the obstacle boundary.

---

# 8 — Boundary integral / BEM slide

Using Green's representation and the Neumann condition on the rigid obstacle, the field can be written as

$$
u(P)
=
\int_\gamma
u(q)
\frac{\partial G_s(P,q)}{\partial n_q}\,ds_q.
$$

When \(P\) approaches the obstacle boundary, we obtain

$$
\frac12u(p)
=
\int_\gamma
u(q)
\frac{\partial G_s(p,q)}{\partial n_q}\,ds_q.
$$

We then parameterize the complete contour with

$$
0\leq\theta<2\pi
$$

and discretize it using midpoint quadrature.

This gives the homogeneous matrix problem

$$
\boxed{
A(k)\mathbf u=0
}
$$

with

$$
A(k)=I-\frac{4\pi}{M}K^w(k).
$$

For each trial value of \(k\), we obtain a different matrix.

A trapped-mode candidate occurs when \(A(k)\) becomes nearly singular, which we detect using its smallest singular value.

Off the diagonal, derivatives of the Green function are evaluated using centred finite differences; on the diagonal, we use the regularized Green function shown on the previous slide.

**Transition:**
The next step is therefore to locate those near-singular values and determine whether they represent genuine localized modes.

---

# 9 — Numerical validation workflow

The numerical validation has several stages.

First, the asymptotic formula provides an expected location for the eigenvalue. Around that region we evaluate

$$
s_{\min}(A(k)),
$$

the smallest singular value of the BEM matrix, and refine the corresponding minimum.

But we do not accept a candidate based only on one small singular value. We also repeat the computation under boundary-mesh refinement and check that the location of the eigenvalue remains stable.

For the **discrete mode**, there is an additional step because the theorem predicts uniqueness.

After resolving the expected mode, we independently scan the complete spectral interval

$$
0<k^2<\Lambda_1.
$$

The sampling is refined near the first cut-off, where the expected eigenvalue becomes extremely close to \(\Lambda*1\). We detect every local minimum of \(s*{\min}(A(k))\) and refine any minimum outside the expected-mode region.

In our computations, the only persistent resolved minimum was the predicted one, so we found no additional discrete modes.

For the **embedded cases**, the additional test is different. Because a propagating transverse mode exists, we reconstruct the field away from the obstacle and measure its projection onto that propagating mode:

$$
R_{\mathrm{prop}}.
$$

A BIC candidate must have

$$
R_{\mathrm{prop}}\approx0.
$$

So the numerical logic is: spectral singularity, mesh stability, and then either a full-band uniqueness check for the discrete mode or radiation cancellation for the BIC.

**Transition:**
With these checks in place, we can compare the numerical values with the asymptotic predictions.

---

## 10

Here we show a representative comparison for

$$
\varepsilon=0.073.
$$

We consider three configurations: the discrete circular obstacle, the \(x\)-symmetric embedded configuration, and the \(y\)-symmetric embedded configuration.

The table compares the numerical BEM value of \(kb\) with the asymptotic prediction.

For the discrete case, the relative discrepancy is approximately

$$
4.4\times10^{-6},
$$

or about \(0.00044\%\).

For the two embedded configurations, it is approximately

$$
1.6\times10^{-4},
$$

or about \(0.016\%\).

So, for these representative small-obstacle cases, the discrepancy is below \(0.017\%\).

**Transition:**
We then repeat the comparison as the obstacle size increases.

---

## 11

These three plots show the numerical and asymptotic values of \(kb\) as \(\varepsilon\) increases.

The discrete case is shown on the left, the \(x\)-symmetric embedded case in the center, and the \(y\)-symmetric case on the right.

For small \(\varepsilon\), the numerical and asymptotic curves are almost indistinguishable.

As the obstacle becomes larger, the curves gradually separate.

This is exactly the qualitative behavior expected from a small-obstacle asymptotic expansion: the leading approximation becomes less accurate as we move away from

$$
\varepsilon=0.
$$

For the embedded cases, however, agreement in \(kb\) alone does not prove localization.

For example, in the \(y\)-symmetric computation, the propagating-mode screening criterion is satisfied for the smaller obstacle but is exceeded for the next reported value.

This does not prove that the BIC disappears; it only means that this numerical candidate no longer satisfies our adopted radiation criterion.

**Transition:**
We also verify that the representative results are stable under BEM mesh refinement.

---

## 12

Here we increase the number of boundary elements.

The computed wavenumbers change only by very small amounts, indicating that their positions are stable under boundary refinement.

For the representative \(x\)-symmetric embedded mode, we additionally obtain

$$
R_{\mathrm{prop}}
=
3.157\times10^{-6},
$$

which is well below the adopted tolerance

$$
10^{-3}.
$$

So the evidence does not come from the small singular value alone.

It comes from the combination of spectral localization, mesh stability, and — for the embedded mode — suppression of the propagating component.

**Transition:**
These numerical checks lead to the following conclusions.

---

## 13

There are three main conclusions.

First, the waveguide Green function allows us to treat the infinite domain directly, so BEM discretizes only the obstacle boundary.

Second, in the small-obstacle regime, the numerical eigenvalues show very strong agreement with the asymptotic predictions. As \(\varepsilon\) increases, the discrepancy grows in the expected way.

For the discrete theorem, the full spectral scan also found no additional resolved minima in

$$
0<k^2<\Lambda_1,
$$

providing numerical support for the predicted uniqueness.

Third, for embedded trapped modes, a nearly singular matrix is not sufficient. The numerical candidate must also remain stable under mesh refinement and exhibit negligible projection onto the propagating transverse mode.

As future work, we want to extend the analysis beyond the small-obstacle regime and investigate whether additional branches of embedded trapped modes can appear.

Thank you very much for your attention.
