---
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python first_theorem_3.py
=== Theorem 2.1 validation with Beyn contour discovery ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
epsilons = (0.05, 0.07, 0.09, 0.11)
Beyn ellipse: real endpoints=[0.00010000, 1.57079623], center=0.78544811, rx=0.78534811, ry=3.000e-02
Beyn settings: M=24, Nq=96, probe_dim=8
local verification M = (16, 24, 32, 40, 48)

Checking complex-kb support required by Beyn...
  complex-kb assembly: PASS

=== epsilon=0.050 ===
  predicted kb = 1.570768582119
  predicted sigma = 9.33604312e-03
  predicted cutoff gap = 2.774e-05
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  Beyn estimated enclosed rank = 2
  singular values of S0: [2.871e-04, 5.023e-09, 2.555e-12, 1.989e-12, 1.658e-12, 1.543e-12, 1.294e-12, 1.057e-12]
  max contour linear-solve residual = 3.751e-16
  raw Beyn eigenvalues:
    [0] +1.570963444967 -3.207e-09i  other
    [1] +1.571856681582 +1.240e-04i  other
  clustered near-real Beyn candidates = 0

  --- Beyn + local-SVD validation result ---
  Beyn estimated contour rank    = 2
  near-real Beyn candidates      = 0
  resolved discrete modes        = 0
  exactly one resolved mode      = FAIL
  kb asymptotic                  = 1.570768582119
  kb BEM                         = nan
  sigma asym                     = 9.33604312e-03
  sigma BEM                      = nan
  final sigma_min(A)             = nan
  final minimum drop             = nan
  final mesh change in sigma     = nan%
  relative sigma error           = nan%
  asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.070 ===
  predicted kb = 1.570689740172
  predicted sigma = 1.82986445e-02
  predicted cutoff gap = 1.066e-04
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  Beyn estimated enclosed rank = 2
  singular values of S0: [5.952e-04, 9.744e-09, 3.313e-12, 2.808e-12, 2.341e-12, 2.069e-12, 1.675e-12, 1.315e-12]
  max contour linear-solve residual = 3.922e-16
  raw Beyn eigenvalues:
    [0] +1.570824196080 -1.781e-09i  other
    [1] +1.571825304560 +5.174e-05i  other
  clustered near-real Beyn candidates = 0

  --- Beyn + local-SVD validation result ---
  Beyn estimated contour rank    = 2
  near-real Beyn candidates      = 0
  resolved discrete modes        = 0
  exactly one resolved mode      = FAIL
  kb asymptotic                  = 1.570689740172
  kb BEM                         = nan
  sigma asym                     = 1.82986445e-02
  sigma BEM                      = nan
  final sigma_min(A)             = nan
  final minimum drop             = nan
  final mesh change in sigma     = nan%
  relative sigma error           = nan%
  asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.090 ===
  predicted kb = 1.570505049848
  predicted sigma = 3.02487797e-02
  predicted cutoff gap = 2.913e-04
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  Beyn estimated enclosed rank = 2
  singular values of S0: [1.293e-03, 1.596e-08, 4.145e-12, 3.921e-12, 2.842e-12, 2.485e-12, 2.292e-12, 1.902e-12]
  max contour linear-solve residual = 4.585e-16
  raw Beyn eigenvalues:
    [0] +1.570605175103 -2.175e-09i  real-candidate
    [1] +1.571739669504 +2.815e-05i  other
  clustered near-real Beyn candidates = 1
  candidate 1: Beyn kb=1.570605175103 -2.175e-09i, local bracket=[1.568605175103, 1.570796226795]
    M=16: kb=1.570529966364, sigma_BEM=2.89261996e-02, sv_min=4.329e-07, drop=4.61e+05, mesh_change=--, interior=yes
    M=24: kb=1.570529969142, sigma_BEM=2.89260488e-02, sv_min=3.654e-07, drop=5.46e+05, mesh_change=0.001%, interior=yes
    M=32: kb=1.570529969807, sigma_BEM=2.89260127e-02, sv_min=3.495e-07, drop=5.71e+05, mesh_change=0.000%, interior=yes
    M=40: kb=1.570529970042, sigma_BEM=2.89259999e-02, sv_min=3.440e-07, drop=5.80e+05, mesh_change=0.000%, interior=yes
    M=48: kb=1.570529970146, sigma_BEM=2.89259943e-02, sv_min=3.411e-07, drop=5.85e+05, mesh_change=0.000%, interior=yes
    resolved=YES

  --- Beyn + local-SVD validation result ---
  Beyn estimated contour rank    = 2
  near-real Beyn candidates      = 1
  resolved discrete modes        = 1
  exactly one resolved mode      = PASS
  kb asymptotic                  = 1.570505049848
  kb BEM                         = 1.570529970146
  sigma asym                     = 3.02487797e-02
  sigma BEM                      = 2.89259943e-02
  final sigma_min(A)             = 3.411e-07
  final minimum drop             = 5.85e+05
  final mesh change in sigma     = 0.000%
  relative sigma error           = 4.373%
  asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.110 ===
  predicted kb = 1.570146262336
  predicted sigma = 4.51864487e-02
  predicted cutoff gap = 6.501e-04
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  Beyn estimated enclosed rank = 2
  singular values of S0: [2.762e-03, 2.363e-08, 5.007e-12, 4.243e-12, 3.526e-12, 3.202e-12, 2.407e-12, 1.745e-12]
  max contour linear-solve residual = 6.575e-16
  raw Beyn eigenvalues:
    [0] +1.570274790691 -2.070e-09i  real-candidate
    [1] +1.571777193693 -8.875e-05i  other
  clustered near-real Beyn candidates = 1
  candidate 1: Beyn kb=1.570274790691 -2.070e-09i, local bracket=[1.568274790691, 1.570796226795]
    M=16: kb=1.570232601205, sigma_BEM=4.20794295e-02, sv_min=5.469e-07, drop=4.45e+05, mesh_change=--, interior=yes
    M=24: kb=1.570232608897, sigma_BEM=4.20791425e-02, sv_min=4.655e-07, drop=5.23e+05, mesh_change=0.001%, interior=yes
    M=32: kb=1.570232610734, sigma_BEM=4.20790739e-02, sv_min=4.453e-07, drop=5.46e+05, mesh_change=0.000%, interior=yes
    M=40: kb=1.570232611386, sigma_BEM=4.20790496e-02, sv_min=4.390e-07, drop=5.54e+05, mesh_change=0.000%, interior=yes
    M=48: kb=1.570232611673, sigma_BEM=4.20790388e-02, sv_min=4.361e-07, drop=5.58e+05, mesh_change=0.000%, interior=yes
    resolved=YES

  --- Beyn + local-SVD validation result ---
  Beyn estimated contour rank    = 2
  near-real Beyn candidates      = 1
  resolved discrete modes        = 1
  exactly one resolved mode      = PASS
  kb asymptotic                  = 1.570146262336
  kb BEM                         = 1.570232611673
  sigma asym                     = 4.51864487e-02
  sigma BEM                      = 4.20790388e-02
  final sigma_min(A)             = 4.361e-07
  final minimum drop             = 5.58e+05
  final mesh change in sigma     = 0.000%
  relative sigma error           = 6.877%
  asymptotic accuracy <= 5.0%: FAIL

=== FINAL SUMMARY ===
Global discovery: Beyn contour method
Local certification: SVD minima + BEM mesh refinement
Interval enclosed: 0 < k^2 < Lambda_1, excluding tiny endpoint margins
Exactly one resolved discrete mode: 2/4
Leading asymptotic sigma within 5.0%: 1/4
  epsilon=0.050: Beyn-rank=2, near-real=0, resolved=0, unique=FAIL, sigma_error=nan%
  epsilon=0.070: Beyn-rank=2, near-real=0, resolved=0, unique=FAIL, sigma_error=nan%
  epsilon=0.090: Beyn-rank=2, near-real=1, resolved=1, unique=PASS, sigma_error=4.373%
  epsilon=0.110: Beyn-rank=2, near-real=1, resolved=1, unique=PASS, sigma_error=6.877%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_1_beyn_validation
