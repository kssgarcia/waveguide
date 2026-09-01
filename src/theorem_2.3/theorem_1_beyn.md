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

---
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python first_theorem_4.py
=== Theorem 2.1 validation with Beyn contour discovery (v2) ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
epsilons = (0.05, 0.07, 0.09, 0.11)
Beyn ellipse: real endpoints=[0.00010000, 1.57078633], center=0.78544316, rx=0.78534316, ry=3.000e-02
Beyn settings: M=24, Nq-levels=(96, 192, 384), probe_dim=8
rank diagnostics: rel_tol=1.0e-08, gap_threshold=1.0e+03
local verification M = (16, 24, 32, 40, 48)

Checking complex-kb support required by Beyn...
  complex-kb assembly: PASS

=== epsilon=0.050 ===
  predicted kb = 1.570768582119
  predicted sigma = 9.33604312e-03
  predicted cutoff gap = 2.774e-05
  Beyn quadrature-convergence study:
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.719e+04
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.825e-04, 4.939e-09, 2.411e-12, 1.946e-12, 1.675e-12, 1.460e-12, 1.204e-12, 9.691e-13]
  max contour linear-solve residual = 3.779e-16
  raw Beyn eigenvalues:
    [0] +1.570963315921 -2.186e-09i  outside-contour/band
    Nq=96: rank=1, candidates=0, S0-change=--, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=192, probe_dim=8
    contour node  32/192
    contour node  64/192
    contour node  96/192
    contour node 128/192
    contour node 160/192
    contour node 192/192
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.688e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.285e-04, 1.353e-09, 1.919e-12, 1.732e-12, 1.499e-12, 1.274e-12, 9.016e-13, 7.497e-13]
  max contour linear-solve residual = 4.364e-16
  raw Beyn eigenvalues:
    [0] +1.570834970811 -6.507e-10i  outside-contour/band
    Nq=192: rank=1, candidates=0, S0-change=2.362e-01, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=384, probe_dim=8
    contour node  64/384
    contour node 128/384
    contour node 192/384
    contour node 256/384
    contour node 320/384
    contour node 384/384
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.058e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.022e-04, 3.997e-10, 1.762e-12, 1.593e-12, 1.296e-12, 1.215e-12, 1.119e-12, 9.452e-13]
  max contour linear-solve residual = 5.810e-16
  raw Beyn eigenvalues:
    [0] +1.570791155740 -1.601e-10i  outside-contour/band
    Nq=384: rank=1, candidates=0, S0-change=1.302e-01, candidate-shift=--
  final clustered near-real Beyn candidates = 0

  --- Beyn-v2 + local-SVD validation result ---
  final Beyn quadrature Nq       = 384
  final Beyn estimated rank      = 1
  Beyn rank stable (last 2 Nq)   = YES
  final near-real candidates     = 0
  locally resolved modes         = 0
  one-mode count supported      = FAIL
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
  Beyn quadrature-convergence study:
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=6.119e+04
  Beyn estimated enclosed rank = 1
  singular values of S0: [5.862e-04, 9.580e-09, 3.160e-12, 2.506e-12, 2.434e-12, 2.079e-12, 1.607e-12, 1.439e-12]
  max contour linear-solve residual = 4.050e-16
  raw Beyn eigenvalues:
    [0] +1.570823960741 -2.099e-09i  outside-contour/band
    Nq=96: rank=1, candidates=0, S0-change=--, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=192, probe_dim=8
    contour node  32/192
    contour node  64/192
    contour node  96/192
    contour node 128/192
    contour node 160/192
    contour node 192/192
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.195e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [5.763e-04, 2.626e-09, 2.702e-12, 2.194e-12, 1.954e-12, 1.830e-12, 1.431e-12, 1.133e-12]
  max contour linear-solve residual = 6.219e-16
  raw Beyn eigenvalues:
    [0] +1.570730866534 +5.850e-11i  real-candidate
    Nq=192: rank=1, candidates=1, S0-change=1.716e-02, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=384, probe_dim=8
    contour node  64/384
    contour node 128/384
    contour node 192/384
    contour node 256/384
    contour node 320/384
    contour node 384/384
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.232e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [6.391e-04, 7.764e-10, 2.464e-12, 2.334e-12, 1.797e-12, 1.505e-12, 1.315e-12, 9.538e-13]
  max contour linear-solve residual = 9.130e-16
  raw Beyn eigenvalues:
    [0] +1.570704461001 +2.110e-10i  real-candidate
    Nq=384: rank=1, candidates=1, S0-change=9.829e-02, candidate-shift=2.641e-05
  final clustered near-real Beyn candidates = 1
  candidate 1: Beyn kb=1.570704461001 +2.110e-10i, local bracket=[1.568704461001, 1.570786326795]
    M=16: kb=1.570694867809, sigma_BEM=1.78530812e-02, sv_min=1.668e-06, drop=6.51e+04, mesh_change=--, interior=yes
    M=24: kb=1.570694869053, sigma_BEM=1.78529718e-02, sv_min=2.158e-06, drop=5.03e+04, mesh_change=0.001%, interior=yes
    M=32: kb=1.570694869347, sigma_BEM=1.78529460e-02, sv_min=2.272e-06, drop=4.78e+04, mesh_change=0.000%, interior=yes
    M=40: kb=1.570694869451, sigma_BEM=1.78529368e-02, sv_min=2.313e-06, drop=4.69e+04, mesh_change=0.000%, interior=yes
    M=48: kb=1.570694869497, sigma_BEM=1.78529328e-02, sv_min=2.331e-06, drop=4.66e+04, mesh_change=0.000%, interior=yes
    resolved=YES

  --- Beyn-v2 + local-SVD validation result ---
  final Beyn quadrature Nq       = 384
  final Beyn estimated rank      = 1
  Beyn rank stable (last 2 Nq)   = YES
  final near-real candidates     = 1
  locally resolved modes         = 1
  one-mode count supported      = PASS
  kb asymptotic                  = 1.570689740172
  kb BEM                         = 1.570694869497
  sigma asym                     = 1.82986445e-02
  sigma BEM                      = 1.78529328e-02
  final sigma_min(A)             = 2.331e-06
  final minimum drop             = 4.66e+04
  final mesh change in sigma     = 0.000%
  relative sigma error           = 2.436%
  asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.090 ===
  predicted kb = 1.570505049848
  predicted sigma = 3.02487797e-02
  predicted cutoff gap = 2.913e-04
  Beyn quadrature-convergence study:
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.144e+04
  Beyn estimated enclosed rank = 1
  singular values of S0: [1.278e-03, 1.569e-08, 4.657e-12, 3.981e-12, 3.622e-12, 2.643e-12, 2.399e-12, 1.756e-12]
  max contour linear-solve residual = 5.174e-16
  raw Beyn eigenvalues:
    [0] +1.570604821989 -1.567e-09i  real-candidate
    Nq=96: rank=1, candidates=1, S0-change=--, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=192, probe_dim=8
    contour node  32/192
    contour node  64/192
    contour node  96/192
    contour node 128/192
    contour node 160/192
    contour node 192/192
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=3.382e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [1.454e-03, 4.300e-09, 3.737e-12, 3.154e-12, 2.598e-12, 2.342e-12, 1.972e-12, 1.692e-12]
  max contour linear-solve residual = 8.628e-16
  raw Beyn eigenvalues:
    [0] +1.570548000754 -4.111e-10i  real-candidate
    Nq=192: rank=1, candidates=1, S0-change=1.214e-01, candidate-shift=5.682e-05
  Beyn contour discovery: M=24, quadrature=384, probe_dim=8
    contour node  64/384
    contour node 128/384
    contour node 192/384
    contour node 256/384
    contour node 320/384
    contour node 384/384
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.260e+06
  Beyn estimated enclosed rank = 1
  singular values of S0: [1.601e-03, 1.271e-09, 3.742e-12, 3.386e-12, 2.903e-12, 2.657e-12, 2.186e-12, 1.895e-12]
  max contour linear-solve residual = 1.385e-15
  raw Beyn eigenvalues:
    [0] +1.570534810125 -3.545e-10i  real-candidate
    Nq=384: rank=1, candidates=1, S0-change=9.179e-02, candidate-shift=1.319e-05
  final clustered near-real Beyn candidates = 1
  candidate 1: Beyn kb=1.570534810125 -3.545e-10i, local bracket=[1.568534810125, 1.570786326795]
    M=16: kb=1.570529958290, sigma_BEM=2.89266380e-02, sv_min=3.453e-06, drop=4.75e+04, mesh_change=--, interior=yes
    M=24: kb=1.570529959047, sigma_BEM=2.89265969e-02, sv_min=4.141e-06, drop=3.96e+04, mesh_change=0.000%, interior=yes
    M=32: kb=1.570529959225, sigma_BEM=2.89265873e-02, sv_min=4.307e-06, drop=3.81e+04, mesh_change=0.000%, interior=yes
    M=40: kb=1.570529959287, sigma_BEM=2.89265839e-02, sv_min=4.366e-06, drop=3.75e+04, mesh_change=0.000%, interior=yes
    M=48: kb=1.570529982676, sigma_BEM=2.89253140e-02, sv_min=4.345e-06, drop=3.77e+04, mesh_change=0.004%, interior=yes
    resolved=YES

  --- Beyn-v2 + local-SVD validation result ---
  final Beyn quadrature Nq       = 384
  final Beyn estimated rank      = 1
  Beyn rank stable (last 2 Nq)   = YES
  final near-real candidates     = 1
  locally resolved modes         = 1
  one-mode count supported      = PASS
  kb asymptotic                  = 1.570505049848
  kb BEM                         = 1.570529982676
  sigma asym                     = 3.02487797e-02
  sigma BEM                      = 2.89253140e-02
  final sigma_min(A)             = 4.345e-06
  final minimum drop             = 3.77e+04
  final mesh change in sigma     = 0.004%
  relative sigma error           = 4.375%
  asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.110 ===
  predicted kb = 1.570146262336
  predicted sigma = 4.51864487e-02
  predicted cutoff gap = 6.501e-04
  Beyn quadrature-convergence study:
  Beyn contour discovery: M=24, quadrature=96, probe_dim=8
    contour node  16/96
    contour node  32/96
    contour node  48/96
    contour node  64/96
    contour node  80/96
    contour node  96/96
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.181e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.744e-03, 2.323e-08, 4.685e-12, 4.382e-12, 3.527e-12, 3.118e-12, 2.955e-12, 2.141e-12]
  max contour linear-solve residual = 6.182e-16
  raw Beyn eigenvalues:
    [0] +1.570274341561 -1.955e-09i  real-candidate
    Nq=96: rank=1, candidates=1, S0-change=--, candidate-shift=--
  Beyn contour discovery: M=24, quadrature=192, probe_dim=8
    contour node  32/192
    contour node  64/192
    contour node  96/192
    contour node 128/192
    contour node 160/192
    contour node 192/192
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.473e+05
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.848e-03, 6.368e-09, 3.934e-12, 3.481e-12, 3.111e-12, 2.573e-12, 2.288e-12, 1.982e-12]
  max contour linear-solve residual = 9.857e-16
  raw Beyn eigenvalues:
    [0] +1.570243629925 -3.585e-10i  real-candidate
    Nq=192: rank=1, candidates=1, S0-change=3.658e-02, candidate-shift=3.071e-05
  Beyn contour discovery: M=24, quadrature=384, probe_dim=8
    contour node  64/384
    contour node 128/384
    contour node 192/384
    contour node 256/384
    contour node 320/384
    contour node 384/384
  rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.505e+06
  Beyn estimated enclosed rank = 1
  singular values of S0: [2.832e-03, 1.882e-09, 4.634e-12, 4.044e-12, 3.593e-12, 3.053e-12, 2.442e-12, 1.875e-12]
  max contour linear-solve residual = 1.421e-15
  raw Beyn eigenvalues:
    [0] +1.570235883656 -1.913e-11i  real-candidate
    Nq=384: rank=1, candidates=1, S0-change=5.667e-03, candidate-shift=7.746e-06
  final clustered near-real Beyn candidates = 1
  candidate 1: Beyn kb=1.570235883656 -1.913e-11i, local bracket=[1.568235883656, 1.570786326795]
    M=16: kb=1.570232590865, sigma_BEM=4.20798153e-02, sv_min=1.660e-06, drop=1.29e+05, mesh_change=--, interior=yes
    M=24: kb=1.570232598950, sigma_BEM=4.20795136e-02, sv_min=1.658e-06, drop=1.30e+05, mesh_change=0.001%, interior=yes
    M=32: kb=1.570232600882, sigma_BEM=4.20794415e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
    M=40: kb=1.570232601567, sigma_BEM=4.20794160e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
    M=48: kb=1.570232601869, sigma_BEM=4.20794047e-02, sv_min=1.656e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
    resolved=YES

  --- Beyn-v2 + local-SVD validation result ---
  final Beyn quadrature Nq       = 384
  final Beyn estimated rank      = 1
  Beyn rank stable (last 2 Nq)   = YES
  final near-real candidates     = 1
  locally resolved modes         = 1
  one-mode count supported      = PASS
  kb asymptotic                  = 1.570146262336
  kb BEM                         = 1.570232601869
  sigma asym                     = 4.51864487e-02
  sigma BEM                      = 4.20794047e-02
  final sigma_min(A)             = 1.656e-06
  final minimum drop             = 1.30e+05
  final mesh change in sigma     = 0.000%
  relative sigma error           = 6.876%
  asymptotic accuracy <= 5.0%: FAIL

=== FINAL SUMMARY ===
Global discovery: Beyn contour method with Nq convergence study
Local certification: SVD minima + BEM mesh refinement
Interval enclosed: 0 < k^2 < Lambda_1, excluding endpoint margins
Stable Beyn rank=1 + one locally resolved mode: 3/4
Leading asymptotic sigma within 5.0%: 2/4
  epsilon=0.050: Nq=384, rank=1, rank_stable=yes, near-real=0, resolved=0, one-mode=FAIL, sigma_error=--
  epsilon=0.070: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=2.436%
  epsilon=0.090: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=4.375%
  epsilon=0.110: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=6.876%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_1_beyn_v2_validation
