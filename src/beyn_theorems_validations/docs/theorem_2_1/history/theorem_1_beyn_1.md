(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations % python first_theorem_beyn.py
Results for this execution: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-07_11-06-56_full_832cf4c8c71d
=== Theorem 2.1 numerical validation: corrected Beyn + BEM + sigma-SVD (v9) ===
b = 1.0
baseline a = 0.6
leading-order a0\* = 0.391826552031
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
detected CPUs = 10
paper workers = 5

Checking complex-kb support required by Beyn...
complex-kb assembly: PASS

=== BASELINE FIXED-M VALIDATION ===
epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
Beyn global discovery fixed M=24; Nq base=(256, 512, 1024), adaptive=(2048, 4096)
fixed BEM/SVD order M = 24
baseline completed 1/10: epsilon=0.0200
baseline completed 2/10: epsilon=0.0600
baseline completed 3/10: epsilon=0.0400
baseline completed 4/10: epsilon=0.0800
baseline completed 5/10: epsilon=0.1000
baseline completed 6/10: epsilon=0.1200
baseline completed 7/10: epsilon=0.1400
baseline completed 8/10: epsilon=0.1600
baseline completed 9/10: epsilon=0.1800
baseline completed 10/10: epsilon=0.2000

=== epsilon=0.02000000 ===
predicted kb = 1.570795616537
predicted sigma = 1.49376690e-03
predicted cutoff gap = 7.103e-07
Beyn ellipse: real=[0.000100000000, 1.570796149230], center=0.785448074615, rx=0.785348074615, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.776e-07
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 2.151e-05
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570795617339
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 2.005e-12
equilibrated vector residual = 2.006e-12
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.570795616537
kb BEM = 1.570795617909
sigma asym = 1.49376690e-03
sigma BEM = 1.49232358e-03
sigma/epsilon^2 = 3.73080894e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 4.61181606e-02
final sigma_min(A) = 8.926e-11
final sigma_max(A) = 4.453e+01
final sigma_min/sigma_max = 2.005e-12
final minimum drop = 2.52e+08 cross-M change = N/A (single fixed M)
relative sigma error = 0.097%
wall Dirichlet residual = 1.017e-09
off-grid BIE residual = 2.420e-06
off-grid BIE residual check = PASS
decay rate left/right = 1.491908e-03/1.491908e-03
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.04000000 ===
predicted kb = 1.570784962635
predicted sigma = 5.97506760e-03
predicted cutoff gap = 1.136e-05
Beyn ellipse: real=[0.000100000000, 1.570793485755], center=0.785446742877, rx=0.785346742877, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 2.841e-06
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 3.175e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570785102978
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 2.151e-11
equilibrated vector residual = 2.155e-11
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.570784962635
kb BEM = 1.570785102850
sigma asym = 5.97506760e-03
sigma BEM = 5.93809190e-03 sigma/epsilon^2 = 3.71130743e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.79486676e-01
final sigma_min(A) = 4.831e-10
final sigma_max(A) = 2.246e+01
final sigma_min/sigma_max = 2.151e-11
final minimum drop = 9.28e+07
cross-M change = N/A (single fixed M)
relative sigma error = 0.619%
wall Dirichlet residual = 7.479e-10
off-grid BIE residual = 6.501e-07
off-grid BIE residual check = PASS
decay rate left/right = 5.936045e-03/5.936044e-03
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.06000000 ===
predicted kb = 1.570738794889
predicted sigma = 1.34439021e-02
predicted cutoff gap = 5.753e-05
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 2.685e-07
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570740717207
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 1.473e-11
equilibrated vector residual = 1.479e-11
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.570738794889
kb BEM = 1.570740717474
sigma asym = 1.34439021e-02
sigma BEM = 1.32173651e-02
sigma/epsilon^2 = 3.67149031e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 3.72779673e-01
final sigma_min(A) = 2.242e-10
final sigma_max(A) = 1.521e+01
final sigma_min/sigma_max = 1.473e-11
final minimum drop = 3.48e+08
cross-M change = N/A (single fixed M)
relative sigma error = 1.685%
wall Dirichlet residual = 4.636e-10
off-grid BIE residual = 6.432e-07
off-grid BIE residual check = PASS
decay rate left/right = 1.321375e-02/1.321375e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.08000000 ===
predicted kb = 1.570614490366
predicted sigma = 2.39002704e-02
predicted cutoff gap = 1.818e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 1.157e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570626410508
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 2.011e-10
equilibrated vector residual = 2.025e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.570614490366
kb BEM = 1.570626410752
sigma asym = 2.39002704e-02
sigma BEM = 2.31036387e-02
sigma/epsilon^2 = 3.60994355e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 6.16028672e-01
final sigma_min(A) = 2.351e-09
final sigma_max(A) = 1.169e+01
final sigma_min/sigma_max = 2.011e-10
final minimum drop = 5.77e+07
cross-M change = N/A (single fixed M)
relative sigma error = 3.333%
wall Dirichlet residual = 2.820e-10
off-grid BIE residual = 1.143e-06
off-grid BIE residual check = PASS
decay rate left/right = 2.309654e-02/2.309654e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.10000000 ===
predicted kb = 1.570352353153
predicted sigma = 3.73441725e-02
predicted cutoff gap = 4.440e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 1.457e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn Aitken diagnostic kb = 1.570400353382
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 1.182e-11
equilibrated vector residual = 1.194e-11
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.570352353153
kb BEM = 1.570400353491
sigma asym = 3.73441725e-02
sigma BEM = 3.52679745e-02
sigma/epsilon^2 = 3.52679745e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 9.01681315e-01
final sigma_min(A) = 1.141e-10
final sigma_max(A) = 9.652e+00
final sigma_min/sigma_max = 1.182e-11
final minimum drop = 1.48e+09
cross-M change = N/A (single fixed M)
relative sigma error = 5.560%
wall Dirichlet residual = 2.188e-10
off-grid BIE residual = 1.785e-06
off-grid BIE residual check = PASS
decay rate left/right = 3.525558e-02/3.525558e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.12000000 ===
predicted kb = 1.569875563291
predicted sigma = 5.37756084e-02
predicted cutoff gap = 9.208e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1 Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 8.330e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570022623017
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 1.590e-10
equilibrated vector residual = 1.611e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.569875563291
kb BEM = 1.570022622889
sigma asym = 5.37756084e-02
sigma BEM = 4.92956782e-02
sigma/epsilon^2 = 3.42331099e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.22274997e+00
final sigma_min(A) = 1.329e-09
final sigma_max(A) = 8.362e+00
final sigma_min/sigma_max = 1.590e-10
final minimum drop = 1.52e+08
cross-M change = N/A (single fixed M)
relative sigma error = 8.331%
wall Dirichlet residual = 9.569e-11
off-grid BIE residual = 2.570e-06
off-grid BIE residual check = PASS
decay rate left/right = 4.927735e-02/4.927735e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.14000000 ===
predicted kb = 1.569090071990
predicted sigma = 7.31945781e-02
predicted cutoff gap = 1.706e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 9.029e-05
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.569462822882
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 2.292e-10
equilibrated vector residual = 2.331e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.569090071990
kb BEM = 1.569462822959
sigma asym = 7.31945781e-02
sigma BEM = 6.47112635e-02
sigma/epsilon^2 = 3.30159508e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.57243616e+00
final sigma_min(A) = 1.719e-09
final sigma_max(A) = 7.499e+00
final sigma_min/sigma_max = 2.292e-10
final minimum drop = 1.36e+08
cross-M change = N/A (single fixed M)
relative sigma error = 11.590%
wall Dirichlet residual = 1.591e-10
off-grid BIE residual = 3.496e-06
off-grid BIE residual check = PASS
decay rate left/right = 6.468494e-02/6.468494e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.16000000 ===
predicted kb = 1.567884413304
predicted sigma = 9.56010815e-02
predicted cutoff gap = 2.912e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 1024
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 5.009e-04
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.568706109285
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 4.529e-10
equilibrated vector residual = 4.625e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.567884413304
kb BEM = 1.568706109251
sigma asym = 9.56010815e-02
sigma BEM = 8.10076729e-02
sigma/epsilon^2 = 3.16436222e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.94416672e+00
final sigma_min(A) = 3.127e-09
final sigma_max(A) = 6.904e+00
final sigma_min/sigma_max = 4.529e-10
final minimum drop = 8.46e+07
cross-M change = N/A (single fixed M)
relative sigma error = 15.265%
wall Dirichlet residual = 1.406e-10
off-grid BIE residual = 4.565e-06
off-grid BIE residual check = PASS
decay rate left/right = 8.097060e-02/8.097060e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.18000000 ===
predicted kb = 1.566129394876
predicted sigma = 1.20995119e-01
predicted cutoff gap = 4.667e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 2048
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 1.345e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.567756530567
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 6.707e-10
equilibrated vector residual = 6.880e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.566129394876
kb BEM = 1.567756530637
sigma asym = 1.20995119e-01
sigma BEM = 9.76757949e-02
sigma/epsilon^2 = 3.01468503e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 2.33176813e+00
final sigma_min(A) = 4.352e-09
final sigma_max(A) = 6.488e+00
final sigma_min/sigma_max = 6.707e-10
final minimum drop = 6.75e+07
cross-M change = N/A (single fixed M)
relative sigma error = 19.273%
wall Dirichlet residual = 7.743e-11
off-grid BIE residual = 5.769e-06
off-grid BIE residual check = PASS
decay rate left/right = 9.762192e-02/9.762192e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.20000000 ===
predicted kb = 1.563677621758
predicted sigma = 1.49376690e-01
predicted cutoff gap = 7.119e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02

--- publication Beyn + local-sigma-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 2048
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final relative S0 change = 2.442e-06
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.566637282120
locally resolved modes = 1
mode existence supported = True
normalized vector residual = 5.584e-10
equilibrated vector residual = 5.756e-10
S0/S1 moments converged = True
uniqueness supported = supported
kb asymptotic = 1.563677621758
kb BEM = 1.566637282173
sigma asym = 1.49376690e-01
sigma BEM = 1.14231022e-01
sigma/epsilon^2 = 2.85577555e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 2.72965391e+00
final sigma_min(A) = 3.461e-09
final sigma_max(A) = 6.199e+00
final sigma_min/sigma_max = 5.584e-10
final minimum drop = 9.27e+07
cross-M change = N/A (single fixed M)
relative sigma error = 23.528%
wall Dirichlet residual = 1.129e-10
off-grid BIE residual = 7.112e-06
off-grid BIE residual check = PASS
decay rate left/right = 1.141665e-01/1.141665e-01
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== PUBLICATION NUMERICAL CONVERGENCE at epsilon=0.100 ===
common reference: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10
finite_difference_step=1e-08: M=24, kb=1.570400353409, sigma=3.52679782e-02, rel_sv=3.707e-11, delta_sigma=0.000%
finite_difference_step=1e-07: M=24, kb=1.570400353492, sigma=3.52679745e-02, rel_sv=7.891e-11, delta_sigma=0.000%
finite_difference_step=1e-06: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%
finite_difference_step=1e-05: M=24, kb=1.570400353497, sigma=3.52679743e-02, rel_sv=1.560e-10, delta_sigma=0.000%
finite_difference_step=0.0001: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.545e-10, delta_sigma=0.000%
lattice_terms=25: M=24, kb=1.570400353501, sigma=3.52679741e-02, rel_sv=1.582e-10, delta_sigma=0.000%
lattice_terms=50: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.499e-10, delta_sigma=0.000%
lattice_terms=100: M=24, kb=1.570400353497, sigma=3.52679742e-02, rel_sv=1.492e-10, delta_sigma=0.000%
lattice_terms=200: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%
harmonic_order=5: M=24, kb=1.570400353493, sigma=3.52679745e-02, rel_sv=1.423e-10, delta_sigma=0.000%
harmonic_order=10: M=24, kb=1.570400353497, sigma=3.52679742e-02, rel_sv=1.518e-10, delta_sigma=0.000%
harmonic_order=15: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%
harmonic_order=20: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%
harmonic_order=30: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%
harmonic_order=40: M=24, kb=1.570400353498, sigma=3.52679742e-02, rel_sv=1.524e-10, delta_sigma=0.000%

=== PUBLICATION CRITICAL-HEIGHT STUDY a\*(epsilon) ===
parallel policy: detected=10, fraction=50%, workers=4
[ 1/4] epsilon=0.1000: unbracketed-or-ambiguous, a_c=--
[ 2/4] epsilon=0.0600: unbracketed-or-ambiguous, a_c=--
[ 3/4] epsilon=0.0800: unbracketed-or-ambiguous, a_c=--
[ 4/4] epsilon=0.0400: unbracketed-or-ambiguous, a_c=--

=== PAPER SWEEP: kb(epsilon) for multiple obstacle heights ===
leading a0\* = 0.391826552031
fixed plotted BEM order M = 24
jobs = 99
CPU policy: detected=10, fraction=50%, workers=5
subcritical policy: numerical empty-contour ladder, no sigma->0 fallback
supercritical policy: numerical cutoff-margin ladder + fixed-M=24 SVD localization
[ 1/99] a=0.350, eps=0.045: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 2/99] a=0.350, eps=0.071: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 3/99] a=0.350, eps=0.096: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 4/99] a=0.350, eps=0.020: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 5/99] a=0.350, eps=0.121: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 6/99] a=0.350, eps=0.146: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 7/99] a=0.350, eps=0.172: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 8/99] a=0.350, eps=0.197: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 9/99] a=0.350, eps=0.222: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 10/99] a=0.350, eps=0.247: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 11/99] a=0.350, eps=0.273: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 12/99] a=0.350, eps=0.298: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 13/99] a=0.350, eps=0.323: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 14/99] a=0.350, eps=0.374: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 15/99] a=0.350, eps=0.348: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 16/99] a=0.350, eps=0.399: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 17/99] a=0.350, eps=0.424: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 18/99] a=0.350, eps=0.449: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 19/99] a=0.350, eps=0.475: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 20/99] a=0.350, eps=0.500: ambiguous, rank=0, margin=1.0e-05, kb_M24=--
[ 21/99] a=0.500, eps=0.071: one, rank=1, margin=1.0e-05, kb_M24=1.570767064750
[ 22/99] a=0.500, eps=0.096: one, rank=1, margin=1.0e-05, kb_M24=1.570699661410
[ 23/99] a=0.500, eps=0.121: one, rank=1, margin=1.0e-05, kb_M24=1.570561520109
[ 24/99] a=0.500, eps=0.172: one, rank=1, margin=1.0e-05, kb_M24=1.569988940075
[ 25/99] a=0.500, eps=0.045: one, rank=1, margin=3.0e-06, kb_M24=1.570791297159
[ 26/99] a=0.500, eps=0.146: one, rank=1, margin=1.0e-05, kb_M24=1.570328765083
[ 27/99] a=0.500, eps=0.197: one, rank=1, margin=1.0e-05, kb_M24=1.569546856604
[ 28/99] a=0.500, eps=0.222: one, rank=1, margin=1.0e-05, kb_M24=1.569026815686
[ 29/99] a=0.500, eps=0.273: one, rank=1, margin=1.0e-05, kb_M24=1.567932243754
[ 30/99] a=0.500, eps=0.247: one, rank=1, margin=1.0e-05, kb_M24=1.568470663695
[ 31/99] a=0.500, eps=0.298: one, rank=1, margin=1.0e-05, kb_M24=1.567469718469
[ 32/99] a=0.500, eps=0.323: one, rank=1, margin=1.0e-05, kb_M24=1.567137023469
[ 33/99] a=0.500, eps=0.348: one, rank=1, margin=1.0e-05, kb_M24=1.566976282679
[ 34/99] a=0.500, eps=0.020: one, rank=1, margin=1.0e-07, kb_M24=1.570796134974
[ 35/99] a=0.500, eps=0.374: one, rank=1, margin=1.0e-05, kb_M24=1.567011871514
[ 36/99] a=0.500, eps=0.399: one, rank=1, margin=1.0e-05, kb_M24=1.567246896006
[ 37/99] a=0.500, eps=0.424: one, rank=1, margin=1.0e-05, kb_M24=1.567662377339
[ 38/99] a=0.500, eps=0.449: one, rank=1, margin=1.0e-05, kb_M24=1.568218892242
[ 39/99] a=0.550, eps=0.045: one, rank=1, margin=1.0e-05, kb_M24=1.570785593712
[ 40/99] a=0.500, eps=0.475: one, rank=1, margin=1.0e-05, kb_M24=1.568860166848
[ 41/99] a=0.550, eps=0.071: one, rank=1, margin=1.0e-05, kb_M24=1.570734605319
[ 42/99] a=0.550, eps=0.096: one, rank=1, margin=1.0e-05, kb_M24=1.570595103335
[ 43/99] a=0.550, eps=0.121: one, rank=1, margin=1.0e-05, kb_M24=1.570314144743
[ 44/99] a=0.550, eps=0.146: one, rank=1, margin=1.0e-05, kb_M24=1.569848632342
[ 45/99] a=0.550, eps=0.172: one, rank=1, margin=1.0e-05, kb_M24=1.569178745442
[ 46/99] a=0.550, eps=0.247: one, rank=1, margin=1.0e-05, kb_M24=1.566216210623
[ 47/99] a=0.550, eps=0.197: one, rank=1, margin=1.0e-05, kb_M24=1.568316185991
[ 48/99] a=0.550, eps=0.222: one, rank=1, margin=1.0e-05, kb_M24=1.567305114537
[ 49/99] a=0.550, eps=0.020: one, rank=1, margin=3.0e-07, kb_M24=1.570795913910
[ 50/99] a=0.550, eps=0.273: one, rank=1, margin=1.0e-05, kb_M24=1.565135961176
[ 51/99] a=0.550, eps=0.298: one, rank=1, margin=1.0e-05, kb_M24=1.564154286723
[ 52/99] a=0.550, eps=0.323: one, rank=1, margin=1.0e-05, kb_M24=1.563352466356
[ 53/99] a=0.550, eps=0.348: one, rank=1, margin=1.0e-05, kb_M24=1.562794143585
[ 54/99] a=0.550, eps=0.374: one, rank=1, margin=1.0e-05, kb_M24=1.562519908140
[ 55/99] a=0.550, eps=0.399: one, rank=1, margin=1.0e-05, kb_M24=1.562546133962
[ 56/99] a=0.600, eps=0.045: one, rank=1, margin=1.0e-05, kb_M24=1.570778007544
[ 57/99] a=0.550, eps=0.424: one, rank=1, margin=1.0e-05, kb_M24=1.562867371105
[ 58/99] a=0.600, eps=0.071: one, rank=1, margin=1.0e-05, kb_M24=1.570691878684
[ 59/99] a=0.550, eps=0.449: one, rank=1, margin=1.0e-05, kb_M24=1.563459921037
[ 60/99] a=0.600, eps=0.096: one, rank=1, margin=1.0e-05, kb_M24=1.570459311782
[ 61/99] a=0.600, eps=0.121: one, rank=1, margin=1.0e-05, kb_M24=1.569997894659
[ 62/99] a=0.600, eps=0.146: one, rank=1, margin=1.0e-05, kb_M24=1.569245326611
[ 63/99] a=0.600, eps=0.172: one, rank=1, margin=1.0e-05, kb_M24=1.568178787889
[ 64/99] a=0.600, eps=0.197: one, rank=1, margin=1.0e-05, kb_M24=1.566823836322
[ 65/99] a=0.600, eps=0.222: one, rank=1, margin=1.0e-05, kb_M24=1.565251160805
[ 66/99] a=0.600, eps=0.247: one, rank=1, margin=1.0e-05, kb_M24=1.563563777331
[ 67/99] a=0.600, eps=0.020: one, rank=1, margin=3.0e-07, kb_M24=1.570795617909
[ 68/99] a=0.600, eps=0.273: one, rank=1, margin=1.0e-05, kb_M24=1.561879243652
[ 69/99] a=0.600, eps=0.298: one, rank=1, margin=1.0e-05, kb_M24=1.560312145911
[ 70/99] a=0.600, eps=0.323: one, rank=1, margin=1.0e-05, kb_M24=1.558959752700
[ 71/99] a=0.600, eps=0.348: one, rank=1, margin=1.0e-05, kb_M24=1.557894768108
[ 72/99] a=0.600, eps=0.374: one, rank=1, margin=1.0e-05, kb_M24=1.557164948365
[ 73/99] a=0.600, eps=0.399: one, rank=1, margin=1.0e-05, kb_M24=1.556797462178
[ 74/99] a=0.650, eps=0.045: one, rank=1, margin=1.0e-05, kb_M24=1.570768980723
[ 75/99] a=0.650, eps=0.071: one, rank=1, margin=1.0e-05, kb_M24=1.570641584495
[ 76/99] a=0.650, eps=0.096: one, rank=1, margin=1.0e-05, kb_M24=1.570301768038
[ 77/99] a=0.650, eps=0.121: one, rank=1, margin=1.0e-05, kb_M24=1.569637349620
[ 78/99] a=0.650, eps=0.146: one, rank=1, margin=1.0e-05, kb_M24=1.568570933050
[ 79/99] a=0.650, eps=0.172: one, rank=1, margin=1.0e-05, kb_M24=1.567084101907
[ 80/99] a=0.650, eps=0.197: one, rank=1, margin=1.0e-05, kb_M24=1.565223765886
[ 81/99] a=0.650, eps=0.020: one, rank=1, margin=1.0e-06, kb_M24=1.570795263377
[ 82/99] a=0.650, eps=0.222: one, rank=1, margin=1.0e-05, kb_M24=1.563090817053
[ 83/99] a=0.650, eps=0.247: one, rank=1, margin=1.0e-05, kb_M24=1.560817404336
[ 84/99] a=0.650, eps=0.273: one, rank=1, margin=1.0e-05, kb_M24=1.558541179902
[ 85/99] a=0.650, eps=0.298: one, rank=1, margin=1.0e-05, kb_M24=1.556385088300
[ 86/99] a=0.650, eps=0.323: one, rank=1, margin=1.0e-05, kb_M24=1.554447236219
[ 87/99] a=0.650, eps=0.348: one, rank=1, margin=1.0e-05, kb_M24=1.552803076680
[ 88/99] a=0.700, eps=0.045: one, rank=1, margin=1.0e-05, kb_M24=1.570759126130
[ 89/99] a=0.700, eps=0.071: one, rank=1, margin=1.0e-05, kb_M24=1.570587429902
[ 90/99] a=0.700, eps=0.096: one, rank=1, margin=1.0e-05, kb_M24=1.570135321373
[ 91/99] a=0.700, eps=0.121: one, rank=1, margin=1.0e-05, kb_M24=1.569265250348
[ 92/99] a=0.700, eps=0.146: one, rank=1, margin=1.0e-05, kb_M24=1.567893274471
[ 93/99] a=0.700, eps=0.172: one, rank=1, margin=1.0e-05, kb_M24=1.566014731056
[ 94/99] a=0.700, eps=0.020: one, rank=1, margin=1.0e-06, kb_M24=1.570794873192
[ 95/99] a=0.700, eps=0.197: one, rank=1, margin=1.0e-05, kb_M24=1.563702470295
[ 96/99] a=0.700, eps=0.222: one, rank=1, margin=1.0e-05, kb_M24=1.561082115436
[ 97/99] a=0.700, eps=0.247: one, rank=1, margin=1.0e-05, kb_M24=1.558297667038
[ 98/99] a=0.700, eps=0.273: one, rank=1, margin=1.0e-05, kb_M24=1.555483548121
[ 99/99] a=0.700, eps=0.298: one, rank=1, margin=1.0e-05, kb_M24=1.552755693304

=== FINAL SUMMARY ===
Baseline: stable rank=1 + exactly one locally resolved mode: 10/10
Baseline leading asymptotic sigma within 5.0%: 4/10
Paper a-sweep M=24: one=79, zero=0, zero-supported=0, ambiguous=20, invalid=0, error=0, total=99
Subcritical 'zero-supported' means that Beyn rank may be contaminated near the first cutoff, but no resolved interior mode was found by the independent whole-band absolute-SVD scan.
Critical-height study: bracketed=0/4 for epsilon=(0.04, 0.06, 0.08, 0.1)
Publication claim audit:

-        supported: Uniqueness for sampled supercritical small-epsilon cases
-        supported: Local mode existence for sampled supercritical small-epsilon cases
- numerical-consistency: Leading asymptotic sigma = C(a) epsilon^2 + O(epsilon^3 |log epsilon|)
-          partial: Non-existence on sampled subcritical small-epsilon geometries
- not-established: Critical height a*(epsilon) = a0* + O(epsilon)
-       documented: Numerical implementation sensitivity at fixed M
- theoretical-only: Analyticity in epsilon and epsilon log epsilon

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-07_11-06-56_full_832cf4c8c71d
Publication audit: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-07_11-06-56_full_832cf4c8c71d/publication_validation_report.md
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations %
