(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations % python second_theorem_X.py
=== Theorem 2.3(iii) numerical validation: odd-sector Beyn + BEM + sigma-SVD ===
b = 1.0
a = 0.0 (required exactly by statement iii)
shape beta = 0.5
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1)b = 1.570796326795
sqrt(Lambda_2)b = 3.141592653590

=== SHAPE / SYMMETRY CHECK ===
reference area = 3.53414792455
min |r'(t)| = 5.000000e-01
X-even defect = 0.000e+00
Y-odd defect = 0.000e+00

=== INDEPENDENT DIPOLE-STRENGTH mu (MFS) ===
circle calibration: mu=1, nu=1.560e-14
N= 60: mu=1.07977772035, nu=-2.535e-16, boundary_res=7.806e-16
N= 80: mu=1.07977772066, nu=5.647e-15, boundary_res=1.229e-14
N= 120: mu=1.07977772066, nu=-2.189e-16, boundary_res=5.565e-15
N= 160: mu=1.07977772066, nu=1.281e-15, boundary_res=4.304e-15
selected mu = 1.07977772066
selected nu = 1.281e-15 (must vanish by x-axis symmetry)
C = pi^3 mu/b^3 = 33.4798867602

complex-kb odd-sector assembly: PASS
target real band ~ [1.570896326795, 3.141585517668]

=== EMBEDDED-MODE EPSILON SWEEP ===
epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
Beyn full M=24 -> odd reduced dimension=12; Nq base=(96, 192, 384), adaptive=(768, 1536)
fixed full BEM/SVD order M = 24

=== epsilon=0.02000000 ===
predicted kb = 3.141564109904
predicted sigma = 1.33919547e-02
predicted Lambda_2 cutoff gap in kb = 2.854e-05
effective Beyn upper margin = 7.136e-06

--- Theorem 2.3(iii) validation result ---
epsilon = 0.020000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 768
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.141564109904
kb BEM = 3.141564257466
sigma asym = 1.33919547e-02
sigma BEM = 1.33572940e-02
sigma/epsilon^2 = 3.33932351e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 1.10750428e+00
relative sigma error = 0.259%
relative singular value = 1.010e-09
minimum drop = 6.865e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 6.132e-10
one embedded mode supported = PASS
wall residual = 1.452e-11
centerline residual = 2.924e-11
odd-parity residual = 6.543e-11
first open-channel amplitude = 1.184e-11
off-grid BIE residual = 4.033e-03
decay left/right = 1.335736e-02/1.335739e-02
physical embedded-mode checks = PASS

=== epsilon=0.04000000 ===
predicted kb = 3.141135923496
predicted sigma = 5.35678188e-02
predicted Lambda_2 cutoff gap in kb = 4.567e-04
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.040000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.141135923496
kb BEM = 3.141150530224
sigma asym = 5.35678188e-02
sigma BEM = 5.27043410e-02
sigma/epsilon^2 = 3.29402131e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 4.19147599e+00
relative sigma error = 1.612%
relative singular value = 7.509e-09
minimum drop = 9.754e+07
cross-M change = N/A (single fixed M)
parity commutator residual = 3.105e-10
one embedded mode supported = PASS
wall residual = 1.689e-11
centerline residual = 3.068e-11
odd-parity residual = 7.891e-11
first open-channel amplitude = 2.310e-11
off-grid BIE residual = 1.011e-03
decay left/right = 5.270461e-02/5.270481e-02
physical embedded-mode checks = PASS

=== epsilon=0.06000000 ===
predicted kb = 3.139279774180
predicted sigma = 1.20527592e-01
predicted Lambda_2 cutoff gap in kb = 2.313e-03
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.060000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.139279774180
kb BEM = 3.139474736791
sigma asym = 1.20527592e-01
sigma BEM = 1.15337670e-01
sigma/epsilon^2 = 3.20382417e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 8.54031647e+00
relative sigma error = 4.306%
relative singular value = 5.883e-09
minimum drop = 1.268e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 2.138e-10
one embedded mode supported = PASS
wall residual = 2.496e-11
centerline residual = 5.986e-11
odd-parity residual = 1.100e-10
first open-channel amplitude = 1.438e-11
off-grid BIE residual = 4.504e-04
decay left/right = 1.153387e-01/1.153396e-01
physical embedded-mode checks = PASS

=== epsilon=0.08000000 ===
predicted kb = 3.134276985476
predicted sigma = 2.14271275e-01
predicted Lambda_2 cutoff gap in kb = 7.316e-03
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.080000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.134276985476
kb BEM = 3.135444038089
sigma asym = 2.14271275e-01
sigma BEM = 1.96456319e-01
sigma/epsilon^2 = 3.06962998e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 1.37761582e+01
relative sigma error = 8.314%
relative singular value = 2.274e-09
minimum drop = 3.375e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 1.644e-10
one embedded mode supported = PASS
wall residual = 2.068e-11
centerline residual = 2.468e-11
odd-parity residual = 7.632e-11
first open-channel amplitude = 1.137e-11
off-grid BIE residual = 2.542e-04
decay left/right = 1.964593e-01/1.964624e-01
physical embedded-mode checks = PASS

=== epsilon=0.10000000 ===
predicted kb = 3.123701989522
predicted sigma = 3.34798868e-01
predicted Lambda_2 cutoff gap in kb = 1.789e-02
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.100000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.123701989522
kb BEM = 3.128197565713
sigma asym = 3.34798868e-01
sigma BEM = 2.89800606e-01
sigma/epsilon^2 = 2.89800606e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 1.95424969e+01
relative sigma error = 13.440%
relative singular value = 9.932e-09
minimum drop = 8.033e+07
cross-M change = N/A (single fixed M)
parity commutator residual = 1.464e-10
one embedded mode supported = PASS
wall residual = 2.755e-11
centerline residual = 2.198e-11
odd-parity residual = 8.391e-11
first open-channel amplitude = 2.720e-11
off-grid BIE residual = 1.633e-04
decay left/right = 2.898069e-01/2.898151e-01
physical embedded-mode checks = PASS

=== epsilon=0.12000000 ===
predicted kb = 3.104379808087
predicted sigma = 4.82110369e-01
predicted Lambda_2 cutoff gap in kb = 3.721e-02
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.120000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.104379808087
kb BEM = 3.117457529950
sigma asym = 4.82110369e-01
sigma BEM = 3.88668175e-01
sigma/epsilon^2 = 2.69908455e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 2.55040674e+01
relative sigma error = 19.382%
relative singular value = 2.503e-09
minimum drop = 3.338e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 1.310e-10
one embedded mode supported = PASS
wall residual = 2.288e-11
centerline residual = 3.426e-11
odd-parity residual = 8.519e-11
first open-channel amplitude = 2.247e-11
off-grid BIE residual = 1.139e-04
decay left/right = 3.886838e-01/3.887045e-01
physical embedded-mode checks = PASS

=== epsilon=0.14000000 ===
predicted kb = 3.072295294194
predicted sigma = 6.56205780e-01
predicted Lambda_2 cutoff gap in kb = 6.930e-02
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.140000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.072295294194
kb BEM = 3.103629407143
sigma asym = 6.56205780e-01
sigma BEM = 4.86917759e-01
sigma/epsilon^2 = 2.48427428e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 3.13786084e+01
relative sigma error = 25.798%
relative singular value = 2.067e-09
minimum drop = 4.244e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 1.279e-10
one embedded mode supported = PASS
wall residual = 3.714e-11
centerline residual = 2.258e-11
odd-parity residual = 3.349e-11
first open-channel amplitude = 2.516e-11
off-grid BIE residual = 8.414e-05
decay left/right = 4.869539e-01/4.869995e-01
physical embedded-mode checks = PASS

=== epsilon=0.16000000 ===
predicted kb = 3.022417828599
predicted sigma = 8.57085101e-01
predicted Lambda_2 cutoff gap in kb = 1.192e-01
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.160000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 3.022417828599
kb BEM = 3.087656600254
sigma asym = 8.57085101e-01
sigma BEM = 5.79638784e-01
sigma/epsilon^2 = 2.26421400e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 3.69620225e+01
relative sigma error = 32.371%
relative singular value = 1.798e-09
minimum drop = 5.040e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 1.298e-10
one embedded mode supported = PASS
wall residual = 4.767e-11
centerline residual = 2.858e-11
odd-parity residual = 2.817e-11
first open-channel amplitude = 2.010e-11
off-grid BIE residual = 6.479e-05
decay left/right = 5.797061e-01/5.797961e-01
physical embedded-mode checks = PASS

=== epsilon=0.18000000 ===
predicted kb = 2.948376749912
predicted sigma = 1.08474833e+00
predicted Lambda_2 cutoff gap in kb = 1.932e-01
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.180000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 2.948376749912
kb BEM = 3.070755229736
sigma asym = 1.08474833e+00
sigma BEM = 6.63375248e-01
sigma/epsilon^2 = 2.04745447e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 4.21343402e+01
relative sigma error = 38.845%
relative singular value = 8.330e-09
minimum drop = 1.076e+08
cross-M change = N/A (single fixed M)
parity commutator residual = 1.424e-10
one embedded mode supported = PASS
wall residual = 4.001e-11
centerline residual = 2.181e-11
odd-parity residual = 3.907e-11
first open-channel amplitude = 5.159e-11
off-grid BIE residual = 5.304e-05
decay left/right = 6.634979e-01/6.636957e-01
physical embedded-mode checks = PASS

=== epsilon=0.20000000 ===
predicted kb = 2.841858527994
predicted sigma = 1.33919547e+00
predicted Lambda_2 cutoff gap in kb = 2.997e-01
effective Beyn upper margin = 1.000e-05

--- Theorem 2.3(iii) validation result ---
epsilon = 0.200000
odd-sector Beyn rank = 1
Beyn rank stable = YES
final Beyn Nq = 384
tighter-cutoff rank check = True
resolved odd modes = 1
kb asymptotic = 2.841858527994
kb BEM = 3.054163381761
sigma asym = 1.33919547e+00
sigma BEM = 7.35996222e-01
sigma/epsilon^2 = 1.83999055e+01
C=pi^3 mu/b^3 = 3.34798868e+01
scaled asymptotic remainder = 4.68485957e+01
relative sigma error = 45.042%
relative singular value = 1.546e-08
minimum drop = 5.717e+07
cross-M change = N/A (single fixed M)
parity commutator residual = 1.475e-10
one embedded mode supported = PASS
wall residual = 2.775e-11
centerline residual = 2.461e-11
odd-parity residual = 7.163e-11
first open-channel amplitude = 6.579e-11
off-grid BIE residual = 6.500e-05
decay left/right = 7.362236e-01/7.365197e-01
physical embedded-mode checks = PASS

=== FINAL SUMMARY ===
stable odd-sector rank=1 + exactly one resolved mode: 10/10
full embedded-field physical diagnostics: 10/10
Publication claim audit:

-        supported: Theorem 2.3(iii) symmetry assumptions (a=0, X even, Y odd, nu=0)
-        supported: One odd embedded trapped mode in [Lambda_1,Lambda_2) for sampled small epsilon
-        supported: sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)
-       documented: Numerical implementation sensitivity at fixed M
- theoretical-only: Analyticity in epsilon and epsilon log epsilon

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/theorem_2_3_iii_embedded_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations %
