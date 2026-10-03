(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations % python second_theorem_Y.py
=== Theorem 2.3(iv) numerical validation: full BEM + complex-root tuning + fixed-M SVD ===
b = 1.0
shape beta = 0.5
fixed BEM order M = 24
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1)b = 1.570796326795
sqrt(Lambda_2)b = 3.141592653590

=== SHAPE / Y-AXIS SYMMETRY CHECK ===
reference area = 3.53429173529
min |r'(t)| = 5.000000e-01
X-odd defect = 0.000e+00
Y-even defect = 0.000e+00

=== INDEPENDENT MFS: mu, nu, Psi, a1 ===
circle calibration: mu=1, nu=-2.546e-16, a1=-3.701e-17, a1(5.47)=-3.701e-17
N= 60: mu=1.1769382794, nu=-3.407e-16, a1(2.7)=-0.229718010292, a1(5.47)=-0.229718010292, mismatch=0.000e+00, boundary_res=1.494e-10
N= 80: mu=1.17693827901, nu=4.362e-17, a1(2.7)=-0.229718010365, a1(5.47)=-0.229718010365, mismatch=0.000e+00, boundary_res=4.516e-14
N= 120: mu=1.17693827901, nu=-1.645e-16, a1(2.7)=-0.229718010365, a1(5.47)=-0.229718010365, mismatch=0.000e+00, boundary_res=3.142e-14
N= 160: mu=1.17693827901, nu=-2.858e-17, a1(2.7)=-0.229718010365, a1(5.47)=-0.229718010365, mismatch=0.000e+00, boundary_res=2.037e-14
selected mu = 1.17693827901
selected nu = -2.858e-17 (must vanish by y-axis symmetry)
selected a1 from formula (2.7) = -0.229718010365
selected a1 from formula (5.47) = -0.229718010365; relative mismatch=0.000e+00
paper Example 5.4 leading diagnostic -beta/12 = -0.0416666666667
C = pi^3 mu/b^3 = 36.4924739147

=== SMALL-BETA EXAMPLE 5.4 DIAGNOSTIC ===
beta=0.0200: a1(2.7)=-0.00999866671113, a1/beta=-0.49993334, -1/12=-0.083333333
beta=0.0500: a1(2.7)=-0.0249791710209, a1/beta=-0.49958342, -1/12=-0.083333333
beta=0.1000: a1(2.7)=-0.0498334740032, a1/beta=-0.49833474, -1/12=-0.083333333

complex-k full-matrix assembly: PASS

=== EMBEDDED-MODE / PLACEMENT SWEEP ===
epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)

=== epsilon=0.02000000 ===
predicted a leading = -4.594360207296e-03
predicted kb = 3.141558741928
predicted sigma = 1.45969896e-02

--- Theorem 2.3(iv) validation result ---
epsilon = 0.020000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -4.594360207296e-03
numerical tuned placement = -5.712549736707e-03
(a-eps*a1)/eps^2 = -2.79547382e+00
placement search interior = YES
kb asymptotic = 3.141558741928
kb BEM = 3.141558879982
sigma asym = 1.45969896e-02
sigma BEM = 1.45672473e-02
sigma/epsilon^2 = 3.64181182e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 0.204%
relative singular value at leading (a,sigma) prediction = 8.515e-04
relative singular value after tuning = 1.425e-08
spectral drop = 2.178e+06
complex root: success=True, residual=1.540e-08, |lambda0|=3.283e-08, seed=coarse-svd
embedded mode resolved = PASS
boundary x-even residual = 8.096e-11
field x-even residual = 7.766e-11
first open-channel amplitude = 1.422e-05
signed open-channel coefficients L/R = (-6.927e-08-2.721e-05j) / (-6.927e-08-2.721e-05j); phase consistency=9.490e-04
wall residual = 1.324e-11
off-grid BIE residual = 2.121e-03
decay left/right = 1.456734e-02/1.456734e-02
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.04000000 ===
predicted a leading = -9.188720414592e-03
predicted kb = 3.141050023069
predicted sigma = 5.83879583e-02
whole-band real-axis screen: 102 points, fixed M=24
20/102
40/102
60/102
80/102
100/102
102/102

--- Theorem 2.3(iv) validation result ---
epsilon = 0.040000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -9.188720414592e-03
numerical tuned placement = -9.181783149249e-03
(a-eps*a1)/eps^2 = 4.33579084e-03
placement search interior = YES
kb asymptotic = 3.141050023069
kb BEM = 3.141065516576
sigma asym = 5.83879583e-02
sigma BEM = 5.75484289e-02
sigma/epsilon^2 = 3.59677681e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 1.438%
relative singular value at leading (a,sigma) prediction = 6.078e-03
relative singular value after tuning = 1.529e-11
spectral drop = 2.063e+09
complex root: success=True, residual=1.634e-11, |lambda0|=3.481e-11, seed=theorem-leading
embedded mode resolved = PASS
boundary x-even residual = 9.798e-11
field x-even residual = 1.074e-10
first open-channel amplitude = 7.783e-07
signed open-channel coefficients L/R = (+7.917e-09+6.543e-07j) / (+7.922e-09+6.543e-07j); phase consistency=2.043e-02
wall residual = 1.682e-11
off-grid BIE residual = 5.294e-04
decay left/right = 5.754884e-02/5.754884e-02
physical embedded-mode checks = PASS
whole-band uniqueness screen = PASS (additional resolved BICs=0)

=== epsilon=0.06000000 ===
predicted a leading = -1.378308062189e-02
predicted kb = 3.138844621932
predicted sigma = 1.31372906e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.060000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -1.378308062189e-02
numerical tuned placement = -1.372994066381e-02
(a-eps*a1)/eps^2 = 1.47610995e-02
placement search interior = YES
kb asymptotic = 3.138844621932
kb BEM = 3.139059977187
sigma asym = 1.31372906e-01
sigma BEM = 1.26122404e-01
sigma/epsilon^2 = 3.50340010e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 3.997%
relative singular value at leading (a,sigma) prediction = 1.705e-02
relative singular value after tuning = 4.384e-10
spectral drop = 7.237e+07
complex root: success=True, residual=4.686e-10, |lambda0|=1.008e-09, seed=epsilon-continuation
embedded mode resolved = PASS
boundary x-even residual = 5.554e-11
field x-even residual = 4.217e-11
first open-channel amplitude = 1.100e-05
signed open-channel coefficients L/R = (+1.260e-07+5.072e-06j) / (+1.260e-07+5.072e-06j); phase consistency=2.829e-02
wall residual = 9.886e-12
off-grid BIE residual = 2.350e-04
decay left/right = 1.261240e-01/1.261240e-01
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.08000000 ===
predicted a leading = -1.837744082918e-02
predicted kb = 3.132899286981
predicted sigma = 2.33551833e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.080000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -1.837744082918e-02
numerical tuned placement = -1.861923911852e-02
(a-eps*a1)/eps^2 = -3.77809827e-02
placement search interior = YES
kb asymptotic = 3.132899286981
kb BEM = 3.134216474034
sigma asym = 2.33551833e-01
sigma BEM = 2.15154584e-01
sigma/epsilon^2 = 3.36179037e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 7.877%
relative singular value at leading (a,sigma) prediction = 3.436e-02
relative singular value after tuning = 7.742e-09
spectral drop = 4.168e+06
complex root: success=True, residual=8.281e-09, |lambda0|=1.800e-08, seed=epsilon-continuation
embedded mode resolved = PASS
boundary x-even residual = 4.435e-11
field x-even residual = 8.785e-11
first open-channel amplitude = 6.961e-05
signed open-channel coefficients L/R = (-7.335e-07-1.980e-05j) / (-7.335e-07-1.980e-05j); phase consistency=1.471e-02
wall residual = 2.490e-11
off-grid BIE residual = 1.330e-04
decay left/right = 2.151585e-01/2.151585e-01
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.10000000 ===
predicted a leading = -2.297180103648e-02
predicted kb = 3.120325998329
predicted sigma = 3.64924739e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.100000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -2.297180103648e-02
numerical tuned placement = -2.284861399208e-02
(a-eps*a1)/eps^2 = 1.23187044e-02
placement search interior = YES
kb asymptotic = 3.120325998329
kb BEM = 3.125470197211
sigma asym = 3.64924739e-01
sigma BEM = 3.17868601e-01
sigma/epsilon^2 = 3.17868601e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 12.895%
relative singular value at leading (a,sigma) prediction = 5.837e-02
relative singular value after tuning = 7.008e-09
spectral drop = 4.755e+06
complex root: success=True, residual=7.491e-09, |lambda0|=1.638e-08, seed=multistart(+0.5,-0.08)
embedded mode resolved = PASS
boundary x-even residual = 4.072e-11
field x-even residual = 7.911e-11
first open-channel amplitude = 1.167e-04
signed open-channel coefficients L/R = (+1.176e-06+1.977e-05j) / (+1.176e-06+1.977e-05j); phase consistency=1.439e-02
wall residual = 3.636e-11
off-grid BIE residual = 8.622e-05
decay left/right = 3.178806e-01/3.178806e-01
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.12000000 ===
predicted a leading = -2.756616124377e-02
predicted kb = 3.097331586028
predicted sigma = 5.25491624e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.120000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -2.756616124377e-02
numerical tuned placement = -3.115596961652e-02
(a-eps*a1)/eps^2 = -2.49292248e-01
placement search interior = YES
kb asymptotic = 3.097331586028
kb BEM = 3.112463376156
sigma asym = 5.25491624e-01
sigma BEM = 4.26820961e-01
sigma/epsilon^2 = 2.96403445e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 18.777%
relative singular value at leading (a,sigma) prediction = 8.930e-02
relative singular value after tuning = 4.582e-06
spectral drop = 7.516e+03
complex root: success=False, residual=4.927e-06, |lambda0|=1.087e-05, seed=epsilon-continuation
embedded mode resolved = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.14000000 ===
predicted a leading = -3.216052145107e-02
predicted kb = 3.059087818036
predicted sigma = 7.15252489e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.140000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -3.216052145107e-02
numerical tuned placement = -3.560088633788e-02
(a-eps*a1)/eps^2 = -1.75528821e-01
placement search interior = YES
kb asymptotic = 3.059087818036
kb BEM = 3.095607272346
sigma asym = 7.15252489e-01
sigma BEM = 5.35555801e-01
sigma/epsilon^2 = 2.73242756e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 25.124%
relative singular value at leading (a,sigma) prediction = 1.272e-01
relative singular value after tuning = 5.762e-06
spectral drop = 6.317e+03
complex root: success=False, residual=6.191e-06, |lambda0|=1.365e-05, seed=multistart(+0.5,+0.08)
embedded mode resolved = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.16000000 ===
predicted a leading = -3.675488165837e-02
predicted kb = 2.999476797963
predicted sigma = 9.34207332e-01

--- Theorem 2.3(iv) validation result ---
epsilon = 0.160000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -3.675488165837e-02
numerical tuned placement = -3.787111656377e-02
(a-eps*a1)/eps^2 = -4.36029260e-02
placement search interior = YES
kb asymptotic = 2.999476797963
kb BEM = 3.075957228029
sigma asym = 9.34207332e-01
sigma BEM = 6.38820423e-01
sigma/epsilon^2 = 2.49539228e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 31.619%
relative singular value at leading (a,sigma) prediction = 1.717e-01
relative singular value after tuning = 7.878e-07
spectral drop = 4.941e+04
complex root: success=False, residual=8.436e-07, |lambda0|=1.847e-06, seed=epsilon-continuation
embedded mode resolved = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.18000000 ===
predicted a leading = -4.134924186566e-02
predicted kb = 2.910607895991
predicted sigma = 1.18235615e+00

--- Theorem 2.3(iv) validation result ---
epsilon = 0.180000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -4.134924186566e-02
numerical tuned placement = -4.146588310573e-02
(a-eps*a1)/eps^2 = -3.60003827e-03
placement search interior = YES
kb asymptotic = 2.910607895991
kb BEM = 3.054939790192
sigma asym = 1.18235615e+00
sigma BEM = 7.32766866e-01
sigma/epsilon^2 = 2.26162613e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 38.025%
relative singular value at leading (a,sigma) prediction = 2.221e-01
relative singular value after tuning = 1.238e-08
spectral drop = 3.366e+06
complex root: success=True, residual=1.325e-08, |lambda0|=2.874e-08, seed=coarse-svd
embedded mode resolved = PASS
boundary x-even residual = 5.536e-11
field x-even residual = 3.814e-11
first open-channel amplitude = 5.880e-04
signed open-channel coefficients L/R = (-4.018e-06-2.264e-05j) / (-4.018e-06-2.264e-05j); phase consistency=8.793e-02
wall residual = 2.565e-11
off-grid BIE residual = 4.358e-05
decay left/right = 7.330171e-01/7.330171e-01
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== epsilon=0.20000000 ===
predicted a leading = -4.594360207296e-02
predicted kb = 2.781884856931
predicted sigma = 1.45969896e+00

--- Theorem 2.3(iv) validation result ---
epsilon = 0.200000
a1 formula (2.7) = -2.297180103648e-01
paper Example 5.4 -beta/12 = -4.166666666667e-02
leading placement eps*a1 = -4.594360207296e-02
numerical tuned placement = -4.589804444325e-02
(a-eps*a1)/eps^2 = 1.13894074e-03
placement search interior = YES
kb asymptotic = 2.781884856931
kb BEM = 3.033959047167
sigma asym = 1.45969896e+00
sigma BEM = 8.15289459e-01
sigma/epsilon^2 = 2.03822365e+01
C=pi^3 mu/b^3 = 3.64924739e+01
relative sigma error = 44.147%
relative singular value at leading (a,sigma) prediction = 2.767e-01
relative singular value after tuning = 1.921e-10
spectral drop = 2.324e+08
complex root: success=True, residual=2.055e-10, |lambda0|=4.407e-10, seed=coarse-svd
embedded mode resolved = PASS
boundary x-even residual = 4.661e-11
field x-even residual = 4.271e-11
first open-channel amplitude = 1.163e-04
signed open-channel coefficients L/R = (-3.365e-07-2.053e-06j) / (-3.365e-07-2.053e-06j); phase consistency=3.440e-01
wall residual = 1.388e-11
off-grid BIE residual = 5.030e-05
decay left/right = 8.157326e-01/8.157326e-01
physical embedded-mode checks = FAIL
whole-band uniqueness screen = not run for this epsilon

=== FINAL SUMMARY ===
resolved full-matrix embedded mode: 7/10
physical BIC diagnostics: 1/10
whole-band uniqueness scans passed: 1/1
Publication claim audit:

-        supported: Theorem 2.3(iv) symmetry assumptions: X odd, Y even, nu=0
-          partial: Placement law a = epsilon a1 + O(epsilon^2)
-          partial: One non-radiating embedded trapped mode for sampled tuned geometries
-        supported: sigma = epsilon^2 pi^3 mu/b^3 + O(epsilon^3 |log epsilon|)
-        supported: Uniqueness on the sampled real embedded band
- theoretical-only: Analyticity in epsilon and epsilon log epsilon

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/theorem_2_3_iv_embedded_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations %
