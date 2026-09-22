(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations % python second_theorem_Y_1.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iv): y-symmetric BIC ===
b = 1.0
beta = 0.1
leading a1 = -beta/12 = -8.333333333333e-03
mu = 1.007140755654e+00
nu = 1.042239055775e-17
X-odd residual = 0.000e+00
Y-even residual = 0.000e+00
Theorem shape symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
refinement M = (8, 16, 24, 32, 40)

=== epsilon=0.020 ===
a = epsilon\*a1 = -1.666666666667e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.141567821035
predicted sigma = 1.24910740e-02
predicted kb gap to sqrt(Lambda_2) = 2.483e-05
expected-BIC mesh refinement:
M= 8: kb=3.141567894364, sigma_BEM=1.24726177e-02, sv_min=9.859e-05, drop=5.59e+03, mesh_change=--, y-axis parity=even(1.203e-10), Rprop=9.067e-06, interior=yes
M=16: kb=3.141567894606, sigma_BEM=1.24725567e-02, sv_min=1.031e-04, drop=5.35e+03, mesh_change=0.000%, y-axis parity=even(8.858e-11), Rprop=9.067e-06, interior=yes
M=24: kb=3.141567894628, sigma_BEM=1.24725511e-02, sv_min=1.035e-04, drop=5.33e+03, mesh_change=0.000%, y-axis parity=even(6.916e-11), Rprop=9.067e-06, interior=yes
M=32: kb=3.141567894634, sigma_BEM=1.24725497e-02, sv_min=1.036e-04, drop=5.32e+03, mesh_change=0.000%, y-axis parity=even(6.832e-11), Rprop=9.067e-06, interior=yes
M=40: kb=3.141567894635, sigma_BEM=1.24725493e-02, sv_min=1.036e-04, drop=5.32e+03, mesh_change=0.000%, y-axis parity=even(4.937e-11), Rprop=9.067e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592491047, sv_min=9.976e-01
kb(M=32)=3.141592491047, sv_min=9.976e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.532e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592613929, sv_min=9.976e-01
kb(M=32)=3.141592613929, sv_min=9.976e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(3.202e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592647935, sv_min=9.976e-01
kb(M=32)=3.141592647935, sv_min=9.976e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.918e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592652046, sv_min=9.976e-01
kb(M=32)=3.141592652046, sv_min=9.976e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.683e-06), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592653169, sv_min=9.976e-01
kb(M=32)=3.141592653169, sv_min=9.976e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(3.547e-06), Rprop=9.999e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.66666667e-04]
expected spectral minimum resolved: FAIL
non-radiating BIC diagnostics: FAIL
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): FAIL
kb asymptotic = 3.141567821035
kb BEM = 3.141567894635
sigma asym = 1.24910740e-02
sigma BEM = 1.24725493e-02
final sv_min(A) = 1.036e-04
final minimum drop = 5.32e+03
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 4.937e-11
max first propagating-mode fraction = 9.067e-06
relative sigma error = 0.148%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.040 ===
a = epsilon\*a1 = -3.333333333333e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.141195309150
predicted sigma = 4.99642959e-02
predicted kb gap to sqrt(Lambda_2) = 3.973e-04
expected-BIC mesh refinement:
M= 8: kb=3.141205381576, sigma_BEM=4.93269890e-02, sv_min=6.305e-06, drop=8.87e+04, mesh_change=--, y-axis parity=even(2.248e-10), Rprop=7.686e-05, interior=yes
M=16: kb=3.141205408833, sigma_BEM=4.93252532e-02, sv_min=4.132e-06, drop=1.35e+05, mesh_change=0.004%, y-axis parity=even(3.805e-11), Rprop=7.685e-05, interior=yes
M=24: kb=3.141205411358, sigma_BEM=4.93250925e-02, sv_min=3.926e-06, drop=1.42e+05, mesh_change=0.000%, y-axis parity=even(4.578e-11), Rprop=7.685e-05, interior=yes
M=32: kb=3.141205411962, sigma_BEM=4.93250540e-02, sv_min=3.877e-06, drop=1.44e+05, mesh_change=0.000%, y-axis parity=even(4.997e-11), Rprop=7.685e-05, interior=yes
M=40: kb=3.141205412176, sigma_BEM=4.93250403e-02, sv_min=3.860e-06, drop=1.45e+05, mesh_change=0.000%, y-axis parity=even(4.192e-11), Rprop=7.685e-05, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592200996, sv_min=9.888e-01
kb(M=32)=3.141592374193, sv_min=9.888e-01, drop=1.00e+00, shift=1.732e-07, y-axis parity=odd(1.159e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592567655, sv_min=9.888e-01
kb(M=32)=3.141592567655, sv_min=9.888e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.296e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592632870, sv_min=9.888e-01
kb(M=32)=3.141592632870, sv_min=9.888e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.394e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592650636, sv_min=9.888e-01
kb(M=32)=3.141592650636, sv_min=9.888e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.152e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592652784, sv_min=9.888e-01
kb(M=32)=3.141592652784, sv_min=9.888e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.409e-06), Rprop=9.999e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-3.33333333e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.141195309150
kb BEM = 3.141205412176
sigma asym = 4.99642959e-02
sigma BEM = 4.93250403e-02
final sv_min(A) = 3.860e-06
final minimum drop = 1.45e+05
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 4.192e-11
max first propagating-mode fraction = 7.685e-05
relative sigma error = 1.279%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.060 ===
a = epsilon\*a1 = -5.000000000000e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.139580580244
predicted sigma = 1.12419666e-01
predicted kb gap to sqrt(Lambda_2) = 2.012e-03
expected-BIC mesh refinement:
M= 8: kb=3.139721672997, sigma_BEM=1.08407644e-01, sv_min=5.106e-06, drop=1.13e+05, mesh_change=--, y-axis parity=even(4.671e-11), Rprop=2.837e-04, interior=yes
M=16: kb=3.139721979273, sigma_BEM=1.08398773e-01, sv_min=4.303e-06, drop=1.34e+05, mesh_change=0.008%, y-axis parity=even(4.186e-11), Rprop=2.836e-04, interior=yes
M=24: kb=3.139721998349, sigma_BEM=1.08398220e-01, sv_min=1.983e-06, drop=2.91e+05, mesh_change=0.001%, y-axis parity=even(4.745e-11), Rprop=2.836e-04, interior=yes
M=32: kb=3.139722003141, sigma_BEM=1.08398082e-01, sv_min=1.609e-06, drop=3.59e+05, mesh_change=0.000%, y-axis parity=even(3.292e-11), Rprop=2.836e-04, interior=yes
M=40: kb=3.139722004876, sigma_BEM=1.08398031e-01, sv_min=1.511e-06, drop=3.82e+05, mesh_change=0.000%, y-axis parity=even(4.823e-11), Rprop=2.836e-04, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592247660, sv_min=9.795e-01
kb(M=32)=3.141592247660, sv_min=9.795e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.849e-08), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592567655, sv_min=9.795e-01
kb(M=32)=3.141592567655, sv_min=9.795e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.090e-07), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592613929, sv_min=9.795e-01
kb(M=32)=3.141592613929, sv_min=9.795e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.930e-07), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592653169, sv_min=9.795e-01
kb(M=32)=3.141592653169, sv_min=9.795e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.935e-06), Rprop=9.998e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-5.00000000e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.139580580244
kb BEM = 3.139722004876
sigma asym = 1.12419666e-01
sigma BEM = 1.08398031e-01
final sv_min(A) = 1.511e-06
final minimum drop = 3.82e+05
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 4.823e-11
max first propagating-mode fraction = 2.836e-04
relative sigma error = 3.577%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.080 ===
a = epsilon\*a1 = -6.666666666667e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.135229099648
predicted sigma = 1.99857184e-01
predicted kb gap to sqrt(Lambda_2) = 6.364e-03
expected-BIC mesh refinement:
M= 8: kb=3.136095591085, sigma_BEM=1.85765564e-01, sv_min=4.438e-06, drop=1.37e+05, mesh_change=--, y-axis parity=even(2.703e-11), Rprop=7.551e-04, interior=yes
M=16: kb=3.136097194410, sigma_BEM=1.85738494e-01, sv_min=4.284e-06, drop=1.42e+05, mesh_change=0.015%, y-axis parity=even(2.441e-11), Rprop=7.550e-04, interior=yes
M=24: kb=3.136097348318, sigma_BEM=1.85735896e-01, sv_min=4.528e-06, drop=1.34e+05, mesh_change=0.001%, y-axis parity=even(3.429e-11), Rprop=7.550e-04, interior=yes
M=32: kb=3.136097321784, sigma_BEM=1.85736344e-01, sv_min=5.712e-06, drop=1.06e+05, mesh_change=0.000%, y-axis parity=even(1.895e-11), Rprop=7.549e-04, interior=yes
M=40: kb=3.136097388274, sigma_BEM=1.85735221e-01, sv_min=4.286e-06, drop=1.42e+05, mesh_change=0.001%, y-axis parity=even(2.684e-11), Rprop=7.549e-04, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592491047, sv_min=9.721e-01
kb(M=32)=3.141592491047, sv_min=9.721e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(9.377e-08), Rprop=9.996e-01, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592580829, sv_min=9.721e-01
kb(M=32)=3.141592580829, sv_min=9.721e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.998e-07), Rprop=9.996e-01, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592642766, sv_min=9.721e-01
kb(M=32)=3.141592642766, sv_min=9.721e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.450e-07), Rprop=9.996e-01, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592652046, sv_min=9.721e-01
kb(M=32)=3.141592652046, sv_min=9.721e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.395e-06), Rprop=9.996e-01, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592653169, sv_min=9.721e-01
kb(M=32)=3.141592653169, sv_min=9.721e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.790e-06), Rprop=9.996e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-6.66666667e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.135229099648
kb BEM = 3.136097388274
sigma asym = 1.99857184e-01
sigma BEM = 1.85735221e-01
final sv_min(A) = 4.286e-06
final minimum drop = 1.42e+05
final mesh change in sigma = 0.001%
y-axis parity diagnostic = even residual 2.684e-11
max first propagating-mode fraction = 7.549e-04
relative sigma error = 7.066%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.100 ===
a = epsilon\*a1 = -8.333333333333e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.126033840269
predicted sigma = 3.12276849e-01
predicted kb gap to sqrt(Lambda_2) = 1.556e-02
expected-BIC mesh refinement:
M= 8: kb=3.129437642286, sigma_BEM=2.76087747e-01, sv_min=9.862e-06, drop=6.58e+04, mesh_change=--, y-axis parity=even(5.908e-11), Rprop=1.686e-03, interior=yes
M=16: kb=3.129442891727, sigma_BEM=2.76028239e-01, sv_min=9.919e-06, drop=6.55e+04, mesh_change=0.022%, y-axis parity=even(7.231e-11), Rprop=1.686e-03, interior=yes
M=24: kb=3.129443383112, sigma_BEM=2.76022668e-01, sv_min=9.928e-06, drop=6.54e+04, mesh_change=0.002%, y-axis parity=even(4.588e-11), Rprop=1.686e-03, interior=yes
M=32: kb=3.129443502183, sigma_BEM=2.76021318e-01, sv_min=9.931e-06, drop=6.54e+04, mesh_change=0.000%, y-axis parity=even(2.986e-11), Rprop=1.686e-03, interior=yes
M=40: kb=3.129443545738, sigma_BEM=2.76020824e-01, sv_min=9.932e-06, drop=6.54e+04, mesh_change=0.000%, y-axis parity=even(1.824e-11), Rprop=1.686e-03, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 4 total; 3 outside expected-BIC window
checking additional minimum 1/3
kb(M=16)=3.141592632870, sv_min=9.683e-01
kb(M=32)=3.141592632870, sv_min=9.683e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(5.713e-07), Rprop=9.991e-01, persistent BIC=no
checking additional minimum 2/3
kb(M=16)=3.141592650636, sv_min=9.683e-01
kb(M=32)=3.141592650636, sv_min=9.683e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(9.126e-07), Rprop=9.991e-01, persistent BIC=no
checking additional minimum 3/3
kb(M=16)=3.141592653169, sv_min=9.683e-01
kb(M=32)=3.141592653169, sv_min=9.683e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.796e-06), Rprop=9.991e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-8.33333333e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.126033840269
kb BEM = 3.129443545738
sigma asym = 3.12276849e-01
sigma BEM = 2.76020824e-01
final sv_min(A) = 9.932e-06
final minimum drop = 6.54e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 1.824e-11
max first propagating-mode fraction = 1.686e-03
relative sigma error = 11.610%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.120 ===
a = epsilon\*a1 = -1.000000000000e-03
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.109243236093
predicted sigma = 4.49678663e-01
predicted kb gap to sqrt(Lambda_2) = 3.235e-02
expected-BIC mesh refinement:
M= 8: kb=3.119314312851, sigma_BEM=3.73473719e-01, sv_min=2.021e-05, drop=3.48e+04, mesh_change=--, y-axis parity=even(6.232e-11), Rprop=3.358e-03, interior=yes
M=16: kb=3.119327495313, sigma_BEM=3.73363600e-01, sv_min=2.043e-05, drop=3.44e+04, mesh_change=0.029%, y-axis parity=even(2.481e-11), Rprop=3.357e-03, interior=yes
M=24: kb=3.119328714200, sigma_BEM=3.73353417e-01, sv_min=2.045e-05, drop=3.44e+04, mesh_change=0.003%, y-axis parity=even(2.968e-11), Rprop=3.357e-03, interior=yes
M=32: kb=3.119328988951, sigma_BEM=3.73351121e-01, sv_min=2.045e-05, drop=3.44e+04, mesh_change=0.001%, y-axis parity=even(3.281e-11), Rprop=3.356e-03, interior=yes
M=40: kb=3.119329085143, sigma_BEM=3.73350318e-01, sv_min=2.046e-05, drop=3.44e+04, mesh_change=0.000%, y-axis parity=even(2.813e-11), Rprop=3.356e-03, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 4 total; 3 outside expected-BIC window
checking additional minimum 1/3
kb(M=16)=3.141592642766, sv_min=9.694e-01
kb(M=32)=3.141592642766, sv_min=9.695e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.194e-07), Rprop=9.959e-01, persistent BIC=no
checking additional minimum 2/3
kb(M=16)=3.141592652046, sv_min=9.694e-01
kb(M=32)=3.141592652046, sv_min=9.695e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.033e-06), Rprop=9.959e-01, persistent BIC=no
checking additional minimum 3/3
kb(M=16)=3.141592653169, sv_min=9.694e-01
kb(M=32)=3.141592653169, sv_min=9.695e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.409e-06), Rprop=9.959e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.00000000e-03]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.109243236093
kb BEM = 3.119329085143
sigma asym = 4.49678663e-01
sigma BEM = 3.73350318e-01
final sv_min(A) = 2.046e-05
final minimum drop = 3.44e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 2.813e-11
max first propagating-mode fraction = 3.356e-03
relative sigma error = 16.974%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.140 ===
a = epsilon\*a1 = -1.166666666667e-03
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.081393149977
predicted sigma = 6.12062625e-01
predicted kb gap to sqrt(Lambda_2) = 6.020e-02
expected-BIC mesh refinement:
M= 8: kb=3.105891270462, sigma_BEM=4.72275150e-01, sv_min=3.754e-05, drop=2.05e+04, mesh_change=--, y-axis parity=even(2.222e-11), Rprop=6.140e-03, interior=yes
M=16: kb=3.105918418557, sigma_BEM=4.72096577e-01, sv_min=3.793e-05, drop=2.03e+04, mesh_change=0.038%, y-axis parity=even(3.487e-11), Rprop=6.136e-03, interior=yes
M=24: kb=3.105920965913, sigma_BEM=4.72079818e-01, sv_min=3.797e-05, drop=2.02e+04, mesh_change=0.004%, y-axis parity=even(4.126e-11), Rprop=6.135e-03, interior=yes
M=32: kb=3.105921568203, sigma_BEM=4.72075855e-01, sv_min=3.798e-05, drop=2.02e+04, mesh_change=0.001%, y-axis parity=even(3.844e-11), Rprop=6.135e-03, interior=yes
M=40: kb=3.105921783749, sigma_BEM=4.72074437e-01, sv_min=3.798e-05, drop=2.02e+04, mesh_change=0.000%, y-axis parity=even(2.629e-11), Rprop=6.135e-03, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 2 total; 1 outside expected-BIC window
checking additional minimum 1/1
kb(M=16)=3.141592652784, sv_min=9.723e-01
kb(M=32)=3.141592652784, sv_min=9.724e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.971e-06), Rprop=8.254e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.16666667e-03]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.081393149977
kb BEM = 3.105921783749
sigma asym = 6.12062625e-01
sigma BEM = 4.72074437e-01
final sv_min(A) = 3.798e-05
final minimum drop = 2.02e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 2.629e-11
max first propagating-mode fraction = 6.135e-03
relative sigma error = 22.872%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.160 ===
a = epsilon\*a1 = -1.333333333333e-03
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.038176772372
predicted sigma = 7.99428734e-01
predicted kb gap to sqrt(Lambda_2) = 1.034e-01
expected-BIC mesh refinement:
M= 8: kb=3.089866098075, sigma_BEM=5.67742809e-01, sv_min=6.486e-05, drop=1.11e+04, mesh_change=--, y-axis parity=even(2.797e-11), Rprop=1.046e-02, interior=yes
M=16: kb=3.089914481372, sigma_BEM=5.67479426e-01, sv_min=6.552e-05, drop=1.09e+04, mesh_change=0.046%, y-axis parity=even(1.977e-11), Rprop=1.045e-02, interior=yes
M=24: kb=3.089919034234, sigma_BEM=5.67454635e-01, sv_min=6.558e-05, drop=1.09e+04, mesh_change=0.004%, y-axis parity=even(3.187e-11), Rprop=1.045e-02, interior=yes
M=32: kb=3.089920123678, sigma_BEM=5.67448703e-01, sv_min=6.559e-05, drop=1.09e+04, mesh_change=0.001%, y-axis parity=even(3.197e-11), Rprop=1.045e-02, interior=yes
M=40: kb=3.089920509178, sigma_BEM=5.67446604e-01, sv_min=6.560e-05, drop=1.09e+04, mesh_change=0.000%, y-axis parity=even(1.935e-11), Rprop=1.045e-02, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 3 total; 2 outside expected-BIC window
checking additional minimum 1/2
kb(M=16)=3.141592652046, sv_min=9.636e-01
kb(M=32)=3.141592652046, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(3.228e-06), Rprop=2.993e-01, persistent BIC=no
checking additional minimum 2/2
kb(M=16)=3.141592653169, sv_min=9.636e-01
kb(M=32)=3.141592653169, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.190e-05), Rprop=2.993e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.33333333e-03]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: FAIL
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): FAIL
kb asymptotic = 3.038176772372
kb BEM = 3.089920509178
sigma asym = 7.99428734e-01
sigma BEM = 5.67446604e-01
final sv_min(A) = 6.560e-05
final minimum drop = 1.09e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 1.935e-11
max first propagating-mode fraction = 1.045e-02
relative sigma error = 29.018%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.180 ===
a = epsilon\*a1 = -1.500000000000e-03
placement condition: PASS (residual=0.000e+00)
predicted kb = 2.974207746672
predicted sigma = 1.01177699e+00
predicted kb gap to sqrt(Lambda_2) = 1.674e-01
expected-BIC mesh refinement:
M= 8: kb=3.072267069245, sigma_BEM=6.56337913e-01, sv_min=1.063e-04, drop=5.47e+03, mesh_change=--, y-axis parity=even(7.510e-11), Rprop=1.676e-02, interior=yes
M=16: kb=3.072344097383, sigma_BEM=6.55977247e-01, sv_min=1.074e-04, drop=5.41e+03, mesh_change=0.055%, y-axis parity=even(1.641e-11), Rprop=1.674e-02, interior=yes
M=24: kb=3.072351333823, sigma_BEM=6.55943353e-01, sv_min=1.074e-04, drop=5.40e+03, mesh_change=0.005%, y-axis parity=even(3.293e-11), Rprop=1.674e-02, interior=yes
M=32: kb=3.072353072460, sigma_BEM=6.55935210e-01, sv_min=1.075e-04, drop=5.40e+03, mesh_change=0.001%, y-axis parity=even(2.746e-11), Rprop=1.673e-02, interior=yes
M=40: kb=3.072353687409, sigma_BEM=6.55932329e-01, sv_min=1.075e-04, drop=5.40e+03, mesh_change=0.000%, y-axis parity=even(3.650e-11), Rprop=1.673e-02, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 3 total; 2 outside expected-BIC window
checking additional minimum 1/2
kb(M=16)=3.141592652046, sv_min=9.507e-01
kb(M=32)=3.141592652046, sv_min=9.508e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.340e-06), Rprop=1.544e-01, persistent BIC=no
checking additional minimum 2/2
kb(M=16)=3.141592653169, sv_min=9.507e-01
kb(M=32)=3.141592653169, sv_min=9.508e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(5.135e-06), Rprop=1.543e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.50000000e-03]
expected spectral minimum resolved: FAIL
non-radiating BIC diagnostics: FAIL
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): FAIL
kb asymptotic = 2.974207746672
kb BEM = 3.072353687409
sigma asym = 1.01177699e+00
sigma BEM = 6.55932329e-01
final sv_min(A) = 1.075e-04
final minimum drop = 5.40e+03
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 3.650e-11
max first propagating-mode fraction = 1.673e-02
relative sigma error = 35.170%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.200 ===
a = epsilon\*a1 = -1.666666666667e-03
placement condition: PASS (residual=0.000e+00)
predicted kb = 2.882591735187
predicted sigma = 1.24910740e+00
predicted kb gap to sqrt(Lambda_2) = 2.590e-01
expected-BIC mesh refinement:
M= 8: kb=3.054228063347, sigma_BEM=7.35727761e-01, sv_min=1.676e-04, drop=2.58e+03, mesh_change=--, y-axis parity=even(4.704e-11), Rprop=2.543e-02, interior=yes
M=16: kb=3.054340428418, sigma_BEM=7.35261143e-01, sv_min=1.691e-04, drop=2.55e+03, mesh_change=0.063%, y-axis parity=even(5.845e-11), Rprop=2.538e-02, interior=yes
M=24: kb=3.054351032009, sigma_BEM=7.35217093e-01, sv_min=1.692e-04, drop=2.55e+03, mesh_change=0.006%, y-axis parity=even(3.639e-11), Rprop=2.537e-02, interior=yes
M=32: kb=3.054353572603, sigma_BEM=7.35206539e-01, sv_min=1.692e-04, drop=2.55e+03, mesh_change=0.001%, y-axis parity=even(4.321e-11), Rprop=2.537e-02, interior=yes
M=40: kb=3.054354487675, sigma_BEM=7.35202737e-01, sv_min=1.693e-04, drop=2.55e+03, mesh_change=0.001%, y-axis parity=even(5.600e-11), Rprop=2.537e-02, interior=yes
whole-band uniqueness screen: 78 points (M=16)
16/78
32/78
48/78
64/78
78/78
sampled local minima: 1 total; 0 outside expected-BIC window

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-1.66666667e-03]
expected spectral minimum resolved: FAIL
non-radiating BIC diagnostics: FAIL
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): FAIL
kb asymptotic = 2.882591735187
kb BEM = 3.054354487675
sigma asym = 1.24910740e+00
sigma BEM = 7.35202737e-01
final sv_min(A) = 1.693e-04
final minimum drop = 2.55e+03
final mesh change in sigma = 0.001%
y-axis parity diagnostic = even residual 5.600e-11
max first propagating-mode fraction = 2.537e-02
relative sigma error = 41.142%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 6/10
Leading asymptotic sigma within 5.0%: 3/10
epsilon=0.020: a=-1.666667e-04, BIC=FAIL, asymptotic=PASS, sv_min=1.036e-04, mesh_change=0.000%, yparity=even(4.937e-11), Rprop=9.067e-06, sigma_error=0.148%
epsilon=0.040: a=-3.333333e-04, BIC=PASS, asymptotic=PASS, sv_min=3.860e-06, mesh_change=0.000%, yparity=even(4.192e-11), Rprop=7.685e-05, sigma_error=1.279%
epsilon=0.060: a=-5.000000e-04, BIC=PASS, asymptotic=PASS, sv_min=1.511e-06, mesh_change=0.000%, yparity=even(4.823e-11), Rprop=2.836e-04, sigma_error=3.577%
epsilon=0.080: a=-6.666667e-04, BIC=PASS, asymptotic=FAIL, sv_min=4.286e-06, mesh_change=0.001%, yparity=even(2.684e-11), Rprop=7.549e-04, sigma_error=7.066%
epsilon=0.100: a=-8.333333e-04, BIC=PASS, asymptotic=FAIL, sv_min=9.932e-06, mesh_change=0.000%, yparity=even(1.824e-11), Rprop=1.686e-03, sigma_error=11.610%
epsilon=0.120: a=-1.000000e-03, BIC=PASS, asymptotic=FAIL, sv_min=2.046e-05, mesh_change=0.000%, yparity=even(2.813e-11), Rprop=3.356e-03, sigma_error=16.974%
epsilon=0.140: a=-1.166667e-03, BIC=PASS, asymptotic=FAIL, sv_min=3.798e-05, mesh_change=0.000%, yparity=even(2.629e-11), Rprop=6.135e-03, sigma_error=22.872%
epsilon=0.160: a=-1.333333e-03, BIC=FAIL, asymptotic=FAIL, sv_min=6.560e-05, mesh_change=0.000%, yparity=even(1.935e-11), Rprop=1.045e-02, sigma_error=29.018%
epsilon=0.180: a=-1.500000e-03, BIC=FAIL, asymptotic=FAIL, sv_min=1.075e-04, mesh_change=0.000%, yparity=even(3.650e-11), Rprop=1.673e-02, sigma_error=35.170%
epsilon=0.200: a=-1.666667e-03, BIC=FAIL, asymptotic=FAIL, sv_min=1.693e-04, mesh_change=0.001%, yparity=even(5.600e-11), Rprop=2.537e-02, sigma_error=41.142%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/svd_theorems_validations/theorem_2_3_y_symmetry_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations %
