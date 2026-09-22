(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations % python second_theorem_X.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iii): x-symmetric BIC ===
b = 1.0
a = 0.0
beta = 0.1
mu = 1.002199318791e+00
nu = -1.687599000907e-17
X-even residual = 0.000e+00
Y-odd residual = 0.000e+00
Theorem symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590
epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
refinement M = (8, 16, 24, 32, 40)
safe Green-function epsilon upper bound = 0.340835817091
Running 10 validation cases with 4 worker processes.
0.13058424044339537
0.13058424044372274
0.13058424044339537
0.13058424044372274
0.13058424044339537
0.13058424044372274
0.13058424044339537
0.13058424044372274

=== epsilon=0.020 ===
predicted kb = 3.141568064115
predicted sigma = 1.24297877e-02
predicted kb gap to sqrt(Lambda_2) = 2.459e-05
expected-BIC mesh refinement:
M= 8: kb=3.141568137865, sigma_BEM=1.24111338e-02, sv_min=1.025e-04, drop=5.41e+03, mesh_change=--, odd_res=4.753e-10, Rprop=8.283e-08, interior=yes
M=16: kb=3.141568138106, sigma_BEM=1.24110728e-02, sv_min=1.071e-04, drop=5.18e+03, mesh_change=0.000%, odd_res=2.490e-10, Rprop=8.282e-08, interior=yes
M=24: kb=3.141568138128, sigma_BEM=1.24110672e-02, sv_min=1.075e-04, drop=5.16e+03, mesh_change=0.000%, odd_res=2.175e-10, Rprop=8.283e-08, interior=yes
M=32: kb=3.141568138133, sigma_BEM=1.24110658e-02, sv_min=1.076e-04, drop=5.16e+03, mesh_change=0.000%, odd_res=2.000e-10, Rprop=8.282e-08, interior=yes
M=40: kb=3.141568138135, sigma_BEM=1.24110653e-02, sv_min=1.076e-04, drop=5.15e+03, mesh_change=0.000%, odd_res=2.192e-10, Rprop=8.282e-08, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592042927, sv_min=9.918e-01
kb(M=32)=3.141592200996, sv_min=9.918e-01, drop=1.00e+00, shift=1.581e-07, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592580829, sv_min=9.918e-01
kb(M=32)=3.141592580829, sv_min=9.918e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592632870, sv_min=9.918e-01
kb(M=32)=3.141592632870, sv_min=9.918e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592653169, sv_min=9.918e-01
kb(M=32)=3.141592653169, sv_min=9.918e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: FAIL
non-radiating BIC diagnostics: FAIL
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): FAIL
kb asymptotic = 3.141568064115
kb BEM = 3.141568138135
sigma asym = 1.24297877e-02
sigma BEM = 1.24110653e-02
final sv_min(A) = 1.076e-04
final minimum drop = 5.15e+03
final mesh change in sigma = 0.000%
odd-parity residual = 2.192e-10
max first propagating-mode fraction = 8.282e-08
relative sigma error = 0.151%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.040 ===
predicted kb = 3.141199198891
predicted sigma = 4.97191510e-02
predicted kb gap to sqrt(Lambda_2) = 3.935e-04
expected-BIC mesh refinement:
M= 8: kb=3.141209213977, sigma_BEM=4.90823300e-02, sv_min=2.282e-06, drop=2.47e+05, mesh_change=--, odd_res=1.179e-10, Rprop=3.547e-07, interior=yes
M=16: kb=3.141209241104, sigma_BEM=4.90805939e-02, sv_min=5.417e-08, drop=1.04e+07, mesh_change=0.004%, odd_res=1.237e-10, Rprop=3.547e-07, interior=yes
M=24: kb=3.141209243604, sigma_BEM=4.90804338e-02, sv_min=2.734e-07, drop=2.06e+06, mesh_change=0.000%, odd_res=1.151e-10, Rprop=3.547e-07, interior=yes
M=32: kb=3.141209244203, sigma_BEM=4.90803955e-02, sv_min=3.261e-07, drop=1.73e+06, mesh_change=0.000%, odd_res=1.546e-10, Rprop=3.547e-07, interior=yes
M=40: kb=3.141209244415, sigma_BEM=4.90803819e-02, sv_min=3.447e-07, drop=1.63e+06, mesh_change=0.000%, odd_res=9.085e-11, Rprop=3.547e-07, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592247660, sv_min=9.830e-01
kb(M=32)=3.141592325862, sv_min=9.830e-01, drop=1.00e+00, shift=7.820e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592567655, sv_min=9.830e-01
kb(M=32)=3.141592567655, sv_min=9.830e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592613929, sv_min=9.830e-01
kb(M=32)=3.141592613929, sv_min=9.830e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592650636, sv_min=9.830e-01
kb(M=32)=3.141592650636, sv_min=9.830e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592652784, sv_min=9.830e-01
kb(M=32)=3.141592652784, sv_min=9.830e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.141199198891
kb BEM = 3.141209244415
sigma asym = 4.97191510e-02
sigma BEM = 4.90803819e-02
final sv_min(A) = 3.447e-07
final minimum drop = 1.63e+06
final mesh change in sigma = 0.000%
odd-parity residual = 9.085e-11
max first propagating-mode fraction = 3.547e-07
relative sigma error = 1.285%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.060 ===
predicted kb = 3.139600282136
predicted sigma = 1.11868090e-01
predicted kb gap to sqrt(Lambda_2) = 1.992e-03
expected-BIC mesh refinement:
M= 8: kb=3.139740307566, sigma_BEM=1.07866594e-01, sv_min=4.409e-06, drop=1.32e+05, mesh_change=--, odd_res=9.407e-11, Rprop=8.817e-07, interior=yes
M=16: kb=3.139740583902, sigma_BEM=1.07858550e-01, sv_min=5.075e-06, drop=1.14e+05, mesh_change=0.007%, odd_res=1.185e-10, Rprop=8.818e-07, interior=yes
M=24: kb=3.139740617259, sigma_BEM=1.07857579e-01, sv_min=3.760e-06, drop=1.55e+05, mesh_change=0.001%, odd_res=1.325e-10, Rprop=8.818e-07, interior=yes
M=32: kb=3.139740625339, sigma_BEM=1.07857344e-01, sv_min=3.420e-06, drop=1.70e+05, mesh_change=0.000%, odd_res=1.467e-10, Rprop=8.818e-07, interior=yes
M=40: kb=3.139740628208, sigma_BEM=1.07857261e-01, sv_min=3.297e-06, drop=1.76e+05, mesh_change=0.000%, odd_res=1.200e-10, Rprop=8.818e-07, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592375422, sv_min=9.737e-01
kb(M=32)=3.141592328758, sv_min=9.737e-01, drop=1.00e+00, shift=4.666e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592567655, sv_min=9.737e-01
kb(M=32)=3.141592567655, sv_min=9.737e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592632870, sv_min=9.737e-01
kb(M=32)=3.141592632870, sv_min=9.737e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592647935, sv_min=9.737e-01
kb(M=32)=3.141592647935, sv_min=9.737e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592652784, sv_min=9.737e-01
kb(M=32)=3.141592652784, sv_min=9.737e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.139600282136
kb BEM = 3.139740628208
sigma asym = 1.11868090e-01
sigma BEM = 1.07857261e-01
final sv_min(A) = 3.297e-06
final minimum drop = 1.76e+05
final mesh change in sigma = 0.000%
odd-parity residual = 1.200e-10
max first propagating-mode fraction = 8.818e-07
relative sigma error = 3.585%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.080 ===
predicted kb = 3.135291453357
predicted sigma = 1.98876604e-01
predicted kb gap to sqrt(Lambda_2) = 6.301e-03
expected-BIC mesh refinement:
M= 8: kb=3.136150577190, sigma_BEM=1.84834949e-01, sv_min=1.267e-07, drop=4.82e+06, mesh_change=--, odd_res=1.824e-10, Rprop=1.769e-06, interior=yes
M=16: kb=3.136152186540, sigma_BEM=1.84807640e-01, sv_min=1.940e-06, drop=3.15e+05, mesh_change=0.015%, odd_res=1.362e-10, Rprop=1.770e-06, interior=yes
M=24: kb=3.136152322804, sigma_BEM=1.84805328e-01, sv_min=1.039e-06, drop=5.88e+05, mesh_change=0.001%, odd_res=1.519e-10, Rprop=1.770e-06, interior=yes
M=32: kb=3.136152344142, sigma_BEM=1.84804966e-01, sv_min=2.596e-07, drop=2.35e+06, mesh_change=0.000%, odd_res=1.658e-10, Rprop=1.770e-06, interior=yes
M=40: kb=3.136152339588, sigma_BEM=1.84805043e-01, sv_min=1.880e-06, drop=3.25e+05, mesh_change=0.000%, odd_res=1.197e-10, Rprop=1.770e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 7 total; 6 outside expected-BIC window
checking additional minimum 1/6
kb(M=16)=3.141592247660, sv_min=9.662e-01
kb(M=32)=3.141592325862, sv_min=9.662e-01, drop=1.00e+00, shift=7.820e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/6
kb(M=16)=3.141592580829, sv_min=9.662e-01
kb(M=32)=3.141592580829, sv_min=9.662e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/6
kb(M=16)=3.141592632870, sv_min=9.662e-01
kb(M=32)=3.141592632870, sv_min=9.662e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/6
kb(M=16)=3.141592647935, sv_min=9.662e-01
kb(M=32)=3.141592647935, sv_min=9.662e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 5/6
kb(M=16)=3.141592652046, sv_min=9.662e-01
kb(M=32)=3.141592652046, sv_min=9.662e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 6/6
kb(M=16)=3.141592653169, sv_min=9.662e-01
kb(M=32)=3.141592653169, sv_min=9.662e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.135291453357
kb BEM = 3.136152339588
sigma asym = 1.98876604e-01
sigma BEM = 1.84805043e-01
final sv_min(A) = 1.880e-06
final minimum drop = 3.25e+05
final mesh change in sigma = 0.000%
odd-parity residual = 1.197e-10
max first propagating-mode fraction = 1.770e-06
relative sigma error = 7.076%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.100 ===
predicted kb = 3.126186516580
predicted sigma = 3.10744694e-01
predicted kb gap to sqrt(Lambda_2) = 1.541e-02
expected-BIC mesh refinement:
M= 8: kb=3.129559405703, sigma_BEM=2.74704072e-01, sv_min=6.895e-07, drop=9.47e+05, mesh_change=--, odd_res=1.805e-10, Rprop=3.147e-06, interior=yes
M=16: kb=3.129564709124, sigma_BEM=2.74643646e-01, sv_min=5.185e-07, drop=1.26e+06, mesh_change=0.022%, odd_res=1.774e-10, Rprop=3.158e-06, interior=yes
M=24: kb=3.129565197957, sigma_BEM=2.74638076e-01, sv_min=3.891e-07, drop=1.68e+06, mesh_change=0.002%, odd_res=1.926e-10, Rprop=3.157e-06, interior=yes
M=32: kb=3.129565315032, sigma_BEM=2.74636742e-01, sv_min=3.563e-07, drop=1.83e+06, mesh_change=0.000%, odd_res=1.659e-10, Rprop=3.157e-06, interior=yes
M=40: kb=3.129565356471, sigma_BEM=2.74636270e-01, sv_min=3.444e-07, drop=1.90e+06, mesh_change=0.000%, odd_res=1.503e-10, Rprop=3.157e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 7 total; 6 outside expected-BIC window
checking additional minimum 1/6
kb(M=16)=3.141592167792, sv_min=9.622e-01
kb(M=32)=3.141592247660, sv_min=9.623e-01, drop=1.00e+00, shift=7.987e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/6
kb(M=16)=3.141592567655, sv_min=9.622e-01
kb(M=32)=3.141592567655, sv_min=9.623e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/6
kb(M=16)=3.141592613929, sv_min=9.622e-01
kb(M=32)=3.141592613929, sv_min=9.623e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/6
kb(M=16)=3.141592642766, sv_min=9.622e-01
kb(M=32)=3.141592642766, sv_min=9.623e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 5/6
kb(M=16)=3.141592650636, sv_min=9.622e-01
kb(M=32)=3.141592650636, sv_min=9.623e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 6/6
kb(M=16)=3.141592653169, sv_min=9.622e-01
kb(M=32)=3.141592653169, sv_min=9.623e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.126186516580
kb BEM = 3.129565356471
sigma asym = 3.10744694e-01
sigma BEM = 2.74636270e-01
final sv_min(A) = 3.444e-07
final minimum drop = 1.90e+06
final mesh change in sigma = 0.000%
odd-parity residual = 1.503e-10
max first propagating-mode fraction = 3.157e-06
relative sigma error = 11.620%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.120 ===
predicted kb = 3.109561526827
predicted sigma = 4.47472359e-01
predicted kb gap to sqrt(Lambda_2) = 3.203e-02
expected-BIC mesh refinement:
M= 8: kb=3.119537324220, sigma_BEM=3.71606356e-01, sv_min=6.955e-07, drop=1.02e+06, mesh_change=--, odd_res=1.925e-10, Rprop=5.145e-06, interior=yes
M=16: kb=3.119550817416, sigma_BEM=3.71493067e-01, sv_min=4.381e-08, drop=1.62e+07, mesh_change=0.030%, odd_res=1.795e-10, Rprop=5.186e-06, interior=yes
M=24: kb=3.119552029854, sigma_BEM=3.71482885e-01, sv_min=3.979e-07, drop=1.78e+06, mesh_change=0.003%, odd_res=1.769e-10, Rprop=5.186e-06, interior=yes
M=32: kb=3.119552359779, sigma_BEM=3.71480115e-01, sv_min=5.242e-07, drop=1.35e+06, mesh_change=0.001%, odd_res=1.838e-10, Rprop=5.186e-06, interior=yes
M=40: kb=3.119552445641, sigma_BEM=3.71479394e-01, sv_min=6.692e-08, drop=1.06e+07, mesh_change=0.000%, odd_res=1.750e-10, Rprop=5.186e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592491047, sv_min=9.633e-01
kb(M=32)=3.141592491047, sv_min=9.634e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592580829, sv_min=9.633e-01
kb(M=32)=3.141592580829, sv_min=9.634e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592647935, sv_min=9.633e-01
kb(M=32)=3.141592647935, sv_min=9.634e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652784, sv_min=9.633e-01
kb(M=32)=3.141592652784, sv_min=9.634e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.109561526827
kb BEM = 3.119552445641
sigma asym = 4.47472359e-01
sigma BEM = 3.71479394e-01
final sv_min(A) = 6.692e-08
final minimum drop = 1.06e+07
final mesh change in sigma = 0.000%
odd-parity residual = 1.750e-10
max first propagating-mode fraction = 5.186e-06
relative sigma error = 16.983%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.140 ===
predicted kb = 3.081988125420
predicted sigma = 6.09059600e-01
predicted kb gap to sqrt(Lambda_2) = 5.960e-02
expected-BIC mesh refinement:
M= 8: kb=3.106247990410, sigma_BEM=4.69923210e-01, sv_min=2.980e-07, drop=2.59e+06, mesh_change=--, odd_res=2.341e-10, Rprop=7.832e-06, interior=yes
M=16: kb=3.106276103192, sigma_BEM=4.69737343e-01, sv_min=1.867e-07, drop=4.14e+06, mesh_change=0.040%, odd_res=1.939e-10, Rprop=7.963e-06, interior=yes
M=24: kb=3.106278638203, sigma_BEM=4.69720580e-01, sv_min=2.368e-07, drop=3.26e+06, mesh_change=0.004%, odd_res=2.147e-10, Rprop=7.962e-06, interior=yes
M=32: kb=3.106279276266, sigma_BEM=4.69716360e-01, sv_min=2.702e-07, drop=2.86e+06, mesh_change=0.001%, odd_res=2.000e-10, Rprop=7.962e-06, interior=yes
M=40: kb=3.106279470798, sigma_BEM=4.69715074e-01, sv_min=6.745e-08, drop=1.15e+07, mesh_change=0.000%, odd_res=1.923e-10, Rprop=7.962e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592567655, sv_min=9.698e-01
kb(M=32)=3.141592567655, sv_min=9.698e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592632870, sv_min=9.698e-01
kb(M=32)=3.141592632870, sv_min=9.698e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592650636, sv_min=9.698e-01
kb(M=32)=3.141592650636, sv_min=9.698e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652784, sv_min=9.698e-01
kb(M=32)=3.141592652784, sv_min=9.698e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.081988125420
kb BEM = 3.106279470798
sigma asym = 6.09059600e-01
sigma BEM = 4.69715074e-01
final sv_min(A) = 6.745e-08
final minimum drop = 1.15e+07
final mesh change in sigma = 0.000%
odd-parity residual = 1.923e-10
max first propagating-mode fraction = 7.962e-06
relative sigma error = 22.879%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.160 ===
predicted kb = 3.039206137054
predicted sigma = 7.95506416e-01
predicted kb gap to sqrt(Lambda_2) = 1.024e-01
expected-BIC mesh refinement:
M= 8: kb=3.090382017741, sigma_BEM=5.64927770e-01, sv_min=1.947e-08, drop=3.70e+07, mesh_change=--, odd_res=2.171e-10, Rprop=1.203e-05, interior=yes
M=16: kb=3.090432552897, sigma_BEM=5.64651253e-01, sv_min=4.538e-07, drop=1.59e+06, mesh_change=0.049%, odd_res=2.205e-10, Rprop=1.151e-05, interior=yes
M=24: kb=3.090437162749, sigma_BEM=5.64626022e-01, sv_min=4.683e-07, drop=1.54e+06, mesh_change=0.004%, odd_res=2.533e-10, Rprop=1.150e-05, interior=yes
M=32: kb=3.090438243901, sigma_BEM=5.64620104e-01, sv_min=4.275e-07, drop=1.68e+06, mesh_change=0.001%, odd_res=2.541e-10, Rprop=1.150e-05, interior=yes
M=40: kb=3.090438626207, sigma_BEM=5.64618011e-01, sv_min=4.100e-07, drop=1.76e+06, mesh_change=0.000%, odd_res=2.168e-10, Rprop=1.150e-05, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 2 total; 1 outside expected-BIC window
checking additional minimum 1/1
kb(M=16)=3.141592652784, sv_min=9.631e-01
kb(M=32)=3.141592652784, sv_min=9.632e-01, drop=1.00e+00, shift=0.000e+00, odd_res=1.364e-05, Rprop=4.285e-05, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.039206137054
kb BEM = 3.090438626207
sigma asym = 7.95506416e-01
sigma BEM = 5.64618011e-01
final sv_min(A) = 4.100e-07
final minimum drop = 1.76e+06
final mesh change in sigma = 0.000%
odd-parity residual = 2.168e-10
max first propagating-mode fraction = 1.150e-05
relative sigma error = 29.024%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.180 ===
predicted kb = 2.975891861568
predicted sigma = 1.00681281e+00
predicted kb gap to sqrt(Lambda_2) = 1.657e-01
expected-BIC mesh refinement:
M= 8: kb=3.072958191962, sigma_BEM=6.53094443e-01, sv_min=8.785e-08, drop=6.65e+06, mesh_change=--, odd_res=2.784e-10, Rprop=1.835e-05, interior=yes
M=16: kb=3.073039525908, sigma_BEM=6.52711631e-01, sv_min=1.806e-07, drop=3.23e+06, mesh_change=0.059%, odd_res=2.739e-10, Rprop=1.570e-05, interior=yes
M=24: kb=3.073046734408, sigma_BEM=6.52677692e-01, sv_min=9.100e-08, drop=6.41e+06, mesh_change=0.005%, odd_res=2.554e-10, Rprop=1.570e-05, interior=yes
M=32: kb=3.073048486232, sigma_BEM=6.52669443e-01, sv_min=1.025e-07, drop=5.69e+06, mesh_change=0.001%, odd_res=2.492e-10, Rprop=1.570e-05, interior=yes
M=40: kb=3.073049078815, sigma_BEM=6.52666653e-01, sv_min=8.486e-08, drop=6.87e+06, mesh_change=0.000%, odd_res=2.483e-10, Rprop=1.570e-05, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 2 total; 1 outside expected-BIC window
checking additional minimum 1/1
kb(M=16)=3.141592653169, sv_min=9.500e-01
kb(M=32)=3.141592653169, sv_min=9.501e-01, drop=1.00e+00, shift=0.000e+00, odd_res=1.255e-05, Rprop=4.642e-05, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 2.975891861568
kb BEM = 3.073049078815
sigma asym = 1.00681281e+00
sigma BEM = 6.52666653e-01
final sv_min(A) = 8.486e-08
final minimum drop = 6.87e+06
final mesh change in sigma = 0.000%
odd-parity residual = 2.483e-10
max first propagating-mode fraction = 1.570e-05
relative sigma error = 35.175%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.200 ===
predicted kb = 2.885239706985
predicted sigma = 1.24297877e+00
predicted kb gap to sqrt(Lambda_2) = 2.564e-01
expected-BIC mesh refinement:
M= 8: kb=3.055100750025, sigma_BEM=7.32095491e-01, sv_min=3.050e-07, drop=1.42e+06, mesh_change=--, odd_res=2.950e-10, Rprop=2.818e-05, interior=yes
M=16: kb=3.055220844812, sigma_BEM=7.31594143e-01, sv_min=1.663e-07, drop=2.60e+06, mesh_change=0.069%, odd_res=2.743e-10, Rprop=2.028e-05, interior=yes
M=24: kb=3.055231441409, sigma_BEM=7.31549889e-01, sv_min=1.762e-07, drop=2.46e+06, mesh_change=0.006%, odd_res=2.817e-10, Rprop=2.028e-05, interior=yes
M=32: kb=3.055234009102, sigma_BEM=7.31539166e-01, sv_min=8.503e-08, drop=5.09e+06, mesh_change=0.001%, odd_res=2.689e-10, Rprop=2.028e-05, interior=yes
M=40: kb=3.055234910113, sigma_BEM=7.31535403e-01, sv_min=1.197e-07, drop=3.62e+06, mesh_change=0.001%, odd_res=2.837e-10, Rprop=2.028e-05, interior=yes
whole-band uniqueness screen: 78 points (M=16)
16/78
32/78
48/78
64/78
78/78
sampled local minima: 1 total; 0 outside expected-BIC window

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 2.885239706985
kb BEM = 3.055234910113
sigma asym = 1.24297877e+00
sigma BEM = 7.31535403e-01
final sv_min(A) = 1.197e-07
final minimum drop = 3.62e+06
final mesh change in sigma = 0.001%
odd-parity residual = 2.837e-10
max first propagating-mode fraction = 2.028e-05
relative sigma error = 41.147%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 9/10
Leading asymptotic sigma within 5.0%: 3/10
epsilon=0.020: BIC=FAIL, asymptotic=PASS, sv_min=1.076e-04, mesh_change=0.000%, odd_res=2.192e-10, Rprop=8.282e-08, sigma_error=0.151%
epsilon=0.040: BIC=PASS, asymptotic=PASS, sv_min=3.447e-07, mesh_change=0.000%, odd_res=9.085e-11, Rprop=3.547e-07, sigma_error=1.285%
epsilon=0.060: BIC=PASS, asymptotic=PASS, sv_min=3.297e-06, mesh_change=0.000%, odd_res=1.200e-10, Rprop=8.818e-07, sigma_error=3.585%
epsilon=0.080: BIC=PASS, asymptotic=FAIL, sv_min=1.880e-06, mesh_change=0.000%, odd_res=1.197e-10, Rprop=1.770e-06, sigma_error=7.076%
epsilon=0.100: BIC=PASS, asymptotic=FAIL, sv_min=3.444e-07, mesh_change=0.000%, odd_res=1.503e-10, Rprop=3.157e-06, sigma_error=11.620%
epsilon=0.120: BIC=PASS, asymptotic=FAIL, sv_min=6.692e-08, mesh_change=0.000%, odd_res=1.750e-10, Rprop=5.186e-06, sigma_error=16.983%
epsilon=0.140: BIC=PASS, asymptotic=FAIL, sv_min=6.745e-08, mesh_change=0.000%, odd_res=1.923e-10, Rprop=7.962e-06, sigma_error=22.879%
epsilon=0.160: BIC=PASS, asymptotic=FAIL, sv_min=4.100e-07, mesh_change=0.000%, odd_res=2.168e-10, Rprop=1.150e-05, sigma_error=29.024%
epsilon=0.180: BIC=PASS, asymptotic=FAIL, sv_min=8.486e-08, mesh_change=0.000%, odd_res=2.483e-10, Rprop=1.570e-05, sigma_error=35.175%
epsilon=0.200: BIC=PASS, asymptotic=FAIL, sv_min=1.197e-07, mesh_change=0.001%, odd_res=2.837e-10, Rprop=2.028e-05, sigma_error=41.147%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/svd_theorems_validations/theorem_2_3_x_symmetry_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations %
