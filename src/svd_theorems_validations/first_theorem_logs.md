(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations % python first_theorem_2.py
=== Focused validation of Theorem 2.1 ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1) b = 1.570796326795
epsilons = (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
refinement M = (8, 16, 24, 32, 40)
Running 10 validation cases with 4 worker processes.

=== epsilon=0.020 ===
predicted kb = 1.570795616537
predicted sigma = 1.49376690e-03
predicted kb cutoff gap = 7.103e-07
expected-mode mesh refinement:
M= 8: kb=1.570795616433, sigma_BEM=1.49387665e-03, sigma_min=4.667e-05, drop=5.33e+02, mesh_change=--, interior=yes
M=16: kb=1.570795616435, sigma_BEM=1.49387467e-03, sigma_min=4.673e-05, drop=5.33e+02, mesh_change=0.000%, interior=yes
M=24: kb=1.570795616435, sigma_BEM=1.49387449e-03, sigma_min=4.673e-05, drop=5.33e+02, mesh_change=0.000%, interior=yes
M=32: kb=1.570795616435, sigma_BEM=1.49387446e-03, sigma_min=4.674e-05, drop=5.33e+02, mesh_change=0.000%, interior=yes
M=40: kb=1.570795616435, sigma_BEM=1.49387444e-03, sigma_min=4.674e-05, drop=5.33e+02, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570795616537
kb BEM = 1.570795616435
sigma asym = 1.49376690e-03
sigma BEM = 1.49387444e-03
final sigma_min(A) = 4.674e-05
final minimum drop = 5.33e+02
final mesh change in sigma = 0.000%
relative sigma error = 0.007%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.040 ===
predicted kb = 1.570784962635
predicted sigma = 5.97506760e-03
predicted kb cutoff gap = 1.136e-05
expected-mode mesh refinement:
M= 8: kb=1.570785105295, sigma_BEM=5.93744516e-03, sigma_min=1.132e-05, drop=4.39e+03, mesh_change=--, interior=yes
M=16: kb=1.570785105658, sigma_BEM=5.93734919e-03, sigma_min=1.181e-05, drop=4.20e+03, mesh_change=0.002%, interior=yes
M=24: kb=1.570785105691, sigma_BEM=5.93734042e-03, sigma_min=1.185e-05, drop=4.19e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570785105699, sigma_BEM=5.93733835e-03, sigma_min=1.186e-05, drop=4.19e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570785105702, sigma_BEM=5.93733757e-03, sigma_min=1.187e-05, drop=4.18e+03, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570784962635
kb BEM = 1.570785105702
sigma asym = 5.97506760e-03
sigma BEM = 5.93733757e-03
final sigma_min(A) = 1.187e-05
final minimum drop = 4.18e+03
final mesh change in sigma = 0.000%
relative sigma error = 0.631%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.060 ===
predicted kb = 1.570738794889
predicted sigma = 1.34439021e-02
predicted kb cutoff gap = 5.753e-05
expected-mode mesh refinement:
M= 8: kb=1.570740720078, sigma_BEM=1.32170557e-02, sigma_min=8.255e-06, drop=8.95e+03, mesh_change=--, interior=yes
M=16: kb=1.570740723230, sigma_BEM=1.32166811e-02, sigma_min=8.844e-06, drop=8.35e+03, mesh_change=0.003%, interior=yes
M=24: kb=1.570740723522, sigma_BEM=1.32166464e-02, sigma_min=8.900e-06, drop=8.30e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570740723591, sigma_BEM=1.32166382e-02, sigma_min=8.912e-06, drop=8.29e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570740723616, sigma_BEM=1.32166352e-02, sigma_min=8.918e-06, drop=8.28e+03, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570738794889
kb BEM = 1.570740723616
sigma asym = 1.34439021e-02
sigma BEM = 1.32166352e-02
final sigma_min(A) = 8.918e-06
final minimum drop = 8.28e+03
final mesh change in sigma = 0.000%
relative sigma error = 1.690%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.080 ===
predicted kb = 1.570614490366
predicted sigma = 2.39002704e-02
predicted kb cutoff gap = 1.818e-04
expected-mode mesh refinement:
M= 8: kb=1.570626393040, sigma_BEM=2.31048427e-02, sigma_min=2.688e-06, drop=3.61e+04, mesh_change=--, interior=yes
M=16: kb=1.570626408042, sigma_BEM=2.31038229e-02, sigma_min=3.036e-06, drop=3.19e+04, mesh_change=0.004%, interior=yes
M=24: kb=1.570626409204, sigma_BEM=2.31037439e-02, sigma_min=2.953e-06, drop=3.28e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570626409482, sigma_BEM=2.31037250e-02, sigma_min=2.934e-06, drop=3.30e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570626409581, sigma_BEM=2.31037183e-02, sigma_min=2.927e-06, drop=3.31e+04, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570614490366
kb BEM = 1.570626409581
sigma asym = 2.39002704e-02
sigma BEM = 2.31037183e-02
final sigma_min(A) = 2.927e-06
final minimum drop = 3.31e+04
final mesh change in sigma = 0.000%
relative sigma error = 3.333%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.100 ===
predicted kb = 1.570352353153
predicted sigma = 3.73441725e-02
predicted kb cutoff gap = 4.440e-04
expected-mode mesh refinement:
M= 8: kb=1.570400264377, sigma_BEM=3.52719423e-02, sigma_min=1.799e-06, drop=6.57e+04, mesh_change=--, interior=yes
M=16: kb=1.570400322563, sigma_BEM=3.52693517e-02, sigma_min=3.102e-08, drop=3.81e+06, mesh_change=0.007%, interior=yes
M=24: kb=1.570400325582, sigma_BEM=3.52692173e-02, sigma_min=4.444e-07, drop=2.66e+05, mesh_change=0.000%, interior=yes
M=32: kb=1.570400326308, sigma_BEM=3.52691849e-02, sigma_min=5.572e-07, drop=2.12e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.570400326566, sigma_BEM=3.52691734e-02, sigma_min=5.967e-07, drop=1.98e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570352353153
kb BEM = 1.570400326566
sigma asym = 3.73441725e-02
sigma BEM = 3.52691734e-02
final sigma_min(A) = 5.967e-07
final minimum drop = 1.98e+05
final mesh change in sigma = 0.000%
relative sigma error = 5.556%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.120 ===
predicted kb = 1.569875563291
predicted sigma = 5.37756084e-02
predicted kb cutoff gap = 9.208e-04
expected-mode mesh refinement:
M= 8: kb=1.570022381932, sigma_BEM=4.93033519e-02, sigma_min=1.684e-06, drop=8.15e+04, mesh_change=--, interior=yes
M=16: kb=1.570022529420, sigma_BEM=4.92986550e-02, sigma_min=1.046e-06, drop=1.31e+05, mesh_change=0.010%, interior=yes
M=24: kb=1.570022542819, sigma_BEM=4.92982283e-02, sigma_min=9.876e-07, drop=1.39e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.570022546024, sigma_BEM=4.92981262e-02, sigma_min=9.737e-07, drop=1.41e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.570022547158, sigma_BEM=4.92980901e-02, sigma_min=9.687e-07, drop=1.42e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.569875563291
kb BEM = 1.570022547158
sigma asym = 5.37756084e-02
sigma BEM = 4.92980901e-02
final sigma_min(A) = 9.687e-07
final minimum drop = 1.42e+05
final mesh change in sigma = 0.000%
relative sigma error = 8.326%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.140 ===
predicted kb = 1.569090071990
predicted sigma = 7.31945781e-02
predicted kb cutoff gap = 1.706e-03
expected-mode mesh refinement:
M= 8: kb=1.569462278930, sigma_BEM=6.47244567e-02, sigma_min=5.187e-07, drop=2.95e+05, mesh_change=--, interior=yes
M=16: kb=1.569462609796, sigma_BEM=6.47164332e-02, sigma_min=1.053e-06, drop=1.45e+05, mesh_change=0.012%, interior=yes
M=24: kb=1.569462638662, sigma_BEM=6.47157331e-02, sigma_min=1.172e-06, drop=1.30e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.569462645558, sigma_BEM=6.47155659e-02, sigma_min=1.201e-06, drop=1.27e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.569462647994, sigma_BEM=6.47155068e-02, sigma_min=1.211e-06, drop=1.26e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.569090071990
kb BEM = 1.569462647994
sigma asym = 7.31945781e-02
sigma BEM = 6.47155068e-02
final sigma_min(A) = 1.211e-06
final minimum drop = 1.26e+05
final mesh change in sigma = 0.000%
relative sigma error = 11.584%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.160 ===
predicted kb = 1.567884413304
predicted sigma = 9.56010815e-02
predicted kb cutoff gap = 2.912e-03
expected-mode mesh refinement:
M= 8: kb=1.568704989728, sigma_BEM=8.10293495e-02, sigma_min=1.189e-06, drop=1.39e+05, mesh_change=--, interior=yes
M=16: kb=1.568705687777, sigma_BEM=8.10158343e-02, sigma_min=7.115e-07, drop=2.32e+05, mesh_change=0.017%, interior=yes
M=24: kb=1.568705750589, sigma_BEM=8.10146181e-02, sigma_min=4.581e-07, drop=3.60e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.568705780738, sigma_BEM=8.10140343e-02, sigma_min=8.170e-07, drop=2.02e+05, mesh_change=0.001%, interior=yes
M=40: kb=1.568705764855, sigma_BEM=8.10143418e-02, sigma_min=8.625e-07, drop=1.91e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.567884413304
kb BEM = 1.568705764855
sigma asym = 9.56010815e-02
sigma BEM = 8.10143418e-02
final sigma_min(A) = 8.625e-07
final minimum drop = 1.91e+05
final mesh change in sigma = 0.000%
relative sigma error = 15.258%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.180 ===
predicted kb = 1.566129394876
predicted sigma = 1.20995119e-01
predicted kb cutoff gap = 4.667e-03
expected-mode mesh refinement:
M= 8: kb=1.567754480597, sigma_BEM=9.77086938e-02, sigma_min=2.483e-08, drop=6.94e+06, mesh_change=--, interior=yes
M=16: kb=1.567755777310, sigma_BEM=9.76878856e-02, sigma_min=5.422e-07, drop=3.18e+05, mesh_change=0.021%, interior=yes
M=24: kb=1.567755900645, sigma_BEM=9.76859062e-02, sigma_min=4.762e-07, drop=3.62e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.567755925665, sigma_BEM=9.76855046e-02, sigma_min=4.486e-07, drop=3.84e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.567755934513, sigma_BEM=9.76853626e-02, sigma_min=4.385e-07, drop=3.93e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.566129394876
kb BEM = 1.567755934513
sigma asym = 1.20995119e-01
sigma BEM = 9.76853626e-02
final sigma_min(A) = 4.385e-07
final minimum drop = 3.93e+05
final mesh change in sigma = 0.000%
relative sigma error = 19.265%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.200 ===
predicted kb = 1.563677621758
predicted sigma = 1.49376690e-01
predicted kb cutoff gap = 7.119e-03
expected-mode mesh refinement:
M= 8: kb=1.566633761475, sigma_BEM=1.14279297e-01, sigma_min=3.992e-08, drop=4.38e+06, mesh_change=--, interior=yes
M=16: kb=1.566636065119, sigma_BEM=1.14247712e-01, sigma_min=8.745e-08, drop=2.00e+06, mesh_change=0.028%, interior=yes
M=24: kb=1.566636251307, sigma_BEM=1.14245159e-01, sigma_min=4.863e-07, drop=3.59e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.566636291663, sigma_BEM=1.14244606e-01, sigma_min=4.238e-07, drop=4.12e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.566636286380, sigma_BEM=1.14244678e-01, sigma_min=5.417e-07, drop=3.22e+05, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.563677621758
kb BEM = 1.566636286380
sigma asym = 1.49376690e-01
sigma BEM = 1.14244678e-01
final sigma_min(A) = 5.417e-07
final minimum drop = 3.22e+05
final mesh change in sigma = 0.000%
relative sigma error = 23.519%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== PAPER FIGURE 2: sweep in a ===
a values = (0.35, 0.5, 0.55, 0.6, 0.65, 0.7)
fixed M = 32
Paper sweep points: 95; resolved modes: 73

=== SUBCRITICAL FULL-BAND CHECK ===
leading-order threshold a0* = 0.391826552031
a=0.350000 (< a0*): PASS; full-band points=20, resolved modes=0, statuses=['no_resolved_mode']

=== FINAL SUMMARY ===
Interval screened: 0 < k^2 < Lambda_1
Existence + no additional resolved modes: 10/10
Leading asymptotic sigma within 5.0%: 4/10
epsilon=0.020: unique=PASS, asymptotic=PASS, sigma_min=4.674e-05, mesh_change=0.000%, sigma_error=0.007%
epsilon=0.040: unique=PASS, asymptotic=PASS, sigma_min=1.187e-05, mesh_change=0.000%, sigma_error=0.631%
epsilon=0.060: unique=PASS, asymptotic=PASS, sigma_min=8.918e-06, mesh_change=0.000%, sigma_error=1.690%
epsilon=0.080: unique=PASS, asymptotic=PASS, sigma_min=2.927e-06, mesh_change=0.000%, sigma_error=3.333%
epsilon=0.100: unique=PASS, asymptotic=FAIL, sigma_min=5.967e-07, mesh_change=0.000%, sigma_error=5.556%
epsilon=0.120: unique=PASS, asymptotic=FAIL, sigma_min=9.687e-07, mesh_change=0.000%, sigma_error=8.326%
epsilon=0.140: unique=PASS, asymptotic=FAIL, sigma_min=1.211e-06, mesh_change=0.000%, sigma_error=11.584%
epsilon=0.160: unique=PASS, asymptotic=FAIL, sigma_min=8.625e-07, mesh_change=0.000%, sigma_error=15.258%
epsilon=0.180: unique=PASS, asymptotic=FAIL, sigma_min=4.385e-07, mesh_change=0.000%, sigma_error=19.265%
epsilon=0.200: unique=PASS, asymptotic=FAIL, sigma_min=5.417e-07, mesh_change=0.000%, sigma_error=23.519%
