(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations % python first_theorem_2.py
=== Focused validation of Theorem 2.1 ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1) b = 1.570796326795
epsilons = (0.01, 0.03111111, 0.05222222, 0.07333333, 0.09444444, 0.11555556, 0.13666667, 0.15777778, 0.17888889, 0.2)
refinement M = (8, 16, 24, 32, 40)
Running 10 validation cases with 4 worker processes.

=== epsilon=0.010 ===
predicted kb = 1.570796282404
predicted sigma = 3.73441725e-04
predicted kb cutoff gap = 4.439e-08
expected-mode mesh refinement:
M= 8: kb=1.570796290981, sigma_BEM=3.35427662e-04, sigma_min=2.290e-03, drop=5.44e+00, mesh_change=--, interior=yes
M=16: kb=1.570796290981, sigma_BEM=3.35427662e-04, sigma_min=2.290e-03, drop=5.44e+00, mesh_change=0.000%, interior=yes
M=24: kb=1.570796290981, sigma_BEM=3.35427662e-04, sigma_min=2.290e-03, drop=5.44e+00, mesh_change=0.000%, interior=yes
M=32: kb=1.570796290981, sigma_BEM=3.35427662e-04, sigma_min=2.290e-03, drop=5.44e+00, mesh_change=0.000%, interior=yes
M=40: kb=1.570796290981, sigma_BEM=3.35427662e-04, sigma_min=2.290e-03, drop=5.44e+00, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: FAIL
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: FAIL
kb asymptotic = 1.570796282404
kb BEM = 1.570796290981
sigma asym = 3.73441725e-04
sigma BEM = 3.35427662e-04
final sigma_min(A) = 2.290e-03
final minimum drop = 5.44e+00
final mesh change in sigma = 0.000%
relative sigma error = 10.179%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: FAIL
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.031 ===
predicted kb = 1.570792168087
predicted sigma = 3.61454681e-03
predicted kb cutoff gap = 4.159e-06
expected-mode mesh refinement:
M= 8: kb=1.570792199490, sigma_BEM=3.60087412e-03, sigma_min=3.897e-05, drop=9.93e+02, mesh_change=--, interior=yes
M=16: kb=1.570792199515, sigma_BEM=3.60086326e-03, sigma_min=3.872e-05, drop=9.99e+02, mesh_change=0.000%, interior=yes
M=24: kb=1.570792199517, sigma_BEM=3.60086224e-03, sigma_min=3.870e-05, drop=1.00e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570792199518, sigma_BEM=3.60086200e-03, sigma_min=3.869e-05, drop=1.00e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570792199518, sigma_BEM=3.60086191e-03, sigma_min=3.869e-05, drop=1.00e+03, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.570792168087
kb BEM = 1.570792199518
sigma asym = 3.61454681e-03
sigma BEM = 3.60086191e-03
final sigma_min(A) = 3.869e-05
final minimum drop = 1.00e+03
final mesh change in sigma = 0.000%
relative sigma error = 0.379%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.052 ===
predicted kb = 1.570763311005
predicted sigma = 1.01843543e-02
predicted kb cutoff gap = 3.302e-05
expected-mode mesh refinement:
M= 8: kb=1.570764096702, sigma_BEM=1.00624442e-02, sigma_min=2.829e-06, drop=2.28e+04, mesh_change=--, interior=yes
M=16: kb=1.570764097021, sigma_BEM=1.00623944e-02, sigma_min=4.383e-06, drop=1.47e+04, mesh_change=0.000%, interior=yes
M=24: kb=1.570764097050, sigma_BEM=1.00623899e-02, sigma_min=4.528e-06, drop=1.43e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570764097057, sigma_BEM=1.00623888e-02, sigma_min=4.563e-06, drop=1.41e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570764097060, sigma_BEM=1.00623884e-02, sigma_min=4.575e-06, drop=1.41e+04, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.570763311005
kb BEM = 1.570764097060
sigma asym = 1.01843543e-02
sigma BEM = 1.00623884e-02
final sigma_min(A) = 4.575e-06
final minimum drop = 1.41e+04
final mesh change in sigma = 0.000%
relative sigma error = 1.198%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.073 ===
predicted kb = 1.570667940347
predicted sigma = 2.00828643e-02
predicted kb cutoff gap = 1.284e-04
expected-mode mesh refinement:
M= 8: kb=1.570674819672, sigma_BEM=1.95374287e-02, sigma_min=4.682e-06, drop=1.91e+04, mesh_change=--, interior=yes
M=16: kb=1.570674810399, sigma_BEM=1.95381742e-02, sigma_min=7.376e-06, drop=1.21e+04, mesh_change=0.004%, interior=yes
M=24: kb=1.570674812881, sigma_BEM=1.95379746e-02, sigma_min=6.242e-06, drop=1.43e+04, mesh_change=0.001%, interior=yes
M=32: kb=1.570674813467, sigma_BEM=1.95379276e-02, sigma_min=5.976e-06, drop=1.50e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570674813677, sigma_BEM=1.95379106e-02, sigma_min=5.880e-06, drop=1.52e+04, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.570667940347
kb BEM = 1.570674813677
sigma asym = 2.00828643e-02
sigma BEM = 1.95379106e-02
final sigma_min(A) = 5.880e-06
final minimum drop = 1.52e+04
final mesh change in sigma = 0.000%
relative sigma error = 2.714%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.094 ===
predicted kb = 1.570443102779
predicted sigma = 3.33100766e-02
predicted kb cutoff gap = 3.532e-04
expected-mode mesh refinement:
M= 8: kb=1.570476709914, sigma_BEM=3.16860204e-02, sigma_min=3.749e-07, drop=3.00e+05, mesh_change=--, interior=yes
M=16: kb=1.570476746204, sigma_BEM=3.16842217e-02, sigma_min=6.924e-07, drop=1.63e+05, mesh_change=0.006%, interior=yes
M=24: kb=1.570476748255, sigma_BEM=3.16841200e-02, sigma_min=1.139e-06, drop=9.88e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570476748729, sigma_BEM=3.16840966e-02, sigma_min=1.252e-06, drop=8.99e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570476748894, sigma_BEM=3.16840884e-02, sigma_min=1.292e-06, drop=8.71e+04, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.570443102779
kb BEM = 1.570476748894
sigma asym = 3.33100766e-02
sigma BEM = 3.16840884e-02
final sigma_min(A) = 1.292e-06
final minimum drop = 8.71e+04
final mesh change in sigma = 0.000%
relative sigma error = 4.881%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.116 ===
predicted kb = 1.570004612194
predicted sigma = 4.98660001e-02
predicted kb cutoff gap = 7.917e-04
expected-mode mesh refinement:
M= 8: kb=1.570121256196, sigma_BEM=4.60471618e-02, sigma_min=1.347e-06, drop=9.89e+04, mesh_change=--, interior=yes
M=16: kb=1.570121370744, sigma_BEM=4.60432558e-02, sigma_min=9.898e-07, drop=1.35e+05, mesh_change=0.008%, interior=yes
M=24: kb=1.570121367002, sigma_BEM=4.60433834e-02, sigma_min=1.686e-06, drop=7.90e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570121390178, sigma_BEM=4.60425931e-02, sigma_min=2.163e-06, drop=6.16e+04, mesh_change=0.002%, interior=yes
M=40: kb=1.570121390176, sigma_BEM=4.60425932e-02, sigma_min=1.995e-06, drop=6.68e+04, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.570004612194
kb BEM = 1.570121390176
sigma asym = 4.98660001e-02
sigma BEM = 4.60425932e-02
final sigma_min(A) = 1.995e-06
final minimum drop = 6.68e+04
final mesh change in sigma = 0.000%
relative sigma error = 7.667%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.137 ===
predicted kb = 1.569246937686
predicted sigma = 6.97506189e-02
predicted kb cutoff gap = 1.549e-03
expected-mode mesh refinement:
M= 8: kb=1.569569154130, sigma_BEM=6.20787458e-02, sigma_min=8.353e-07, drop=1.80e+05, mesh_change=--, interior=yes
M=16: kb=1.569569441664, sigma_BEM=6.20714755e-02, sigma_min=3.438e-08, drop=4.38e+06, mesh_change=0.012%, interior=yes
M=24: kb=1.569569467463, sigma_BEM=6.20708231e-02, sigma_min=3.005e-08, drop=5.01e+06, mesh_change=0.001%, interior=yes
M=32: kb=1.569569473634, sigma_BEM=6.20706671e-02, sigma_min=4.535e-08, drop=3.32e+06, mesh_change=0.000%, interior=yes
M=40: kb=1.569569475818, sigma_BEM=6.20706118e-02, sigma_min=5.043e-08, drop=2.99e+06, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.569246937686
kb BEM = 1.569569475818
sigma asym = 6.97506189e-02
sigma BEM = 6.20706118e-02
final sigma_min(A) = 5.043e-08
final minimum drop = 2.99e+06
final mesh change in sigma = 0.000%
relative sigma error = 11.011%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.158 ===
predicted kb = 1.568042986053
predicted sigma = 9.29639401e-02
predicted kb cutoff gap = 2.753e-03
expected-mode mesh refinement:
M= 8: kb=1.568798942026, sigma_BEM=7.91895181e-02, sigma_min=9.222e-07, drop=1.78e+05, mesh_change=--, interior=yes
M=16: kb=1.568799583860, sigma_BEM=7.91768019e-02, sigma_min=9.268e-07, drop=1.77e+05, mesh_change=0.016%, interior=yes
M=24: kb=1.568799637637, sigma_BEM=7.91757363e-02, sigma_min=7.771e-07, drop=2.11e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.568799653362, sigma_BEM=7.91754248e-02, sigma_min=9.797e-07, drop=1.67e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.568799658590, sigma_BEM=7.91753212e-02, sigma_min=1.024e-06, drop=1.60e+05, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.568042986053
kb BEM = 1.568799658590
sigma asym = 9.29639401e-02
sigma BEM = 7.91753212e-02
final sigma_min(A) = 1.024e-06
final minimum drop = 1.60e+05
final mesh change in sigma = 0.000%
relative sigma error = 14.832%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.179 ===
predicted kb = 1.566243730998
predicted sigma = 1.19505964e-01
predicted kb cutoff gap = 4.553e-03
expected-mode mesh refinement:
M= 8: kb=1.567812045701, sigma_BEM=9.67806263e-02, sigma_min=7.165e-08, drop=2.40e+06, mesh_change=--, interior=yes
M=16: kb=1.567813309825, sigma_BEM=9.67601458e-02, sigma_min=3.391e-08, drop=5.07e+06, mesh_change=0.021%, interior=yes
M=24: kb=1.567813404855, sigma_BEM=9.67586060e-02, sigma_min=4.861e-07, drop=3.54e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.567813427595, sigma_BEM=9.67582375e-02, sigma_min=6.088e-07, drop=2.83e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.567813435642, sigma_BEM=9.67581071e-02, sigma_min=6.517e-07, drop=2.64e+05, mesh_change=0.000%, interior=yes
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
kb asymptotic = 1.566243730998
kb BEM = 1.567813435642
sigma asym = 1.19505964e-01
sigma BEM = 9.67581071e-02
final sigma_min(A) = 6.517e-07
final minimum drop = 2.64e+05
final mesh change in sigma = 0.000%
relative sigma error = 19.035%
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
Paper sweep points: 120; resolved modes: 77

=== SUBCRITICAL FULL-BAND CHECK ===
leading-order threshold a0* = 0.391826552031
a=0.350000 (< a0*): PASS; full-band points=20, resolved modes=0, statuses=['no_resolved_mode']

=== FINAL SUMMARY ===
Interval screened: 0 < k^2 < Lambda_1
Existence + no additional resolved modes: 9/10
Leading asymptotic sigma within 5.0%: 4/10
epsilon=0.010: unique=FAIL, asymptotic=FAIL, sigma_min=2.290e-03, mesh_change=0.000%, sigma_error=10.179%
epsilon=0.031: unique=PASS, asymptotic=PASS, sigma_min=3.869e-05, mesh_change=0.000%, sigma_error=0.379%
epsilon=0.052: unique=PASS, asymptotic=PASS, sigma_min=4.575e-06, mesh_change=0.000%, sigma_error=1.198%
epsilon=0.073: unique=PASS, asymptotic=PASS, sigma_min=5.880e-06, mesh_change=0.000%, sigma_error=2.714%
epsilon=0.094: unique=PASS, asymptotic=PASS, sigma_min=1.292e-06, mesh_change=0.000%, sigma_error=4.881%
epsilon=0.116: unique=PASS, asymptotic=FAIL, sigma_min=1.995e-06, mesh_change=0.000%, sigma_error=7.667%
epsilon=0.137: unique=PASS, asymptotic=FAIL, sigma_min=5.043e-08, mesh_change=0.000%, sigma_error=11.011%
epsilon=0.158: unique=PASS, asymptotic=FAIL, sigma_min=1.024e-06, mesh_change=0.000%, sigma_error=14.832%
epsilon=0.179: unique=PASS, asymptotic=FAIL, sigma_min=6.517e-07, mesh_change=0.000%, sigma_error=19.035%
epsilon=0.200: unique=PASS, asymptotic=FAIL, sigma_min=5.417e-07, mesh_change=0.000%, sigma_error=23.519%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/svd_theorems_validations/theorem_2_1_refined_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/svd_theorems_validations %
