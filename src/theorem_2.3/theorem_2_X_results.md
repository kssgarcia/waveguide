---
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python second_theorem_1.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iii): x-symmetric BIC ===
b = 1.0
a = 0.0
beta = 0.1
mu = 1.002199318791e+00
nu = -1.687599000907e-17
X-even residual = 0.000e+00
Y-odd residual  = 0.000e+00
Theorem symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590
epsilons = (0.05, 0.07, 0.09, 0.11)
refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
  predicted kb = 3.140631984100
  predicted sigma = 7.76861734e-02
  predicted kb gap to sqrt(Lambda_2) = 9.607e-04
  expected-BIC mesh refinement:
    M=16: kb=3.140675347570, sigma_BEM=7.59128596e-02, sv_min=4.455e-07, drop=1.28e+06, mesh_change=--, odd_res=1.035e-10, Rprop=5.806e-07, interior=yes
    M=24: kb=3.140675357638, sigma_BEM=7.59124431e-02, sv_min=5.307e-07, drop=1.07e+06, mesh_change=0.001%, odd_res=1.150e-10, Rprop=5.806e-07, interior=yes
    M=32: kb=3.140675360048, sigma_BEM=7.59123434e-02, sv_min=5.509e-07, drop=1.04e+06, mesh_change=0.000%, odd_res=1.062e-10, Rprop=5.806e-07, interior=yes
    M=40: kb=3.140675360901, sigma_BEM=7.59123081e-02, sv_min=5.580e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.104e-10, Rprop=5.806e-07, interior=yes
    M=48: kb=3.140675361277, sigma_BEM=7.59122925e-02, sv_min=5.611e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.255e-10, Rprop=5.806e-07, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141592488189, sv_min=9.782e-01
      kb(M=32)=3.141592441525, sv_min=9.782e-01, drop=1.00e+00, shift=4.666e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592580829, sv_min=9.782e-01
      kb(M=32)=3.141592580829, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592642766, sv_min=9.782e-01
      kb(M=32)=3.141592642766, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592652046, sv_min=9.782e-01
      kb(M=32)=3.141592652046, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

  --- validation result ---
  symmetry conditions (a=0, X even, Y odd, nu~0): PASS
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.140631984100
  kb BEM        = 3.140675361277
  sigma asym    = 7.76861734e-02
  sigma BEM     = 7.59122925e-02
  final sv_min(A) = 5.611e-07
  final minimum drop = 1.02e+06
  final mesh change in sigma = 0.000%
  odd-parity residual = 1.255e-10
  max first propagating-mode fraction = 5.806e-07
  relative sigma error = 2.283%
  asymptotic accuracy <= 5.0%: PASS
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.070 ===
  predicted kb = 3.137900540385
  predicted sigma = 1.52264900e-01
  predicted kb gap to sqrt(Lambda_2) = 3.692e-03
  expected-BIC mesh refinement:
    M=16: kb=3.138273756691, sigma_BEM=1.44368380e-01, sv_min=1.959e-06, drop=3.03e+05, mesh_change=--, odd_res=1.651e-10, Rprop=1.273e-06, interior=yes
    M=24: kb=3.138273862981, sigma_BEM=1.44366069e-01, sv_min=3.835e-06, drop=1.55e+05, mesh_change=0.002%, odd_res=1.160e-10, Rprop=1.273e-06, interior=yes
    M=32: kb=3.138273870562, sigma_BEM=1.44365905e-01, sv_min=2.448e-06, drop=2.43e+05, mesh_change=0.000%, odd_res=1.716e-10, Rprop=1.273e-06, interior=yes
    M=40: kb=3.138273873527, sigma_BEM=1.44365840e-01, sv_min=2.001e-06, drop=2.97e+05, mesh_change=0.000%, odd_res=1.247e-10, Rprop=1.273e-06, interior=yes
    M=48: kb=3.138273874883, sigma_BEM=1.44365811e-01, sv_min=1.811e-06, drop=3.28e+05, mesh_change=0.000%, odd_res=1.306e-10, Rprop=1.273e-06, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141592167792, sv_min=9.696e-01
      kb(M=32)=3.141592247660, sv_min=9.696e-01, drop=1.00e+00, shift=7.987e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592491047, sv_min=9.696e-01
      kb(M=32)=3.141592491047, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592580829, sv_min=9.696e-01
      kb(M=32)=3.141592580829, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592652046, sv_min=9.696e-01
      kb(M=32)=3.141592652046, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

  --- validation result ---
  symmetry conditions (a=0, X even, Y odd, nu~0): PASS
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.137900540385
  kb BEM        = 3.138273874883
  sigma asym    = 1.52264900e-01
  sigma BEM     = 1.44365811e-01
  final sv_min(A) = 1.811e-06
  final minimum drop = 3.28e+05
  final mesh change in sigma = 0.000%
  odd-parity residual = 1.306e-10
  max first propagating-mode fraction = 1.273e-06
  relative sigma error = 5.188%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.090 ===
  predicted kb = 3.131493237939
  predicted sigma = 2.51703202e-01
  predicted kb gap to sqrt(Lambda_2) = 1.010e-02
  expected-BIC mesh refinement:
    M=16: kb=3.133273351258, sigma_BEM=2.28478680e-01, sv_min=7.463e-07, drop=8.45e+05, mesh_change=--, odd_res=1.403e-10, Rprop=2.392e-06, interior=yes
    M=24: kb=3.133273637521, sigma_BEM=2.28474755e-01, sv_min=1.827e-07, drop=3.45e+06, mesh_change=0.002%, odd_res=1.518e-10, Rprop=2.392e-06, interior=yes
    M=32: kb=3.133273683641, sigma_BEM=2.28474122e-01, sv_min=1.474e-06, drop=4.28e+05, mesh_change=0.000%, odd_res=1.435e-10, Rprop=2.392e-06, interior=yes
    M=40: kb=3.133273744923, sigma_BEM=2.28473282e-01, sv_min=9.288e-07, drop=6.79e+05, mesh_change=0.000%, odd_res=1.286e-10, Rprop=2.392e-06, interior=yes
    M=48: kb=3.133273751158, sigma_BEM=2.28473196e-01, sv_min=6.664e-07, drop=9.46e+05, mesh_change=0.000%, odd_res=1.342e-10, Rprop=2.392e-06, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141590702613, sv_min=9.636e-01
      kb(M=32)=3.141591166236, sv_min=9.637e-01, drop=1.00e+00, shift=4.636e-07, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592567655, sv_min=9.636e-01
      kb(M=32)=3.141592567655, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592613929, sv_min=9.636e-01
      kb(M=32)=3.141592613929, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592650636, sv_min=9.636e-01
      kb(M=32)=3.141592650636, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

  --- validation result ---
  symmetry conditions (a=0, X even, Y odd, nu~0): PASS
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.131493237939
  kb BEM        = 3.133273751158
  sigma asym    = 2.51703202e-01
  sigma BEM     = 2.28473196e-01
  final sv_min(A) = 6.664e-07
  final minimum drop = 9.46e+05
  final mesh change in sigma = 0.000%
  odd-parity residual = 1.342e-10
  max first propagating-mode fraction = 2.392e-06
  relative sigma error = 9.229%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.110 ===
  predicted kb = 3.119010674785
  predicted sigma = 3.76001079e-01
  predicted kb gap to sqrt(Lambda_2) = 2.258e-02
  expected-BIC mesh refinement:
    M=16: kb=3.124989777385, sigma_BEM=3.22557425e-01, sv_min=7.528e-07, drop=9.02e+05, mesh_change=--, odd_res=1.852e-10, Rprop=4.084e-06, interior=yes
    M=24: kb=3.124990573688, sigma_BEM=3.22549710e-01, sv_min=8.084e-07, drop=8.40e+05, mesh_change=0.002%, odd_res=1.642e-10, Rprop=4.083e-06, interior=yes
    M=32: kb=3.124990761648, sigma_BEM=3.22547889e-01, sv_min=7.336e-07, drop=9.26e+05, mesh_change=0.001%, odd_res=1.571e-10, Rprop=4.083e-06, interior=yes
    M=40: kb=3.124990826953, sigma_BEM=3.22547256e-01, sv_min=6.671e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.640e-10, Rprop=4.083e-06, interior=yes
    M=48: kb=3.124990855494, sigma_BEM=3.22546980e-01, sv_min=6.298e-07, drop=1.08e+06, mesh_change=0.000%, odd_res=1.597e-10, Rprop=4.083e-06, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141592567655, sv_min=9.621e-01
      kb(M=32)=3.141592567655, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592642766, sv_min=9.621e-01
      kb(M=32)=3.141592642766, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592650636, sv_min=9.621e-01
      kb(M=32)=3.141592650636, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592652784, sv_min=9.621e-01
      kb(M=32)=3.141592652784, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

  --- validation result ---
  symmetry conditions (a=0, X even, Y odd, nu~0): PASS
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.119010674785
  kb BEM        = 3.124990855494
  sigma asym    = 3.76001079e-01
  sigma BEM     = 3.22546980e-01
  final sv_min(A) = 6.298e-07
  final minimum drop = 1.08e+06
  final mesh change in sigma = 0.000%
  odd-parity residual = 1.597e-10
  max first propagating-mode fraction = 4.083e-06
  relative sigma error = 14.216%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 4/4
Leading asymptotic sigma within 5.0%: 1/4
  epsilon=0.050: BIC=PASS, asymptotic=PASS, sv_min=5.611e-07, mesh_change=0.000%, odd_res=1.255e-10, Rprop=5.806e-07, sigma_error=2.283%
  epsilon=0.070: BIC=PASS, asymptotic=FAIL, sv_min=1.811e-06, mesh_change=0.000%, odd_res=1.306e-10, Rprop=1.273e-06, sigma_error=5.188%
  epsilon=0.090: BIC=PASS, asymptotic=FAIL, sv_min=6.664e-07, mesh_change=0.000%, odd_res=1.342e-10, Rprop=2.392e-06, sigma_error=9.229%
  epsilon=0.110: BIC=PASS, asymptotic=FAIL, sv_min=6.298e-07, mesh_change=0.000%, odd_res=1.597e-10, Rprop=4.083e-06, sigma_error=14.216%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_3_x_symmetry_validation
