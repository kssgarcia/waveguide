---
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python second_theorem_Y.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iv): y-symmetric BIC ===
b = 1.0
beta = 0.1
leading a1 = -beta/12 = -8.333333333333e-03
mu = 1.007140755654e+00
nu = 1.042239055775e-17
X-odd residual  = 0.000e+00
Y-even residual = 0.000e+00
Theorem shape symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590
epsilons = (0.05, 0.07, 0.09, 0.11)
refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
  a = epsilon*a1 = -4.166666666667e-04
  placement condition: PASS (residual=0.000e+00)
  predicted kb = 3.140622485938
  predicted sigma = 7.80692123e-02
  predicted kb gap to sqrt(Lambda_2) = 9.702e-04
  expected-BIC mesh refinement:
    M=16: kb=3.140666142067, sigma_BEM=7.62927596e-02, sv_min=7.244e-07, drop=7.83e+05, mesh_change=--, y-axis parity=even(4.197e-11), Rprop=1.564e-04, interior=yes
    M=24: kb=3.140666151929, sigma_BEM=7.62923536e-02, sv_min=7.003e-07, drop=8.10e+05, mesh_change=0.001%, y-axis parity=even(3.322e-11), Rprop=1.564e-04, interior=yes
    M=32: kb=3.140666154290, sigma_BEM=7.62922564e-02, sv_min=6.950e-07, drop=8.16e+05, mesh_change=0.000%, y-axis parity=even(5.053e-11), Rprop=1.564e-04, interior=yes
    M=40: kb=3.140666155126, sigma_BEM=7.62922220e-02, sv_min=6.932e-07, drop=8.18e+05, mesh_change=0.000%, y-axis parity=even(2.793e-11), Rprop=1.564e-04, interior=yes
    M=48: kb=3.140666155494, sigma_BEM=7.62922068e-02, sv_min=6.924e-07, drop=8.19e+05, mesh_change=0.000%, y-axis parity=even(2.584e-11), Rprop=1.564e-04, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141592441525, sv_min=9.841e-01
      kb(M=32)=3.141592375422, sv_min=9.841e-01, drop=1.00e+00, shift=6.610e-08, y-axis parity=odd(8.789e-08), Rprop=9.999e-01, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592580829, sv_min=9.841e-01
      kb(M=32)=3.141592580829, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.734e-07), Rprop=9.999e-01, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592642766, sv_min=9.841e-01
      kb(M=32)=3.141592642766, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.469e-07), Rprop=9.999e-01, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592652046, sv_min=9.841e-01
      kb(M=32)=3.141592652046, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.059e-06), Rprop=9.999e-01, persistent BIC=no

  --- validation result ---
  shape symmetry (X odd, Y even, nu~0): PASS
  placement a=epsilon*a1 (leading order): PASS [a1=-8.33333333e-03, a=-4.16666667e-04]
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.140622485938
  kb BEM        = 3.140666155494
  sigma asym    = 7.80692123e-02
  sigma BEM     = 7.62922068e-02
  final sv_min(A) = 6.924e-07
  final minimum drop = 8.19e+05
  final mesh change in sigma = 0.000%
  y-axis parity diagnostic = even residual 2.584e-11
  max first propagating-mode fraction = 1.564e-04
  relative sigma error = 2.276%
  asymptotic accuracy <= 5.0%: PASS
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.070 ===
  a = epsilon*a1 = -5.833333333333e-04
  placement condition: PASS (residual=0.000e+00)
  predicted kb = 3.137864020328
  predicted sigma = 1.53015656e-01
  predicted kb gap to sqrt(Lambda_2) = 3.729e-03
  expected-BIC mesh refinement:
    M=16: kb=3.138240305549, sigma_BEM=1.45093714e-01, sv_min=3.417e-06, drop=1.73e+05, mesh_change=--, y-axis parity=even(3.234e-11), Rprop=4.759e-04, interior=yes
    M=24: kb=3.138240338969, sigma_BEM=1.45092991e-01, sv_min=3.947e-06, drop=1.50e+05, mesh_change=0.000%, y-axis parity=even(5.922e-11), Rprop=4.759e-04, interior=yes
    M=32: kb=3.138240344147, sigma_BEM=1.45092879e-01, sv_min=5.428e-06, drop=1.09e+05, mesh_change=0.000%, y-axis parity=even(1.713e-11), Rprop=4.759e-04, interior=yes
    M=40: kb=3.138240413163, sigma_BEM=1.45091386e-01, sv_min=5.385e-06, drop=1.10e+05, mesh_change=0.001%, y-axis parity=even(3.287e-11), Rprop=4.759e-04, interior=yes
    M=48: kb=3.138240413756, sigma_BEM=1.45091373e-01, sv_min=5.115e-06, drop=1.16e+05, mesh_change=0.000%, y-axis parity=even(3.469e-11), Rprop=4.759e-04, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 6 total; 5 outside expected-BIC window
    checking additional minimum 1/5
      kb(M=16)=3.141592247660, sv_min=9.755e-01
      kb(M=32)=3.141592325862, sv_min=9.755e-01, drop=1.00e+00, shift=7.820e-08, y-axis parity=odd(6.729e-08), Rprop=9.998e-01, persistent BIC=no
    checking additional minimum 2/5
      kb(M=16)=3.141592491047, sv_min=9.755e-01
      kb(M=32)=3.141592491047, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.307e-07), Rprop=9.998e-01, persistent BIC=no
    checking additional minimum 3/5
      kb(M=16)=3.141592632870, sv_min=9.755e-01
      kb(M=32)=3.141592632870, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.822e-07), Rprop=9.998e-01, persistent BIC=no
    checking additional minimum 4/5
      kb(M=16)=3.141592647935, sv_min=9.755e-01
      kb(M=32)=3.141592647935, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.950e-07), Rprop=9.998e-01, persistent BIC=no
    checking additional minimum 5/5
      kb(M=16)=3.141592652046, sv_min=9.755e-01
      kb(M=32)=3.141592652046, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.133e-06), Rprop=9.998e-01, persistent BIC=no

  --- validation result ---
  shape symmetry (X odd, Y even, nu~0): PASS
  placement a=epsilon*a1 (leading order): PASS [a1=-8.33333333e-03, a=-5.83333333e-04]
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.137864020328
  kb BEM        = 3.138240413756
  sigma asym    = 1.53015656e-01
  sigma BEM     = 1.45091373e-01
  final sv_min(A) = 5.115e-06
  final minimum drop = 1.16e+05
  final mesh change in sigma = 0.000%
  y-axis parity diagnostic = even residual 3.469e-11
  max first propagating-mode fraction = 4.759e-04
  relative sigma error = 5.179%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.090 ===
  a = epsilon*a1 = -7.500000000000e-04
  placement condition: PASS (residual=0.000e+00)
  predicted kb = 3.131393237609
  predicted sigma = 2.52944248e-01
  predicted kb gap to sqrt(Lambda_2) = 1.020e-02
  expected-BIC mesh refinement:
    M=16: kb=3.133189193309, sigma_BEM=2.29629876e-01, sv_min=6.697e-06, drop=9.36e+04, mesh_change=--, y-axis parity=even(7.267e-11), Rprop=1.147e-03, interior=yes
    M=24: kb=3.133189445989, sigma_BEM=2.29626428e-01, sv_min=6.561e-06, drop=9.56e+04, mesh_change=0.002%, y-axis parity=even(2.897e-11), Rprop=1.147e-03, interior=yes
    M=32: kb=3.133189499024, sigma_BEM=2.29625705e-01, sv_min=6.651e-06, drop=9.43e+04, mesh_change=0.000%, y-axis parity=even(2.945e-11), Rprop=1.147e-03, interior=yes
    M=40: kb=3.133189558673, sigma_BEM=2.29624891e-01, sv_min=6.657e-06, drop=9.42e+04, mesh_change=0.000%, y-axis parity=even(3.635e-11), Rprop=1.147e-03, interior=yes
    M=48: kb=3.133189564015, sigma_BEM=2.29624818e-01, sv_min=6.610e-06, drop=9.48e+04, mesh_change=0.000%, y-axis parity=even(2.801e-11), Rprop=1.147e-03, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 4 total; 3 outside expected-BIC window
    checking additional minimum 1/3
      kb(M=16)=3.141592375422, sv_min=9.697e-01
      kb(M=32)=3.141592375422, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.257e-08), Rprop=9.995e-01, persistent BIC=no
    checking additional minimum 2/3
      kb(M=16)=3.141592642766, sv_min=9.697e-01
      kb(M=32)=3.141592642766, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.442e-07), Rprop=9.995e-01, persistent BIC=no
    checking additional minimum 3/3
      kb(M=16)=3.141592652784, sv_min=9.697e-01
      kb(M=32)=3.141592652784, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.896e-06), Rprop=9.995e-01, persistent BIC=no

  --- validation result ---
  shape symmetry (X odd, Y even, nu~0): PASS
  placement a=epsilon*a1 (leading order): PASS [a1=-8.33333333e-03, a=-7.50000000e-04]
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.131393237609
  kb BEM        = 3.133189564015
  sigma asym    = 2.52944248e-01
  sigma BEM     = 2.29624818e-01
  final sv_min(A) = 6.610e-06
  final minimum drop = 9.48e+04
  final mesh change in sigma = 0.000%
  y-axis parity diagnostic = even residual 2.801e-11
  max first propagating-mode fraction = 1.147e-03
  relative sigma error = 9.219%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.110 ===
  a = epsilon*a1 = -9.166666666667e-04
  placement condition: PASS (residual=0.000e+00)
  predicted kb = 3.118786624544
  predicted sigma = 3.77854988e-01
  predicted kb gap to sqrt(Lambda_2) = 2.281e-02
  expected-BIC mesh refinement:
    M=16: kb=3.124821514981, sigma_BEM=3.24183437e-01, sv_min=1.446e-05, drop=4.67e+04, mesh_change=--, y-axis parity=even(3.585e-11), Rprop=2.408e-03, interior=yes
    M=24: kb=3.124822313310, sigma_BEM=3.24175741e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.002%, y-axis parity=even(3.300e-11), Rprop=2.407e-03, interior=yes
    M=32: kb=3.124822516708, sigma_BEM=3.24173781e-01, sv_min=1.447e-05, drop=4.67e+04, mesh_change=0.001%, y-axis parity=even(3.151e-11), Rprop=2.407e-03, interior=yes
    M=40: kb=3.124822588659, sigma_BEM=3.24173087e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.000%, y-axis parity=even(3.662e-11), Rprop=2.407e-03, interior=yes
    M=48: kb=3.124822620361, sigma_BEM=3.24172782e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.000%, y-axis parity=even(2.617e-11), Rprop=2.407e-03, interior=yes
  whole-band uniqueness screen: 79 points (M=16)
     16/79
     32/79
     48/79
     64/79
     79/79
  sampled local minima: 5 total; 4 outside expected-BIC window
    checking additional minimum 1/4
      kb(M=16)=3.141592567655, sv_min=9.682e-01
      kb(M=32)=3.141592567655, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.752e-07), Rprop=9.983e-01, persistent BIC=no
    checking additional minimum 2/4
      kb(M=16)=3.141592613929, sv_min=9.682e-01
      kb(M=32)=3.141592613929, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(3.686e-07), Rprop=9.983e-01, persistent BIC=no
    checking additional minimum 3/4
      kb(M=16)=3.141592642766, sv_min=9.682e-01
      kb(M=32)=3.141592642766, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.500e-07), Rprop=9.983e-01, persistent BIC=no
    checking additional minimum 4/4
      kb(M=16)=3.141592652046, sv_min=9.682e-01
      kb(M=32)=3.141592652046, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.278e-06), Rprop=9.983e-01, persistent BIC=no

  --- validation result ---
  shape symmetry (X odd, Y even, nu~0): PASS
  placement a=epsilon*a1 (leading order): PASS [a1=-8.33333333e-03, a=-9.16666667e-04]
  expected spectral minimum resolved: PASS
  non-radiating BIC diagnostics: PASS
  additional resolved BICs: 0
  uniqueness screen: PASS
  exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
  kb asymptotic = 3.118786624544
  kb BEM        = 3.124822620361
  sigma asym    = 3.77854988e-01
  sigma BEM     = 3.24172782e-01
  final sv_min(A) = 1.448e-05
  final minimum drop = 4.66e+04
  final mesh change in sigma = 0.000%
  y-axis parity diagnostic = even residual 2.617e-11
  max first propagating-mode fraction = 2.407e-03
  relative sigma error = 14.207%
  asymptotic accuracy <= 5.0%: FAIL
  EXISTENCE / BIC / UNIQUENESS CHECK: PASS
  ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 4/4
Leading asymptotic sigma within 5.0%: 1/4
  epsilon=0.050: a=-4.166667e-04, BIC=PASS, asymptotic=PASS, sv_min=6.924e-07, mesh_change=0.000%, yparity=even(2.584e-11), Rprop=1.564e-04, sigma_error=2.276%
  epsilon=0.070: a=-5.833333e-04, BIC=PASS, asymptotic=FAIL, sv_min=5.115e-06, mesh_change=0.000%, yparity=even(3.469e-11), Rprop=4.759e-04, sigma_error=5.179%
  epsilon=0.090: a=-7.500000e-04, BIC=PASS, asymptotic=FAIL, sv_min=6.610e-06, mesh_change=0.000%, yparity=even(2.801e-11), Rprop=1.147e-03, sigma_error=9.219%
  epsilon=0.110: a=-9.166667e-04, BIC=PASS, asymptotic=FAIL, sv_min=1.448e-05, mesh_change=0.000%, yparity=even(2.617e-11), Rprop=2.407e-03, sigma_error=14.207%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_3_y_symmetry_validation
