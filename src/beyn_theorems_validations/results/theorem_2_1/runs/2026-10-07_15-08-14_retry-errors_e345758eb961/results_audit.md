# Targeted field-Green correction

Source run: `2026-10-07_12-16-27_pending_8a0a461161b2`. Original files preserved.
Assembly/Green version: `theorem-2.1-periodic-green-v12`; M=24.
Only prior error geometries were recomputed. No full sweep or critical study.

| a | epsilon | status | mode supported | uniqueness supported | physical passed |
|---|---|---|---|---|---|
| 0.55 | 0.449474 | one | True | True | True |
| 0.5 | 0.449474 | one | True | True | True |
| 0.5 | 0.474737 | one | True | True | True |

44 regression tests passed. All three rows were independently checked against the original configured tolerances, including rank agreement, S0/S1 convergence, tighter contour, residuals and physical checks.

| Metric | Maximum over the three cases |
|---|---|
| normalized_residual | 1.7043415e-10 |
| equilibrated_residual | 2.5571197e-10 |
| beyn_s0_relative_change | 4.8508285e-05 |
| beyn_s1_relative_change | 4.8508293e-05 |
| contour_solve_residual | 1.3837052e-14 |
| wall_relative_residual | 4.3073865e-10 |
| boundary_integral_relative_residual | 0.00013929034 |
| decay_relative_error_left | 1.2199868e-05 |
| decay_relative_error_right | 1.2199612e-05 |

The 20 ambiguous rows and critical-height conclusions in the source run are unchanged. The remaining source rows have not been re-certified as v12 by this targeted execution.
