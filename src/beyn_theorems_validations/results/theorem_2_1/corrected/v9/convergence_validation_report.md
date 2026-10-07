# Convergencia de Beyn corregida — v9

Estado: 29 regresiones aprobadas, puerta representativa aprobada y estudio dirigido completado (código de salida 0). M=24; formulación física, elipse y umbrales originales conservados.

La integración angular uniforme se reemplazó por Gauss-Legendre en paneles geométricos concentrados en los extremos reales. La estrecha capa próxima al corte ya se resuelve sin reducir tolerancias. La continuación también se aplica a búsquedas sin candidatos; agotarla no implica certificación.

| epsilon | Nq efectivo | cambio S0 | cambio S1 | residual A | unicidad |
|---:|---:|---:|---:|---:|---|
| 0.02 | 1024 | 2.151e-05 | 2.151e-05 | 8.924e-12 | supported |
| 0.04 | 1024 | 3.175e-06 | 3.175e-06 | 7.409e-12 | supported |
| 0.06 | 1024 | 2.684e-07 | 2.685e-07 | 2.729e-11 | supported |
| 0.08 | 1024 | 1.157e-06 | 1.157e-06 | 1.046e-10 | supported |
| 0.1 | 1024 | 1.457e-06 | 1.457e-06 | 2.131e-11 | supported |
| 0.12 | 1024 | 8.330e-06 | 8.330e-06 | 2.034e-10 | supported |
| 0.14 | 1024 | 9.029e-05 | 9.029e-05 | 2.227e-10 | supported |
| 0.16 | 1024 | 5.009e-04 | 5.009e-04 | 4.448e-10 | supported |
| 0.18 | 2048 | 1.345e-06 | 1.345e-06 | 6.729e-10 | supported |
| 0.2 | 2048 | 2.442e-06 | 2.442e-06 | 5.776e-10 | supported |

Los diez casos tienen rango por umbral=gap=1, estabilidad del candidato/rango, momentos convergentes, contorno más cercano consistente y checks físicos aprobados. Los contrastes H20/H40 pasan en los diez casos. Esto apoya unicidad numérica dentro del contorno y en la discretización configurada; no es una prueba del teorema ni de convergencia en M.

## Subcríticos y límites

Los cuatro casos a=.35, epsilon=.02,.04,.06,.10 tienen rango final cero estable. Sin embargo, continúan ambiguous: la búsqueda independiente atrae mínimos hacia el corte en los casos pequeños y mantiene un valle no resuelto en .10. Además, candidatos de cuadraturas gruesas se evalúan y no se borran silenciosamente. No se activa nonexistence_supported. La altura crítica no se establece por estos resultados.

## Alcance y conservación

Esta ejecución incluye los diez baseline y cuatro subcríticos. No se repitió el barrido de 99 geometrías ni el estudio completo de altura crítica con v9; sus resultados v8 permanecen intactos y no deben atribuirse a v9. La nueva regla está implementada para esas rutas, pero sus nuevas conclusiones necesitan una ejecución específica.

Los 131 archivos legacy mantienen sus hashes originales. El fingerprint de la evidencia coincide con el código científico vigente. El piloto se conserva separado y no se reutiliza como resultado final.

Datos: convergence_study/baseline_summary.csv, convergence_study/results.json y representative/epsilon_*/. Pruebas: regression_tests.log. Reproducción: validate_first_theorem_beyn.py seguido de validate_beyn_convergence.py.
