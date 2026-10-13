# Auditoría final de resultados corregidos — Teorema 2.3(iii), X

Versión: `theorem-2.3-iii-periodic-green-v3`. M=24. **Ejecución completa y verificada: 10/10 modos físicos y 10/10 conteos convergentes del sector impar.** El fingerprint coincide con el código actual.

69 pruebas de regresión aprobadas. La puerta .02/.10/.20 pasó antes del sweep completo. Se verificaron las tablas, certificados por caso, residuales, rango, momentos, candidatos y checks físicos contra las tolerancias guardadas.

| Métrica | Máximo en los diez casos |
|---|---|
| normalized_residual | 2.0523753e-10 |
| equilibrated_residual | 2.1127553e-10 |
| equilibrated_sigma_min | 2.3870919e-10 |
| full_vector_residual | 1.3555374e-10 |
| parity_leakage_residual | 8.4215003e-11 |
| operator_parity_commutator_residual | 3.7666014e-10 |
| beyn_final_s0_relative_change | 0.00071230013 |
| beyn_final_s1_relative_change | 0.00071230013 |
| wall_relative_residual | 3.5299284e-11 |
| centerline_relative_residual | 5.5778758e-11 |
| odd_parity_relative_residual | 1.6078194e-11 |
| first_open_channel_relative_amplitude | 6.5549037e-11 |
| boundary_integral_relative_residual | 0.0040334376 |
| dispersion_error | 1.7763568e-15 |

## Referencia modal independiente

Las matrices y vectores de .02/.10/.20 se compararon con la representación modal de Linton, acelerada por resta/suma de la cola estática y con 12000/24000 modos. Se mantiene M=24. Las diferencias de matriz y los residuales del vector son menores que 5e-10; la variación del truncamiento modal es menor que 5e-10. Estos valores son incertidumbres medidas, no cotas rigurosas del operador continuo.

## Salidas anteriores

Las salidas históricas están intactas y no se reutilizaron como certificados. Las raíces anteriores son coherentes con las corregidas a las diferencias medidas siguientes; sus antiguas afirmaciones de unicidad siguen siendo pre-fix y no adquieren validez por esta comparación.

| epsilon | sigma corregida | cambio relativo frente a legacy |
|---|---|---|
| 0.02 | 0.0133572940625 | 1.21431e-09 |
| 0.04 | 0.0527043414185 | 7.76913e-09 |
| 0.06 | 0.115337669373 | 6.09854e-09 |
| 0.08 | 0.196456318999 | 2.23951e-09 |
| 0.10 | 0.289800602731 | 9.59112e-09 |
| 0.12 | 0.388668173823 | 2.32279e-09 |
| 0.14 | 0.486917758274 | 1.86889e-09 |
| 0.16 | 0.579638784896 | 1.42185e-09 |
| 0.18 | 0.66337524402 | 6.54457e-09 |
| 0.20 | 0.735996230213 | 1.13417e-08 |

El estudio H20/H40 en epsilon=.10 muestra cambio relativo en sigma de 0. El estudio interno de h/lattice/H está en `internal_convergence.csv`; no constituye una cota uniforme de discretización espacial para todas las geometrías.

## Interpretación matemática y física

Las frecuencias cumplen k²+sigma²=Lambda2 y permanecen en [Lambda1,Lambda2). La geometría y mu satisfacen los controles de simetría, calibración circular, residual Neumann off-grid e identidad de frontera del artículo. Los vectores elevados satisfacen la matriz completa y la fuga de paridad; los campos cumplen paredes, centro, paridad, supresión del canal abierto, BIE off-grid y decaimiento con tasa sigma. Estos controles apoyan modos físicamente consistentes del problema discretizado.

**Unicidad apoyada sólo en el sector impar y en las elipses muestreadas.** `full_problem_uniqueness_supported=False` en todos los casos. El sector par y los extremos omitidos no se han excluido; por ello no se declara validado numéricamente todo el contenido de unicidad del teorema. El flag histórico `unique_embedded_mode_supported` es un alias de la unicidad impar en esta versión.

La tendencia de sigma/epsilon² es compatible con pi³*mu/b³ al disminuir epsilon. Un conjunto finito de restos no prueba un Big-O; el informe usa `consistent-with-leading-term`. La analiticidad permanece teórica. La existencia física y el conteo impar están verificados numéricamente, no demostrados para el operador continuo ni para cualquier obstáculo simétrico.

Las ejecuciones v2 rechazadas o detenidas están separadas de esta ejecución v3; se preservan para trazabilidad.
