# Certificación numérica de Y

Estado: completado; histórico. Responsable: Codex. Última revisión: 2026-10-13.
Alcance: implementar los hallazgos de la auditoría 2.3(iv), verificar física y matemática a resolución fija. No incluye despliegues ni prueba analítica del teorema.

Se reutilizan `second_theorem_Y.py`, `execution_harness.py` y las regresiones científicas existentes de X/2.1. Se mantienen M=24, geometría, Green de Linton, ecuación integral y umbrales físicos. Se preservan salidas históricas y se invalida caché por versión/código/configuración.

## Registro de descubrimientos

| Hallazgo | Evidencia | Fase |
| --- | --- | --- |
| NaN y candidatos pendientes pueden aprobar el scan | auditoría Y, `whole_band_uniqueness_screen` | P1 |
| root_success ignora estado del optimizador | `_solve_complex_eigenvalue_root` | P1 |
| X aporta ensamblaje por arrays, certificados y modal lejano | `second_theorem_X.py` | P2 |
| Y requiere cancelar canal propagante además del residual espectral | tablas físicas históricas, Remark 2.5 | P2 |
| sensibilidad necesita reajustar a y sigma | `run_internal_convergence_study` | P3 |

## P1 — Estados de certificación honestos

- [x] R1: revisar dataclasses, selección/refinamiento y publicación de Y; la auditoría reproduce fallos de optimización y NaN.
- [x] Tras R1, separar existencia, resultado del screen y unicidad; `ScreenResult` conserva fallos, límites y cobertura; las nuevas tablas agregan campos, sin alterar las históricas.
- [x] Tras R1, exigir estados explícitos del optimizador y finitud; intentos registrados y excepciones inválidas no se convierten en residual normal.
- [x] Regresiones de NaN, optimizador fallido, candidatos pendientes y límite de candidatos pasan en `test_second_theorem_Y_science.py`; diff inspeccionado y `git diff --check` limpio.

Gate: ningún diagnóstico incompleto aprueba unicidad o inexistencia.

## P1.1 — Corrección de desduplicación del screen

- [x] R1.1: la primera ejecución smoke ubica el mismo modo de .10 a distancias de 3e-11/5e-11 en dos grids; la tolerancia anterior era menor que el suelo de parada de scipy.
- [x] Tras R1.1, propagar la resolución de sigma a frecuencia y exigir overlap del vector para reconocer la raíz principal; registrar el nivel para no contar dos veces sus candidatos.
- [x] Regresión de resolución/vector pasa; smoke nuevo `2026-10-13_16-16-02_smoke_f9d5aa0ddda9` termina con 3/3 modos físicos y sin duplicados; suite 86 tests aprobados.

Enmienda A-01: el primer smoke pasa los checks físicos independientes, pero revela este error de conteo. Se conservan sus resultados como historial; P1.1 corrige hacia adelante y P3 usa una ejecución nueva.

## P2 — Raíz y campo físicamente certificados

- [x] R2: revisar ensamblaje/certificado/modal de X y operadores de Y; no se aplica reducción impar de X a Y.
- [x] R3: revisar MFS y (2.7)/(5.47); la sustitución impresa da -beta/2, no -beta/12.
- [x] Tras R2, patrón por arrays, separación directa, certificado normalizado/equilibrado y Green modal lejano implementados en Y.
- [x] Tras R2, condición espectral y canal firmado implementados con seguimiento por overlap, fase fija, Jacobiano centrado y certificado COMPLEJO posterior.
- [x] Tras R3, identidad de mu, residual MFS fuera de nodos, gates y regresión del ejemplo discrepante implementados.
- [x] Green/modal, periodicidad, armonicos20/40, matriz escalar/arrays, residuos engañosos y canal independiente pasan; 84 tests en 2.739 s. Se extendió el checker modal existente preservando su comportamiento de X.

Gate: regresiones físicas/matemáticas pasan sin relajar tolerancias ni cambiar M.

## P3 — Resultados nuevos y convergencia

- [x] R4: inspeccionados publicación, gráficos, `run_internal_convergence_study`, `_run_validation` y harness; se conservan campos previos y el harness como dueño de versiones/caché.
- [x] Tras R4, reajuste conjunto verificado en 15 variantes; certificados físicos detallados e internal_sensitivity_gate.json (7 variantes de referencia aprobadas). Los 3 diagnósticos insuficientes quedan ambiguos.
- [x] Casos .02/.10/.20 versionados y certificados: matrices modales 12000/24000 difieren <5e-10, variantes físicas pasan. Desde a=0 con ventana amplia se recupera .02/.10 con diferencias <1e-12 en a (`independent_placement.json`). Gate permite el barrido completo.
- [x] Validación completa 2026-10-13_16-21-58_full_c622152bc891: 10/10 modos físicos, 5/5 screens muestreados, gate modal/placement/sensibilidad aprobado; results_audit.md distingue unicidad global no establecida y consistencia finita.
- [x] 86 tests pasan (tests.log); git diff --check limpio; hashes de tablas históricas iguales a evidence.json; harness sin cambios; checker modal extendido conserva X y pasa su regresión.

Gate: reporte trazable con resultados válidos o explícitamente ambiguos; no se confunde consistencia finita con prueba de unicidad/Big-O.

## Aceptación

- [x] Hallazgos implementados o resueltos con evidencia: unicidad global explícitamente no establecida, ejemplo aritmético documentado, certificados independientes y estados fail-closed.
- [x] Tests, casos representativos y reporte vinculados en CORRECTED_VALIDATION.md; fingerprint coincide con código actual.
- [x] Sólo cambios de código/docs y salidas nuevas; tablas históricas preservadas; sin operación destructiva.
