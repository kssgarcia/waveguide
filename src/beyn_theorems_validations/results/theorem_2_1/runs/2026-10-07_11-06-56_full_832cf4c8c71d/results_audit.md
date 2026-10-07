# Auditoría de resultados — ejecución 2026-10-07 11:06:56

El log theorem_1_beyn_1.md corresponde a esta carpeta; el manifiesto indica full-validation-completed. Se revisaron CSV y checksums de la rejilla y los diez checkpoints. No se repitió el cálculo ni se modificaron las salidas originales.

## Baseline

10/10 casos apoyan existencia y unicidad numérica dentro del contorno, con M=24. En todos, rango por umbral=gap=1, estabilidad y contorno más cercano aprobados; los cambios de S0/S1 están por debajo de 1e-3. El mayor cambio es 5.009e-4 en epsilon=.16. Residual normalizado máximo 6.707e-10 (<1e-8). Los sistemas del contorno cumplen su tolerancia.

Los diez modos pasan paredes, BIE fuera de nodos y decaimiento. Máximos: pared 1.017e-9; BIE 7.112e-6; error relativo de decaimiento 5.646e-4 (0.0565%). Estos checks físicos independientes están guardados para baseline.

## Rejilla de geometrías

99 filas: 79 one, 20 ambiguous; sin error ni geometría inválida. Los 79 tienen rango uno estable, check de contorno más cercano, diagnóstico completo y flags coherentes. Residual máximo 7.784e-10.

Límite de trazabilidad: el CSV de la rejilla no conserva cambios numéricos de S0/S1 ni ambos estimadores de rango; el código los utiliza para activar unicidad, pero no se puede auditar independientemente su convergencia por fila desde este CSV. Tampoco la rejilla ejecuta los checks físicos completos de baseline: evaluate_paper_sweep_geometry desactiva run_physical_diagnostics. Por ello, los 79 one son certificados BEM/Beyn declarados por el programa, no 79 validaciones físicas independientes equivalentes al baseline.

## Ausencia y altura crítica

Los 20 subcríticos a=.35 permanecen ambiguous, con diagnóstico incompleto y nonexistence_supported=False. No encontrar un modo no establece ausencia. Los cuatro epsilon de altura crítica tienen cota superior de existencia .511826552, pero ninguna cota inferior; no se establece a*(epsilon).

La etiqueta partial para ausencia del reporte automático no equivale a apoyo parcial certificado: cero casos certifican ausencia. La descripción legacy de zero-supported en el log no refleja por sí sola el requisito actual de diagnóstico convergente; además no hay filas zero-supported en esta ejecución.

## Asintótica y sensibilidad

En epsilon<=.1, error relativo del término principal: .0966%, .6188%, 1.6851%, 3.3331%, 5.5596%. El umbral arbitrario de 5% falla en .1; eso no contradice el teorema de epsilon suficientemente pequeño. Los datos son consistentes con la expansión, no prueba del Big-O.

En epsilon=.1, variaciones relativas máximas de sigma: paso de diferencias finitas 1.113e-7; harmonic_order 5.957e-9; lattice_terms 4.268e-9. Es estabilidad de implementación a M fijo, no convergencia de malla. La analiticidad sigue siendo una afirmación teórica.

## Conclusión

La ejecución respalda existencia y unicidad numérica de los diez baseline; reporta 79 certificados de rejilla coherentes con sus campos guardados. No valida completamente ausencia subcrítica, altura crítica ni analiticidad. Para fortalecer la auditabilidad futura, guardar S0/S1 y estimadores de rango por geometría y distinguir checks físicos realizados de no realizados.

El código fue reorganizado después de terminar esta ejecución para usar execution_harness.py; su fingerprint histórico no debe compararse como si correspondiera al fuente actual. Los checksums permiten verificar integridad de las salidas, pero no reconstruyen un snapshot del fuente ejecutado.
