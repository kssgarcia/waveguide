# Validación corregida del apartado X, Teorema 2.3(iii)

Implementación `theorem-2.3-iii-periodic-green-v3`, M=24 (dimensión impar 12). Conserva paredes Dirichlet, obstáculo Neumann, geometría, coeficiente asintótico y elipse. La versión anterior y su evidencia permanecen intactas en `results/theorem_2_3_iii/embedded_validation/`; no son resultados corregidos.

## Certificados

Un modo local requiere mínimo interior y caída, residual normalizado del vector singular, residual tras equilibrar filas y valor singular mínimo en escala controlada. También se verifica el vector elevado contra la matriz completa y la fuga de paridad. No se acepta una minimización fallida ni valores no finitos.

`mode_exists_supported` exige el certificado local, paridad del operador y todos los checks físicos: campo en las paredes y centro, paridad, cancelación del canal abierto, residual BIE off-grid y decaimiento. La cancelación del canal y la paridad se esperan por construcción; son comprobaciones de consistencia.

`odd_uniqueness_supported` exige además convergencia relativa de S0 y S1 (1e-3), estabilidad del candidato y rango, acuerdo entre threshold/gap (1), residual de solve <=1e-10 y consistencia de candidato/rangos/solve al acercar el contorno al corte. Se evalúan candidatos estrictos de todos los niveles; los duplicados resueltos no se cuentan dos veces. Un candidato sin bracket o con refinamiento fallido se registra y bloquea unicidad.

`unique_embedded_mode_supported` se conserva por compatibilidad y significa únicamente el mismo certificado de unicidad del sector impar. **`full_problem_uniqueness_supported` permanece False**: no se ha excluido el sector par ni los márgenes omitidos de los umbrales. El reporte indica ese alcance. Una búsqueda de resonancias complejas del sector par no constituye por sí sola un conteo de autovalores atrapados.

## Green, cuadratura y eficiencia

Green ordinario y regularizado reducen ambas separaciones verticales exactamente al mismo período 4b utilizado por el campo. Se preserva la singularidad directa. El ensamblaje evalúa el mismo kernel escalar y diferencias centrales de forma batched; una regresión contrasta ambos, incluyendo k complejo.

La elipse se integra con paneles Gauss concentrados en los extremos reales. Presupuestos base 256/512/1024 y adaptativos 2048/4096; no se cambia el contorno. La reconstrucción lejana usa la representación modal Dirichlet, con rama saliente para los canales abiertos y control de cola. La proyección impar se contrasta con el Green del semicanal 0<y<b y con conjugación compleja del operador reducido.

## mu y asintótica

Se conserva el MFS del obstáculo inflado y su calibración circular. Se añaden el residual Neumann en puntos que no son de colocación y la identidad integral (2.3) del artículo; falta de convergencia o controles insuficientes producen error. El orden de MFS no es M del BEM.

El informe marca la tendencia asintótica como `consistent-with-leading-term`, no como prueba de Big-O. La finitud de un conjunto de restos no establece una ley asintótica. Analiticidad permanece teórica.

## Ejecución

Desde la raíz:

```sh
python -m unittest discover -s src/beyn_theorems_validations/tests
python src/beyn_theorems_validations/runners/validate_second_theorem_X.py
```

El runner crea una carpeta nueva, verifica epsilon=.02,.10,.20 y sólo si pasan continúa el sweep completo. Reutiliza únicamente los checkpoints firmados de esos tres casos dentro de la misma ejecución; no importa certificados anteriores. La puerta exige existencia física y unicidad impar. Se conservan `representative_gate.json`, `geometry_diagnostics/`, resultados CSV, informe y manifiesto.

Una comparación independiente de matrices con Green modal acelerado, sin repetir Beyn:

```sh
python src/beyn_theorems_validations/runners/check_X_modal_operator.py RUTA_A_LA_EJECUCION
```

Esta comparación mantiene M=24 y usa truncaciones modales independientes 12000/24000. Contrasta matriz y residual del vector en la representación alternativa. Su error medido debe declararse como incertidumbre numérica; no es una cota rigurosa del operador continuo.

## Puerta de validación y refinamiento adicional

La primera ejecución corregida v2 aprobó .02 y .10 pero rechazó .20: cambio S0/S1=0.002698>0.001 y residual local=1.62e-8>1e-8. Se preserva con estado fallido; no se lanzó su sweep completo.

v3 conserva la elipse y añade paneles en la parte regular del ángulo, además de los paneles geométricos del corte. Esto resuelve modos alejados del umbral sin concentrar todo el presupuesto cerca de los extremos. La minimización bounded se repite en una coordenada centrada: evita que la regla interna sqrt(eps)*abs(x) de SciPy impida alcanzar el xatol absoluto solicitado. Se conserva el objetivo singular absoluto y se acepta el pulido sólo si lo mejora; ambas minimizaciones deben terminar correctamente. En el control local .20 el residual pasó de 1.62e-8 a 1.38e-12, sin cambiar umbrales ni M.

La referencia modal usa `expm1` y calcula la separación singular directamente antes del traslado de coordenadas; trasladar primero y restar después introduce pérdida de cifras en la diagonal. La referencia refinada, con 12000/24000 modos, concuerda con las matrices .02/.10 de v2 a 4.2e-10 y 1.4e-10 y confirma sus vectores. Las referencias iniciales de la auditoría no se interpretan como cotas del error de producción.

## Resultado verificado de v3

La ejecución `2026-10-13_14-59-29_full_0a01e51d6f48` terminó: 10/10 modos físicos, 10/10 conteos convergentes del sector impar, 69 pruebas aprobadas y fingerprint contrastado con el código actual. Residual normalizado máximo 2.06e-10, fuga de paridad máxima 8.43e-11 y canal abierto máximo 6.56e-11. Los representantes pasan además la matriz modal independiente a menos de 5e-10. Los valores sigma históricos difieren de los corregidos en como máximo 1.14e-8 relativo; eso respalda sus raíces pero no sus certificados de unicidad antiguos.

Ver [informe final](../../results/theorem_2_3_iii/runs/2026-10-13_14-59-29_full_0a01e51d6f48/results_audit.md). La unicidad del problema completo y una prueba analítica de Big-O/analiticidad no se declaran establecidas por esta ejecución.
