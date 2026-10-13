# Implementación actual del Teorema 2.1

M=24; formulación física, elipse y criterios de certificación conservados.

- Green ordinario y regularizado usan la misma reducción periódica que la reconstrucción del campo.
- Un modo exige residual explícito normalizado, certificado sobre la matriz equilibrada, mínimo interior, caída y checks físicos configurados.
- Existencia y unicidad se reportan por separado. Unicidad requiere convergencia de S0/S1, estabilidad de candidato/rango, acuerdo entre estimadores y contorno más cercano.
- Fallos, NaN, candidatos pendientes y diagnósticos incompletos permanecen ambiguos; nunca apoyan no existencia.
- La cuadratura concentra nodos en los extremos reales de la elipse para resolver la capa cercana al corte. No cambia el contorno ni relaja umbrales.
- Las optimizaciones reutilizan geometría inmutable y constantes; el baseline puede usar hasta cuatro procesos.
- Los checkpoints son atómicos y exigen fingerprint y checksum. No se usan formatos pickle.

## Salidas por ejecución

Cada comando principal crea automáticamente una carpeta con fecha de Bogotá y un identificador único bajo `results/theorem_2_1/runs/`. La versión del ensamblaje se registra en el manifiesto, en vez de obligar al usuario a cambiar nombres v8/v9/v10.

El runner por etapas comparte una sola carpeta para sus pruebas, representantes y barrido completo. Una nueva invocación crea otra carpeta y no reutiliza resultados de otra ejecución.

Ver [README](../../README.md) para los comandos y estructura actuales. La limpieza inicial eliminó los resultados previos; las ejecuciones posteriores se conservan en `runs/` y deben identificarse por su manifiesto. El cambio organizativo no ejecuta un nuevo barrido ni implica que las conclusiones del teorema estén establecidas.

## Auditoría por geometría (assembly v11)

La rejilla conserva en su CSV los rangos por umbral y gap, cambios S0/S1, convergencia de momentos, estabilidad del candidato, residual del sistema del contorno y checks físicos ejecutados/aprobados con sus residuales. La caché CSV conserva esos campos al recargar.

Cuando `run_physical_diagnostics=True` (configuración por defecto), los modos de la rejilla reciben los mismos checks físicos que baseline. No se certifica existencia física si fallan.

Cada geometría deja bajo `geometry_diagnostics/a_<a>_epsilon_<epsilon>/` los niveles de Beyn, scans independientes, candidatos evaluados (incluidos los no resueltos), errores/candidatos sin tarea, diagnósticos físicos y el certificado final. Estos archivos también se generan para los probes de altura crítica; son salidas diagnósticas, no una caché de certificados.

Un scan independiente incompleto se repite a 2 y 4 veces la densidad original, manteniendo M=24 y tolerancias. Se conservan los estados anteriores: una evaluación fallida o un candidato pendiente no desaparece al refinar. Por tanto, el refinamiento aporta evidencia y puede localizar modos omitidos, pero no garantiza establecer ausencia o altura crítica. No demuestra ausencia arbitrariamente cerca del corte.

La versión/fingerprint cambia; las salidas anteriores no se reutilizan como evidencia v11. La ejecución v11 `2026-10-07_12-16-27_pending_8a0a461161b2` terminó con 76 casos `one`, 20 ambiguos y 3 errores de evaluación del campo lejano; su informe conserva esos resultados. La rejilla ahora tarda más por los checks físicos y los refinamientos de diagnósticos incompletos.

## Campo lejano (Green/assembly v12)

La reducción periódica de Y no reduce la separación horizontal X: por eso tres geometrías v11 intentaban evaluar la serie lattice fuera de `r <= 0.99*period`. Para `|X| >= b`, la reconstrucción usa la expansión modal de Dirichlet de Linton (4.9), con el mismo signo, normalización, paredes y número de onda. Se comparan las sumas de 64 y 128 modos mediante la suma absoluta de la cola; una cola relativa mayor que 1e-12 o valores no finitos provocan error. Para esta separación los modos altos decaen exponencialmente. Este número de modos de reconstrucción no es M del BEM: M sigue siendo 24.

El Green del ensamblaje y su regularización no cambian. Cambia la versión/fingerprint, y la repetición dirigida crea una ejecución nueva. La suite de 44 pruebas pasa, incluyendo expansión modal independiente cerca del corte, puntos fuera del disco lattice y coincidencia de ambas representaciones en su región común. Los resultados de la repetición se guardan en su propio `results_audit.md`.

## Extremo inferior del scan (v13)

El scan independiente incluye el extremo configurado `beyn_low_k_margin`. Antes, sus objetivos lo excluían con una desigualdad estricta, y el viaje k→sigma→k introducía redondeo adicional. Ese infinito no era un fallo del Green: era una inconsistencia entre intervalo y evaluación. Los objetivos ahora evalúan el extremo exactamente en el k configurado; fuera del intervalo y en el corte siguen rechazando. La conversión sigma→k usa directamente la raíz lejos del corte y conserva la forma racionalizada cerca de él.

No se borran fallos ni se relajan tolerancias. Los diagnósticos nuevos se calculan desde cero en otra ejecución. 47 tests aprobados, incluidos evaluación exacta del extremo, rechazo fuera del intervalo y propagación de fallos reales. Los resultados antiguos no se reinterpretan automáticamente como certificados.

## Candidatos tentativos y exclusión independiente (v14)

Un candidato generado por momentos sin converger puede desaparecer al refinar Beyn. Si su evaluación termina correctamente pero no certifica un modo, ya no se confunde automáticamente con un fallo de ejecución. Tampoco se descarta por el residual únicamente: para cerrar esa situación se requiere una comprobación independiente del determinante de la matriz completa, sin probing aleatorio.

La comprobación aplica el principio del argumento en una elipse con la misma altura y extremo inferior y un extremo superior a la resolución explícita del scan (`cutoff_gap=1e-12` por defecto). Se usa `slogdet` de la matriz con filas normalizadas positivamente, lo que preserva fase y ceros. Dos discretizaciones consecutivas deben dar el mismo conteo entero, con error de entero <=1e-6, salto máximo de fase <=pi/4, menor singular value relativo >1e-10 y residual del solve <=1e-10. Presupuestos: 256, 512, 1024, 2048. Los parámetros están en Config y fingerprint. Ningún criterio de certificación de modos se modifica.

Sólo un conteo cero resuelto puede cerrar los candidatos tentativos evaluados. Se exige además rank final cero de ambos estimadores, convergencia de momentos, scan independiente completo y sin modos, y ausencia de excepciones/candidatos sin evaluar. Un fallo real de refinamiento, NaN o solve no se borra con esta ruta. Los niveles y condiciones se guardan en `geometry_diagnostics/.../contour_exclusion.json`.

Esto es apoyo numérico a resolución finita para el operador BEM con M=24; no es una demostración analítica del problema continuo ni excluye modos en la franja entre el corte y el extremo registrado. Suite: 52 tests aprobados, incluidos los casos de conteo cero/uno, NaN, matriz mal condicionada y fase sin resolver.

## Valles del scan evaluados correctamente (v15)

La exclusión independiente también puede cerrar valles del scan que se evaluaron correctamente pero no certificaron un modo. Sus diagnósticos originales ambiguos se conservan; el cierre procede del conteo independiente de la matriz completa, nunca del mero fallo del certificado. No se permite esta ruta si hay NaN en un diagnóstico requerido, excepción, refinamiento fallido, cap de candidatos sin evaluar, falta de resultados de una tarea o momento final incompleto.

La condición del contorno usa el singular value relativo de la matriz equilibrada >1e-10, junto al residual del solve <=1e-10. Es un criterio de fiabilidad del conteo y no modifica los umbrales de certificación de un modo. Los presupuestos de candidatos aumentan de 12 a 512; el runner dirigido registra el presupuesto fuente y el actual. Tras un conteo cero resuelto se evitan los scans extra a 2x y 4x; si no se resuelve, se conservan los refinamientos y fallos anteriores. Los modos del scan, incluidos los no certificados, ahora también aparecen en `candidate_modes.csv`.

54 tests aprobados. Las ejecuciones detenidas durante la auditoría llevan `stopped-for-additional-audit` y conservan sólo evidencia parcial; la ejecución final debe consultarse por su propio manifiesto.

## Cierre de la validación dirigida

La ejecución v15 `2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b` terminó y su fingerprint coincide con el código actual. Se verificaron 20/20 exclusiones numéricas, 4/4 intervalos de altura crítica con ancho 0.0009375, los 40 probes de la bisección, tres certificados locales/físicos recalculados y tres controles del determinante con modos conocidos. El informe final distingue el alcance a resolución finita y la evidencia anterior conservada. Ver [auditoría final](../../results/theorem_2_1/runs/2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b/results_audit.md).
