# Auditoría final del Teorema 2.1

Ejecución: `2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b`. Estado: **terminada y verificada**. Implementación `theorem-2.1-periodic-green-v15`; fingerprint comprobado contra los archivos actuales.

Se mantuvieron M=24, la formulación, el Green del ensamblaje y los criterios de aceptación de modos. Se corrigió el extremo inferior del scan y se añadió la exclusión independiente por el determinante de la matriz completa para candidatos evaluados correctamente que no certifican un modo. Fallos, NaN de evaluación y candidatos sin evaluar siguen impidiendo apoyar ausencia. El presupuesto de candidatos aumenta de 12 a 512 y queda registrado.

**Resultados verificados:**

- 20/20 geometrías antes ambiguas: `zero-supported`, con evidencia de conteo cero convergente y sin fallos de evaluación.
- 4/4 alturas críticas: `bracketed`, con ocho pasos de bisección y extremos respaldados por exclusión numérica y existencia certificada.
- 3/3 representantes (epsilon=0.02, 0.10, 0.20): refinamiento local y checks físicos recalculados y aprobados.
- 3/3 controles con modos conocidos: el determinante cuenta uno y rechaza exclusión, incluida la adaptación de resolución de fase.
- 54 tests de regresión aprobados.

| epsilon | cota inferior sin modo detectado | cota superior con modo certificado | ancho |
|---|---|---|---|
| 0.04 | 0.3899515520 | 0.3908890520 | 0.0009375 |
| 0.06 | 0.3890140520 | 0.3899515520 | 0.0009375 |
| 0.08 | 0.3880765520 | 0.3890140520 | 0.0009375 |
| 0.10 | 0.3871390520 | 0.3880765520 | 0.0009375 |

La columna histórica `a_upper_one_mode` representa la cota de existencia. La unicidad permanece separada: no se fuerza a partir de esa cota ni del conteo del determinante; conserva los requisitos de momentos/Beyn.

**Alcance:** apoyo numérico para el operador discretizado con M=24, k b >=0.0001 y una resolución de distancia al corte de 1e-12 en el scan/exclusión. No es una prueba analítica del problema continuo ni excluye modos dentro de la franja no explorada junto al corte. Los momentos vacíos se verifican con el umbral absoluto configurado; cambios relativos grandes cerca de cero no se usan como evidencia de unicidad.

**Procedencia:** los resultados anteriores están intactos. No se repitió la rejilla completa, ni la búsqueda Beyn de los diez baseline. Los roots anteriores sólo se usaron como semillas para los tres certificados locales/físicos nuevos. La evidencia de unicidad del baseline original pertenece a su propia ejecución y versión.

Archivos de evidencia: `nonexistence_audit.json`, `critical_height.csv`, `critical_height_probes.csv`, `representative_rechecks.json`, `regression_tests.log` y `geometry_diagnostics/`. Las ejecuciones parciales de auditoría conservan su estado detenido y no deben confundirse con ésta.
