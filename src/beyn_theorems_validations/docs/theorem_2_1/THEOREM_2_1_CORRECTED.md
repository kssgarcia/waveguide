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

Ver [README](../../README.md) para los comandos y estructura actuales. Todos los resultados anteriores se eliminaron a petición del usuario; los informes históricos no son evidencia disponible en disco. El cambio organizativo no ejecuta un nuevo barrido ni implica que las conclusiones del teorema estén establecidas.

## Auditoría por geometría (assembly v11)

La rejilla conserva en su CSV los rangos por umbral y gap, cambios S0/S1, convergencia de momentos, estabilidad del candidato, residual del sistema del contorno y checks físicos ejecutados/aprobados con sus residuales. La caché CSV conserva esos campos al recargar.

Cuando `run_physical_diagnostics=True` (configuración por defecto), los modos de la rejilla reciben los mismos checks físicos que baseline. No se certifica existencia física si fallan.

Cada geometría deja bajo `geometry_diagnostics/a_<a>_epsilon_<epsilon>/` los niveles de Beyn, scans independientes, candidatos evaluados (incluidos los no resueltos), errores/candidatos sin tarea, diagnósticos físicos y el certificado final. Estos archivos también se generan para los probes de altura crítica; son salidas diagnósticas, no una caché de certificados.

Un scan independiente incompleto se repite a 2 y 4 veces la densidad original, manteniendo M=24 y tolerancias. Se conservan los estados anteriores: una evaluación fallida o un candidato pendiente no desaparece al refinar. Por tanto, el refinamiento aporta evidencia y puede localizar modos omitidos, pero no garantiza establecer ausencia o altura crítica. No demuestra ausencia arbitrariamente cerca del corte.

La versión/fingerprint cambia; las salidas anteriores no se reutilizan como evidencia v11. Suite conjunta: 40 pruebas aprobadas, incluida persistencia de todos los nuevos campos. No se ha ejecutado un barrido largo de esta versión; las conclusiones científicas anteriores permanecen las de su ejecución original. La rejilla ahora tarda más por los checks físicos y los refinamientos de diagnósticos incompletos.
