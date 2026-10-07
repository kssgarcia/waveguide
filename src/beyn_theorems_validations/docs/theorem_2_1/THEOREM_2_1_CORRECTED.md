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
