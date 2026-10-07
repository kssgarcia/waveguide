# Verificación de la implementación corregida v8

Estado: validación completa terminada con código de salida 0; regresiones, etapa representativa y auditoría final de salidas aprobadas. run_manifest.json registra full-validation-completed.

La formulación física y matemática no cambia. M=24; harmonic_order=20, con contraste a 40 en la etapa representativa.

| Requisito | Implementación y evidencia |
|---|---|
| Green periódico ordinario y regularizado | Reducción común al período 4b en ensamblaje escalar/por lotes y campo; Linton y comparación modal aprobados. |
| No existencia ante fallos | Fail-closed para NaN, excepción, refinamiento fallido, valle no resuelto, candidato omitido/perdido y diagnóstico global incompleto. |
| Certificado del modo | Residual explícito sobre A, residual sobre B equilibrada por filas y s_min(B), más interioridad, caída y checks físicos. |
| Existencia/unicidad separadas | Flags independientes; unicidad requiere S0/S1, candidato/rango, acuerdo de estimadores y contorno más cercano. |
| Caché | Fingerprint de configuración completa, geometría de caso, versión y fuentes; checksum CSV y manifiesto obligatorios. |
| Legacy intacto | Protección contra escritura y snapshot SHA-256 de 131 archivos; verificación final aprobada: los 131 archivos conservan sus hashes originales. |
| Regresiones pequeñas | 26 tests aprobados; ver regression_tests.log. |
| Ejecución por etapas | validation_gate.json aprobado antes de iniciar el barrido completo; logs y outputs en directorio corregido. |

## Casos representativos

| epsilon | sigma | existencia | unicidad | residual A | cambio relativo H20/H40 |
|---:|---:|---|---|---:|---:|
| 0.02 | 0.00149232355577 | True | ambiguous | 1.024e-12 | 1.047e-08 |
| 0.10 | 0.0352679743791 | True | ambiguous | 5.962e-11 | 4.390e-09 |
| 0.20 | 0.114231021886 | True | ambiguous | 6.155e-10 | 9.636e-09 |

Contraste adicional de Green: 99 geometrías, 297 pares; máximo error relativo modal 1.401e-11.

Los casos pasan la puerta de implementación por existencia local, checks físicos y estabilidad armónica. La unicidad ambigua se conserva como tal: no es un fallo oculto ni una prueba de multiplicidad.

No se identifica una muestra finita con una prueba del Big-O, la analiticidad o ausencia arbitrariamente cerca de sigma=0. El reporte científico final se genera al completar el barrido.

## Resultado del barrido completo

- Baseline: 10/10 casos con existencia apoyada; residual normalizado máximo 6.155e-10.
- Rejilla de geometrías: 99 casos con M=24; 79 mode-found y 20 ambiguous, sin errores de ejecución.
- Unicidad: 0 casos certificados. No existencia: 0 casos certificados. Las conclusiones permanecen ambiguous/not-established.
- Altura crítica: los cuatro valores de epsilon quedaron unbracketed-or-ambiguous; sólo se obtuvo la cota superior de existencia.
- Legacy: SHA-256 idénticos para los 131 archivos anteriores.
- Invariantes de versiones, M, residuales y banderas aprobados; evidencia en final_output_audit.json.
- Ensamblaje por lotes frente a escalar para M=24: error relativo máximo 1.674e-11; ver batch_vs_scalar_M24.json.

El informe científico generado es publication_validation_report.md. Su etiqueta partial en las filas de unicidad y no existencia indica cobertura incompleta del estudio; con 0 casos certificados, esas afirmaciones no están establecidas por esta ejecución. La implementación corregida no equivale a una validación completa de todas las conclusiones del teorema.
