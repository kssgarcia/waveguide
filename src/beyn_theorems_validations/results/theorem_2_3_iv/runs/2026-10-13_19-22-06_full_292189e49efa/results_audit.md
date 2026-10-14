# Resultado verificado de Y, teorema 2.3(iv)

Ejecución `2026-10-13_19-22-06_full_292189e49efa`; implementación `2026-10-13-Y-signed-channel-certified-v3`; M=24.

**Existencia física apoyada numéricamente: 10/10 casos.**
**Unicidad de la banda completa: ambiguous/not-established.** El screen finito no es un conteo certificado.

Screen refinado aprobado en 5/5 casos pequeños; no significa unicidad global.

| epsilon | a | sigma | Residual normalizado | Canal abierto | BIE off-grid | Modo apoyado |
| --- | --- | --- | --- | --- | --- | --- |
| 0.02 | -0.004595457722 | 0.01456764013 | 9.545e-12 | 1.588e-13 | 2.117e-03 | True |
| 0.04 | -0.009195163775 | 0.05754841685 | 5.298e-12 | 2.497e-13 | 5.294e-04 | True |
| 0.06 | -0.0137999325 | 0.1261222069 | 1.293e-11 | 1.290e-12 | 2.350e-04 | True |
| 0.08 | -0.01840838918 | 0.215156075 | 3.117e-12 | 1.617e-12 | 1.330e-04 | True |
| 0.1 | -0.02301722068 | 0.3178667646 | 5.275e-12 | 7.658e-12 | 8.623e-05 | True |
| 0.12 | -0.02762159687 | 0.4269467681 | 5.592e-12 | 2.698e-12 | 6.124e-05 | True |
| 0.14 | -0.03221597161 | 0.5357095478 | 4.619e-12 | 6.471e-12 | 4.661e-05 | True |
| 0.16 | -0.03679502099 | 0.6388640508 | 2.665e-12 | 4.620e-12 | 3.809e-05 | True |
| 0.18 | -0.04135435608 | 0.7327716238 | 7.791e-12 | 3.275e-11 | 4.360e-05 | True |
| 0.2 | -0.0458908098 | 0.8152898278 | 1.375e-12 | 3.537e-11 | 5.030e-05 | True |

## Contrastes independientes

Matrices modales independientes 12000/24000: diferencia lattice/modal máxima 4.025e-10; residual del vector en matriz modal máximo 1.912e-10.
Harmonic order 20/40: cambio relativo máximo 2.107e-12.
Los tres representantes .02/.10/.20 se contrastan además con cuadratura 64/128, fuentes BIE sobremuestreadas por 8 y secciones desplazadas. Estos órdenes no modifican M=24.
Búsqueda de colocación desde a=0, sin usar a1 como centro: diferencia máxima en a 8.737e-13 (.02/.10).

La incertidumbre del operador modal y el residual BIE se declaran como errores medidos, no cotas rigurosas del operador continuo. El resultado es evidencia numérica a resolución fija.

## Sensibilidad y asintótica

La sensibilidad reajusta a y sigma conjuntamente; sus tablas y `internal_physical_certificates.json` conservan también canales complejos firmados, paredes, BIE y decaimiento. Los casos de baja precisión no se convierten en inexistencia.
Los ajustes de a/epsilon y sigma/epsilon² usan modos físicos. Consistencia con el término principal no prueba Big-O ni analiticidad. El ejemplo 5.4 del artículo es internamente discrepante: su sustitución en (5.47) da -beta/2, frente al -beta/12 impreso; se conserva (2.7).

## Procedencia

Resultados corrected/v3; carpetas embedded_validation y embedded_validation_v2 son legacy/pre-fix y permanecen intactas. Cada ejecución crea una versión nueva; el fingerprint incluye configuración, código científico, Green y harness. Los contrastes rechazan una procedencia incompatible.
Ver run_manifest.json, verification_gate.json, certificates/, publication_claims.csv, independent_validation.json e independent_placement.json.
