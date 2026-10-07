# Auditoría de los resultados corregidos del Teorema 2.1

Se analizaron las salidas existentes; no se repitió el barrido ni se modificó el código científico. M=24.

## Existencia local

Los diez casos baseline apoyan existencia local y pasan los checks físicos configurados. Residual normalizado máximo: 6.155e-10; equilibrado: 6.346e-10. Residual de pared máximo: 9.555e-10; residual BIE fuera de nodos máximo: 7.112e-6. Error relativo máximo de decaimiento: 5.646e-4 (0.0565%). Son certificados numéricos del problema discretizado y evidencia física; no constituyen prueba ni estudio de convergencia en M.

La rejilla contiene 99 geometrías admisibles: 79 mode-found y 20 ambiguous. Residual normalizado máximo de los 79: 1.750e-9; equilibrado: 1.799e-9, ambos inferiores a 1e-8. No hay certificados con residual no finito.

## Conteo global y unicidad

| epsilon | cambio S0 | cambio S1 | rango umbral/gap | contorno más cercano consistente |
|---:|---:|---:|---|---|
| 0.02 | 0.311802 | 0.311812 | 2/1 | False |
| 0.04 | 0.0319366 | 0.0319416 | 2/1 | True |
| 0.06 | 0.140439 | 0.140437 | 2/1 | True |
| 0.08 | 0.00635371 | 0.00635288 | 2/1 | True |
| 0.1 | 0.00177719 | 0.00177773 | 2/1 | True |
| 0.12 | 0.000708013 | 0.000708398 | 2/1 | True |
| 0.14 | 0.000325718 | 0.000326012 | 2/1 | True |
| 0.16 | 0.000169514 | 0.000169749 | 2/1 | True |
| 0.18 | 9.78071e-05 | 9.8002e-05 | 2/1 | True |
| 0.2 | 6.15448e-05 | 6.17115e-05 | 2/1 | True |

Con tolerancia de momentos 1e-3, sólo epsilon=.12,.14,.16,.18,.20 convergen en la comparación final. Los diez casos mantienen desacuerdo 2/1 de rango; por tanto ninguno certifica unicidad. epsilon=.02 además falla estabilidad del candidato y contorno más cercano; epsilon=.04 falla estabilidad del candidato. La estabilidad del rango seleccionado no basta.

El segundo valor singular de S0 desciende al aumentar Nq: en epsilon=.20 pasa de 7.588e-8 a 5.046e-10, mientras el primero se estabiliza cerca de .0121823. Esto es compatible con contaminación de cuadratura en el segundo valor singular; no demuestra un segundo modo ni permite descartarlo con el presupuesto actual. Los residuales de los sistemas del contorno están cerca de precisión máquina (hasta 2.45e-15 en los ejemplos .02,.06,.10,.20); no explican por sí solos la falta de convergencia de los momentos.

## Ausencia y altura crítica

Los 20 casos a=.35 tienen cero modos resueltos, pero diagnóstico_complete=False y nonexistence_supported=False. Esto significa ausencia no establecida, no ausencia demostrada. No se encontró una bandera de ausencia apoyada con diagnóstico incompleto. Los cuatro estudios de altura crítica sólo tienen cota superior de existencia a=.511826552; sin cota inferior, la altura crítica permanece not-established.

## Asintótica y estabilidad

Para epsilon<=.10, el error relativo frente al término principal crece de 0.0966% a 5.5596%; el resto escalado Q va de .0461 a .9017. Es coherencia numérica con la expansión en esta muestra, no prueba del Big-O ni de analiticidad. Los contrastes H20/H40 de los tres representantes son inferiores a 1.05e-8 relativos en sigma; la estabilidad armónica no resuelve el conteo de Beyn.

## Conclusión y próximo experimento

Las salidas son internamente coherentes con los certificados corregidos. Apoyan existencia local y consistencia asintótica; no establecen unicidad, ausencia subcrítica ni altura crítica. El estado partial del informe automático para unicidad/ausencia no debe leerse como certificación: hay cero casos certificados para esas afirmaciones.

Antes de otro barrido, conviene un experimento de convergencia de cuadratura en pocos casos (.02,.10,.20 y un subcrítico), registrando S0/S1, sus valores singulares, estabilidad de candidatos y contornos. Mantener M=24 y todos los umbrales actuales; no relajar el umbral de rango para obtener una conclusión deseada. El descenso observado del segundo valor singular justifica estudiar Nq, pero no garantiza que aumentar Nq baste.

Legacy: 131 archivos con hashes originales conservados. Evidencia de invariantes en final_output_audit.json.
