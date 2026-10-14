# Auditoría final de resultados y figuras — Teorema 2.1

Se revisaron visualmente las 99 figuras PNG del directorio principal y se contrastaron las tablas finales, certificados de ausencia y código de clasificación. No se repitió el barrido.

Los 10 casos baseline apoyan existencia física y unicidad en el dominio buscado: momentos S0/S1 convergentes, rango estable, candidato estable, contorno más estrecho y certificados normalizados. La identidad sigma²=(pi/2)²-(kb)² concuerda con las tablas. Los 79 modos del barrido pasan los controles físicos y de unicidad configurados.

Los 20 casos sin modo tienen rango cero. En 17 hay candidatos sospechosos del scan SVD: su exclusión no se basa en ignorarlos, sino en certificados independientes de conteo por determinante (`contour_exclusion.json`), completos y con conteo cero a dos resoluciones. Los NaN de kb/sigma/residual en filas sin modo significan que no hay un modo al que asignar esas cantidades; no constituyen la evidencia de ausencia. El NaN de cambio de malla baseline corresponde a M fijo, sin estudio de malla.

Las curvas asintóticas concuerdan al disminuir epsilon: error 0.0966% en .02 frente a 5.56% en .10 y 23.53% en .20. La divergencia de curvas para obstáculos grandes no invalida una expansión de epsilon pequeño; los giros de ramas fuera de ese régimen están respaldados por los controles registrados, pero no son conclusiones del teorema asintótico. H20 y H40 dan el mismo sigma a precisión guardada en el estudio de sensibilidad.

Cuatro alturas críticas están acotadas por brackets de ancho 9.375e-4. La figura del desplazamiento normalizado usa sus puntos medios: no debe interpretarse su ligera curvatura como un coeficiente de orden superior preciso, pues el error del bracket dividido por epsilon es apreciable.

Las figuras de localización a M=24 son puntos únicos, correctamente rotulados: no demuestran convergencia de malla. Las figuras de contorno completo no resuelven visualmente la pequeña distancia entre modo y corte; las tablas dan esa separación. Los espectros de S0 ilustran rango, mientras la convergencia de los momentos completos se verifica en las tablas.

No se encontró una contradicción física o matemática en estos resultados finitos. Se corrigió una descripción inexacta del mecanismo de exclusión en el informe de publicación. El alcance sigue siendo M=24 y dominios numéricos muestreados, con resolución finita cerca del corte; no es una prueba analítica de ausencia, unicidad, Big-O o analiticidad.

## Máximos baseline

- normalized_residual: 6.70749e-10
- equilibrated_residual: 6.87997e-10
- beyn_final_s0_relative_change: 0.000500918
- beyn_final_s1_relative_change: 0.000500919
- wall_relative_residual: 7.65267e-10
- boundary_integral_relative_residual: 7.11223e-06
- decay_relative_error_left: 4.37863e-07
- decay_relative_error_right: 4.3785e-07
