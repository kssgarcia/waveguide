# Revisión conjunta de resultados y gráficas

Fecha: 2026-10-13. Revisión sin ejecutar barridos ni modificar algoritmos/resultados anteriores.

## Dictamen

Las ejecuciones de existencia/localización están terminadas y sus resultados son numéricamente coherentes a M=24. No corresponde afirmar que todo el contenido de los tres apartados está certificado: X sólo apoya unicidad en el sector impar/contornos estudiados; Y no apoya unicidad global. La evidencia finita tampoco demuestra Big-O o analiticidad ni constituye una prueba del operador continuo.

Se inspeccionaron visualmente 97 PNG del baseline 2.1, 9 PNG de X y 17 PNG de Y mediante hojas de contacto, junto con tablas y rutinas de generación. Se contrastaron directamente banda espectral, relación k²+sigma²=Lambda2, sigma/epsilon² y coeficiente pi³mu/b³ en los 20 casos X/Y; error de dispersión máximo 1.78e-15. `check_tables.py` reproduce estos controles y conserva `evidence.json`.

## Resultados

- 2.1: suplemento final v15 `2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b`, 20 casos con ausencia apoyada, cuatro alturas críticas bracketed con ancho .0009375, tres representantes recalculados y físicamente certificados. El baseline de diez modos y sus figuras pertenecen a la ejecución v10 del 7 de octubre, no a un nuevo barrido completo v15. La unicidad de ese baseline conserva su procedencia histórica.
- X: v3 `2026-10-13_14-59-29_full_0a01e51d6f48`, 10 modos físicos y 10 conteos convergentes en el sector impar. Unicidad completa=False.
- Y: v3 `2026-10-13_16-21-58_full_c622152bc891`, 10 modos físicos y cinco screens muestreados estables sin otro candidato resuelto. Unicidad completa=False.

No se detectó en esta revisión una contradicción nueva entre las raíces finales, sus tablas y las curvas principales. Esto es una revisión de evidencia y presentación, no una nueva auditoría exhaustiva de todas las rutas del solver.

## Hallazgos de las figuras: pendientes reales de presentación

1. **Convergencia X/Y oculta ceros.** Los PNG `internal_convergence_*` usan eje y logarítmico. X tiene cinco cambios exactamente cero y Y seis: esos puntos desaparecen y las líneas parecen cortadas/verticales. Se debe usar symlog o anotar los ceros sin atribuirles un error positivo artificial. No es divergencia del método.
2. **Certificado físico no distinguido en sensibilidad Y.** El paso h=1e-5/1e-4 puede mostrar cambio muy pequeño de sigma y aun así fallar física; H=5 falla el certificado. Las figuras actuales no marcan ese estado. Deben distinguir puntos certificados/ambiguos y mostrar también canal abierto/residual, no sólo frecuencia.
3. **Figura Y `relative_singular_value.png` muestra la tolerancia antigua parcial 1e-4.** Es una medida combinada auxiliar correcta, pero omite los umbrales nuevos 1e-8 y el certificado equilibrado/campo. Debe indicarse que es diagnóstico auxiliar o añadirse un panel del certificado efectivo.
4. **Los screens Y no muestran el mínimo refinado.** Las curvas usan sólo la malla y llegan aproximadamente a 1e-2/1e-3, aunque la raíz refinada tiene residual cercano a 1e-11. No son valores contradictorios: la malla no cae exactamente sobre la raíz. Añadir el mínimo refinado y un zoom en sigma, marcar márgenes omitidos y cambiar “Whole-band” por “Sampled-band; uniqueness not established”.
5. **2.1 carece de figuras finales consolidadas.** La ejecución final v15 no tiene PNG. Las figuras del baseline y de la auditoría parcial conservan texto “ambiguous” para a=.35; el suplemento final ya apoya ausencia en esos 20 casos. No reutilizar ese texto como conclusión actual. Hace falta una salida nueva de figuras que combine tablas con procedencia explícita y muestre alturas críticas con barras del bracket, conservando originales.
6. **Títulos poco informativos de 2.1.** `candidate_1_sigma_epsilon_*.png` dice “sigma convergence”, pero sólo contiene un punto M=24; no es un estudio de convergencia espacial. Los 40 PNG de candidato a M fijo pueden sustituirse en la presentación por tablas/paneles de certificado; no deben apoyar convergencia en M.
7. **Mu/MFS usa offset automático.** Los ejes de X/Y despliegan desplazamientos de aproximadamente 1e-10 alrededor de mu, y pueden parecer cambios grandes si se omite el offset. Conviene mostrar explícitamente |mu_N-mu_ref| o valores con contexto.
8. **Referencia discrepante a1.** El PNG Y de convergencia traza -beta/12 como “Example 5.4” sin advertir la incompatibilidad con sus fórmulas previas. El informe sí la explica. Añadir “printed example, inconsistent reference; not acceptance target”. Para beta=.5, tampoco es una referencia exacta aun si el coeficiente impreso fuera correcto.

No se observaron figuras vacías inesperadas, ejes intercambiados o curvas de frecuencia que incumplan la banda. La figura vacía de a=.35 del baseline es intencional, pero su estado está desactualizado respecto del suplemento final.

## Qué parece raro pero tiene sentido

- La aproximación asintótica y BEM se separan en epsilon grande: la ley es para epsilon→0. El error de sigma frente al término principal es .259% (X)/.201% (Y) en .02 y ~45%/~44% en .20. Esto no es error del solver ni valida la fórmula como aproximación precisa en todo el intervalo.
- La curva a/epsilon de Y no es monótona y cruza a1 cerca de .18. Esa monotonía no es exigida por el teorema. El ajuste pequeño recupera a1 con .0173% de diferencia; las raíces son físicamente certificadas. No interpretar el cruce como una segunda rama sin evidencia.
- Restos escalados crecientes con epsilon no refutan Big-O cuando disminuyen hacia el extremo pequeño; tampoco demuestran la cota a partir de diez puntos.
- NaN en kb/sigma de filas 2.1 con ausencia apoyada es una salida sin modo, no un fallo de evaluación. Los certificados de exclusión distinguen esas situaciones.

## Cierre honesto

Se puede presentar la existencia/localización y los checks físicos como resultados numéricos terminados a resolución fija. Antes de una presentación final conviene consolidar/regenerar las figuras indicadas desde las tablas, sin repetir los barridos. Para afirmar validación completa del contenido de unicidad de X/Y falta evidencia matemática adicional de exclusión del resto del espectro; no se resuelve cambiando títulos o flags.

## Correcciones aplicadas

Los hallazgos de presentación anteriores quedaron corregidos en una salida nueva de figuras: [índice de figuras revisadas](../../../src/beyn_theorems_validations/results/reviewed_figures/2026-10-13_17-46-59/README.md). Se conservaron las imágenes originales. Los tres scripts también incorporan los cambios de sus funciones plot para futuras ejecuciones.

86 tests aprobados; todos los PNG nuevos se abren correctamente y los paneles modificados fueron inspeccionados visualmente. La comparación AST excluyendo funciones plot confirma que el código numérico de los tres scripts permanece igual a la versión fuente anterior. El fingerprint de archivo completo cambia por las funciones plot; no se falsifican ni actualizan los fingerprints de las ejecuciones históricas. El manifiesto del paquete registra las tablas fuente y sus hashes.

No se ejecutó ningún barrido, refinamiento ni reconstrucción física nueva. Los límites de unicidad global, resolución fija y Big-O permanecen explícitos.
