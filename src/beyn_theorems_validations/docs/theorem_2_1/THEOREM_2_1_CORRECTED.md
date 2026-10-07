# Theorem 2.1: implementación corregida v8

Se conserva el problema físico, la formulación de doble capa, el círculo, la cuadratura midpoint y las diferencias finitas centrales. El perfil de validación usa M=24. Los cambios afectan la evaluación numérica y la evidencia exigida para aceptar una conclusión.

## Green y ensamblaje

Las separaciones verticales de los términos directo e imagen se reducen al representante en [-2b,2b), utilizando la periodicidad exacta de beta=0 y período 4b. Esto se aplica al Green ordinario, al regularizado y a la reconstrucción del campo. Solo se elimina la singularidad directa en el Green regularizado. La restricción de radio de la expansión se comprueba después de reducir la separación.

El ensamblaje agrupa las mismas cuatro evaluaciones de diferencias centrales para todos los pares de nodos. La diagonal conserva el Green regularizado y el término de curvatura originales. La regresión compara este ensamblaje por lotes con la versión escalar de la misma regla Nyström.

## Certificación local y unicidad

`mode_exists_supported` y `uniqueness_supported` son independientes. El primero requiere un mínimo interior con caída suficiente, razón singular pequeña, certificado del vector singular y los checks físicos configurados. Un modo localizado conserva su valor aunque `uniqueness_status` sea `ambiguous`.

El certificado registra explícitamente

    ||A v||_2 / (||A||_2 ||v||_2).

Para un vector singular exacto, este número coincide con s_min/s_max: no constituye evidencia independiente por sí solo. Por eso se forma también B=D A, donde D normaliza cada fila no nula a norma 2 igual a uno. La norma espectral de B está acotada por sqrt(M). Se exige un residual pequeño sobre B y un s_min(B) pequeño, además del residual original. Una fila nula o datos no finitos requieren otro diagnóstico y no se certifican automáticamente. Este equilibrio solo afecta los checks; no cambia el operador físico, la localización ni sus ceros.

La unicidad exige estabilidad del rango y del candidato, convergencia de ambos momentos S0 y S1, acuerdo entre el rango por umbral y el seleccionado por gap (configurable), residuo del solve del contorno y verificación con contorno más cercano al cutoff. Agotar las cuadraturas con un modo local certificado conserva la existencia y deja la unicidad ambigua.

## No existencia

Una evaluación no finita, excepción, refinamiento fallido, valle no resuelto o candidato omitido impide apoyar no existencia. Los candidatos estrictos de todos los niveles de cuadratura se conservan para evaluación, aunque desaparezcan al cambiar de nivel. En la clasificación de geometrías, un barrido SVD limpio tampoco anula un diagnóstico incompleto de Beyn: se exige además un conteo vacío estable y momentos convergidos. Para momentos vacíos próximos a cero se utiliza el piso absoluto configurado; el cambio relativo de ruido cercano a cero no es un criterio útil.

Toda ausencia apoyada continúa siendo no detección por encima del piso de resolución configurado, no prueba matemática de ausencia hasta sigma=0.

La altura crítica se acota usando no detección en el extremo inferior y existencia certificada en el superior. Esto no afirma unicidad en el extremo superior. Se mantiene `a_upper_one_mode` como nombre de columna por compatibilidad; el informe explica que representa la cota de existencia. Los probes contienen ambos flags por separado.

## Resultados y caché

- `theorem_2_1_publication_ready/`: legacy/pre-fix; no se escribe ni reutiliza.
- `theorem_2_1_corrected_v8/`: resultados corregidos, con manifiestos y fingerprints.
- `theorem_2_1_corrected_v8_development_attempt1/`: intento previo interrumpido durante la revisión final; no es evidencia de validación terminada ni se reutiliza.

El fingerprint incorpora la configuración completa, la geometría de cada caso, versión explícita de Green/assembly y SHA-256 de ambos fuentes científicos. La caché exige además SHA-256 del CSV, versión de sus filas, M y rejilla exacta. Una caché sin manifiesto, alterada o incompatible se recalcula. La ruta legacy está protegida contra escritura desde los entrypoints corregidos.

`unique_mode_verified` se conserva únicamente como alias de compatibilidad de `uniqueness_supported`. Las ramas locales se pueden dibujar aunque la unicidad esté pendiente. El informe denomina la comparación asintótica `numerical-consistency`; una muestra finita no prueba el Big-O ni la analiticidad.

## Ejecución por etapas

Dependencias: numpy, scipy, sympy y matplotlib. Los comandos se ejecutan desde la raíz del repositorio con el Python del entorno científico.

```sh
python -m unittest discover -s src/beyn_theorems_validations/tests -p test_first_theorem_beyn.py -v
python src/beyn_theorems_validations/runners/validate_first_theorem_beyn.py
```

El segundo comando ejecuta regresiones y los tres casos epsilon=.02,.10,.20 con el perfil completo de Beyn, refinamiento, certificado y checks físicos. También refina cada modo con harmonic_order=40 y exige estabilidad frente a 20.

Para autorizar el barrido completo únicamente después de pasar esa puerta:

```sh
python src/beyn_theorems_validations/runners/validate_first_theorem_beyn.py --full-after-gate
```

La puerta exige existencia local, checks físicos, estabilidad H20/H40 y consistencia de las banderas; una unicidad honradamente ambigua no se transforma en fallo de existencia. Si falla una regresión o un caso, no arranca el barrido largo. El resultado está en `validation_gate.json`, y la inmutabilidad de legacy se verifica con SHA-256 antes y después. `run_manifest.json` distingue running/completed/failed. No debe confundirse un proceso en ejecución con una validación completada.

## Continuación v9: convergencia del conteo de Beyn

La implementación actual escribe por defecto en `theorem_2_1_corrected_v9/`. V8 conserva su barrido completo y sus conclusiones ambiguas; no se sobrescribe. `theorem_2_1_corrected_v9_quadrature_pilot/` conserva el piloto previo a extender la continuación a las búsquedas subcríticas.

Se mantiene la elipse y el operador BEM. La regla `graded-gauss` integra su parametrización con Gauss-Legendre en paneles geométricos concentrados cerca de los extremos reales. No cambia el contorno ni los umbrales; resuelve la estrecha capa angular próxima al corte que el muestreo uniforme no resolvía. La regla anterior sigue disponible como `beyn_quadrature_rule="periodic-midpoint"`. La regla y los presupuestos se incluyen automáticamente en el fingerprint.

El presupuesto de nodos se redondea al siguiente orden completo de paneles y se registra el número efectivo. Los niveles baseline son 256,512,1024 y la continuación 2048,4096. Las búsquedas de geometrías/altura crítica también continúan hasta 4096 cuando rango, momentos o candidatos no se estabilizan, incluso sin candidatos; los candidatos encontrados en cualquier nivel se conservan para evaluación. Alcanzar el presupuesto máximo no certifica nada por sí solo.

Las tres regresiones nuevas verifican pesos positivos y cobertura, integración suave y conteo/momento de un polo muy próximo al extremo. Las 26 regresiones originales permanecen.

Después de la puerta representativa vigente, el estudio dirigido se ejecuta con:

```sh
python src/beyn_theorems_validations/runners/validate_beyn_convergence.py
```

Verifica los diez baseline y cuatro subcríticos; no reemplaza ni declara repetido el barrido de 99 geometrías o el estudio de altura crítica. Exige una puerta aprobada con fingerprint vigente y conserva logs, certificados y manifiesto de finalización.
