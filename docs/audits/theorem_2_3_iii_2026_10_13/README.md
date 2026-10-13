# Auditoría del Teorema 2.3(iii): simetría X

Fecha: 13 de octubre de 2026. Se interpreta «parte 3» como el apartado (iii) y se revisa `src/beyn_theorems_validations/second_theorem_X.py`. No se modificó código científico ni se ejecutó el barrido largo. Las comprobaciones pequeñas y los resultados históricos se guardan en `evidence.json`.

**Conclusión: la formulación de la reducción impar es pertinente, pero la implementación todavía no queda aprobada para una validación completa. Hay defectos concretos en los certificados y faltan controles de convergencia y alcance.** El harness común aporta organización, fingerprints y manejo de ejecución; no incorpora por sí mismo las correcciones científicas del Teorema 2.1.

## Qué exige el artículo y qué está bien planteado

En `docs/references/WG_03.03.25.pdf`, página 4, el Teorema 2.3(iii) exige a=0, X par, Y impar, y afirma para epsilon suficientemente pequeño un único autovalor en [Lambda1,Lambda2), con k²=Lambda2−sigma² y sigma=epsilon²*pi³*mu/b³+O(epsilon³*ln epsilon). La Observación 2.5, página 5, justifica restringir a funciones impares en y. El problema tiene paredes Dirichlet y obstáculo Neumann; el Green Dirichlet de la guía no convierte el obstáculo en Dirichlet.

La parametrización implementada cumple esa simetría y su escala epsilon está separada del obstáculo inflado utilizado para mu. Las matrices Q de elevación y R de reflexión son correctas: Q*Q=I a 7.7e-16 y RQ=−Q exactamente en el control. Para kb=2.8, epsilon=.20, M=24, el commutador relativo AR−RA fue 1.26e-10 y la parte imaginaria relativa del operador reducido 1.50e-11. Son controles puntuales, no cotas uniformes.

Linton, `docs/references/linton1998greens.pdf`, página PDF 16 (página impresa 392), ecuaciones (4.5) y (4.9), respalda la diferencia directo−imagen y la representación modal Dirichlet. La estructura y los signos del Green implementado corresponden a esa construcción. Los pesos del contorno y la matriz reducida de Beyn siguen la construcción habitual de momentos; los defectos encontrados son de resolución y certificación, no un motivo para cambiar el problema físico.

## 1. Certificación local insuficiente — crítico

`run_candidate_refinement`, líneas 1406–1409, sólo exige mínimo interior, razón smin/smax pequeña y caída. Un control inyectado con smin=1, smax=1e10, razón 1e-10 y caída 1e4 devuelve `resolved=True`. Ese valle no es singular en escala absoluta.

Falta calcular explícitamente ||Aodd*v||/(||Aodd||*||v||), exigir su tolerancia y añadir certificado de escala controlada por equilibrado, manteniendo mínimo interior y caída. Para un singular vector exacto, el residual normalizado coincide con la razón singular: no constituye una prueba independiente. Por eso se necesita además el control de escala y sensibilidad. Verificar también el vector elevado en la matriz completa: ||Afull*Qv||/(||Afull||*||Qv||) y la fuga al sector complementario.

## 2. El optimizador fallido puede aceptarse — crítico

`refine_candidate_for_M`, líneas 1312–1320, usa `result.x` sin comprobar `result.success`. El control inyectado con `success=False` y los umbrales por defecto termina con `resolved=True` cuando se suministran métricas pequeñas. No implica que un fallo así ocurriera en las salidas históricas; demuestra que la ruta de aceptación no lo impide.

Falta rechazar explícitamente optimización fallida, no finitos en el grid/bracket y diagnóstico incompleto. Los errores de brackets actualmente ignorados con `except ValueError: continue` en `validate_epsilon` deben registrarse y bloquear el conteo completo. No convertir un candidato omitido en evidencia de unicidad o ausencia.

## 3. Unicidad por rango repetido — crítico

`run_beyn_convergence_study` detiene la escalada cuando hay candidato estricto y rango repetido, sin exigir convergencia de momentos. Sólo guarda cambios de S0; no expone ni certifica S1. El tighter contour compara sólo el rango. `one_mode_supported` no exige acuerdo threshold/gap, candidato estable o residual de los solves dentro de una tolerancia.

En los resultados históricos `results/theorem_2_3_iii/embedded_validation`, 9/10 casos tienen threshold_rank=2 y gap_rank=estimated_rank=1 y aun así declaran unicidad. El caso epsilon=.02 sí termina con ambos rangos 1; tampoco tiene la verificación de S1 necesaria. En epsilon=.04, S0 todavía cambia 1.2541%; con el criterio 1e-3 usado en el Teorema 2.1 ese nivel no estaría convergido. La discrepancia 2/1 no demuestra que exista un segundo modo: demuestra que no hay un conteo consensuado.

Falta trasladar los criterios científicos del Teorema 2.1: S0 y S1 convergentes, rank estable y estimadores concordantes, candidato estable, solves controlados y tighter contour con candidato/residuales. Adaptar cuadratura al corte; la regla uniforme actual puede resolver mal la capa angular cerca de Lambda2. Aitken sólo puede suministrar una semilla. La analiticidad del operador reducido debe comprobarse con equivalencia al Green del semicanal, conjugación y continuidad en los contornos usados; finitud compleja sola no basta.

## 4. Separación de existencia, BIC físico y unicidad — crítico

`unique_embedded_mode_supported`, líneas 1876–1880, no incorpora `physical_verified` ni exige el commutador de paridad pequeño. `build_publication_claims` cuenta los checks físicos pero no los incluye en la condición de su afirmación de modo atrapado.

Un control con physical=False, cambio S0=.8, rangos 2/1 y commutador=.3 devuelve unique=True. Esto verifica un defecto de lógica, no falsedad de los modos históricos.

Faltan `mode_exists_supported`, `physical_embedded_mode_verified` y `uniqueness_supported` separados y coherentes. La cancelación del canal abierto y la paridad son esperadas por construcción del subespacio; constituyen controles de consistencia, no un segundo método independiente para probar existencia. Exigir campo no trivial, residual BIE, paredes, decaimiento y canal abierto controlados antes de afirmar un BIC físico.

## 5. Green del ensamblaje y del campo — corrección necesaria

El ensamblaje usa `greens_dirichlet`/`greens_dirichlet_reg` sin reducir periódicamente y−eta e y+eta. El campo sí reduce ambas separaciones. La periodicidad exacta no garantiza igualdad de dos series truncadas evaluadas en representantes diferentes.

Se compararon frecuencias reales y complejas, H=20/40 y tres pares de puntos con el modal de Linton. Para los dos pares cercanos al obstáculo pequeño, los errores fueron del orden de 1e-10; no se encontró aquí un daño grande para la configuración por defecto. En un control de estrés con y=.8, eta=.75 y separación horizontal .5 —fuera del obstáculo del sweep actual— el error relativo del ensamblaje llegó a 5.71e-3 con H=20 y 3.49e-5 con H=40; el Green del campo reducido quedó por debajo de 7.48e-10. No se atribuyen esos errores de estrés a los autovalores históricos.

Falta aplicar la misma reducción exacta en Green ordinario y regularizado del ensamblaje, sin alterar la singularidad directa, y comparar derivadas/matrices reducidas con una referencia modal acelerada, especialmente X=0, puntos próximos, paredes y cortes. La reconstrucción de X todavía usa exclusivamente la serie dentro de su disco; falta la representación modal de campo lejano añadida al Teorema 2.1 para puntos que puedan salir de ese dominio. No se reprodujo una excepción de radio en el punto de campo usado en esta auditoría; se señala el límite de implementación, no un fallo observado del sweep actual.

Una referencia especialmente útil es el semicanal 0<y<b, con Dirichlet en y=0 y y=b: permite ensamblar directamente el sector impar y evitar cancelaciones numéricas del canal abierto. Comparar ese operador con Q*Afull*Q antes de considerar sustituir una representación.

## 6. Qué significa unicidad en todo el apartado (iii)

El conteo actual sólo examina el sector impar en y y una elipse que omite márgenes alrededor de Lambda1 y Lambda2. Puede apoyar un modo impar único dentro de esa región, no por sí solo ausencia de otros BIC en el sector par ni en los extremos omitidos. La primera banda propagante permanece en el sector par: sus ceros del operador saliente pueden ser resonancias, no autovalores atrapados. No basta buscar ceros complejos y contarlos como modos físicos.

Falta declarar el alcance del conteo explícitamente. Para apoyar unicidad en el problema completo, combinar exclusión del sector complementario con controles de radiación y resolución del umbral, o indicar qué parte de esa exclusión se toma del resultado teórico (i) y no se valida numéricamente. Un conteo completo del sector impar mediante determinante puede reforzar el probing, pero no resuelve el sector par. Ninguna malla finita prueba ausencia arbitrariamente cerca del corte.

## 7. mu y asintótica — mejorar la evidencia

El MFS usa la normal hacia el obstáculo, restricción de carga total cero y coeficiente dipolar consistente con (2.1)–(2.3). La calibración circular mu=1 y nu≈0 es adecuada. Controles nuevos con source_scale=.5,.65,.8 y orden 160 dan mu entre 1.0797777206648569 y 1.0797777206648989 (dispersión relativa 3.89e-14): buena evidencia local de estabilidad. No se detectó un signo incorrecto de mu.

Sin embargo, la convergencia débil sólo imprime una advertencia y no bloquea la afirmación asintótica; no se exige residual Neumann off-grid ni se contrasta mu con la identidad de frontera (2.3). Añadir esos controles y variar rcond/source_scale, sin confundir el orden MFS con M del BEM.

Los datos históricos satisfacen k²+sigma²=Lambda2 a ~1.8e-15 y mejoran hacia epsilon pequeño. El error relativo del término principal es 0.259%, 1.612%, 4.306%, 8.314%, 13.440% para epsilon=.02,.04,.06,.08,.10. Esto es compatible con una aproximación asintótica que mejora al disminuir epsilon, no con igualdad exacta del término principal.

La regla de reporte actual sólo exige finitud de Q=|sigma−C*epsilon²|/(epsilon³*|log epsilon|) y compara los errores de los extremos. Cualquier conjunto finito de Q finitos es acotado; esa regla no valida por sí sola un Big-O. Falta una secuencia decreciente de epsilon con incertidumbre numérica menor que el resto, estabilidad del coeficiente C y análisis de tendencia/escala. No forzar un coeficiente logarítmico no nulo: el resto puede ser más pequeño que el bound. La analiticidad continúa siendo teórica.

## 8. Resolución espacial y regresiones

El estudio interno mantiene M=24 y cambia h, lattice y H en epsilon=.10. El dato H20/H40 es estable en sigma (~1.88e-11 relativo), pero no es una validación de todo el ensamblaje ni de la resolución de la frontera para epsilon→0. Con M fijo se puede estar viendo un coeficiente del operador discretizado. Falta una comparación independiente de alta precisión o una verificación de resolución espacial puntual, conservando M=24 para producción. No es necesario convertir M en eje del barrido largo.

Las pruebas de X/Y existentes son de harness, fingerprint y geometría; no cubren estos certificados científicos. Añadir regresiones para Green/Linton en la banda embebida, Q/R, semicanal, mu del círculo/identidad de frontera, NaN/failure/brackets, valle no singular, rank estable sin momentos, checks físicos fallidos y conteos 0/1/2. Que los tests del Teorema 2.1 pasen no certifica el código X.

## Secuencia recomendada

1. Corregir los puntos 1–5 y añadir tests científicos pequeños. Versionar implementación y conservar salidas anteriores.
2. Validar Green y operador impar frente al semicanal/modal; comprobar mu y los certificados con fallos controlados.
3. Recalcular epsilon=.02,.10,.20 con M=24: existencia, física y unicidad separadas, todos los niveles/candidatos registrados.
4. Cuando esos controles pasen, ejecutar el sweep nuevo y evaluar el régimen asintótico; ampliar conteo/exclusión si se quiere afirmar unicidad del problema completo.

No se requiere rehacer la formulación física ni cambiar M para reparar la lógica. Lo pendiente es convertir evidencia plausible de una rama impar en certificados numéricos coherentes y delimitar las afirmaciones del informe.

## Reproducir esta auditoría

Desde la raíz, con las dependencias del repositorio:

```sh
python docs/audits/theorem_2_3_iii_2026_10_13/audit_checks.py
python docs/audits/theorem_2_3_iii_2026_10_13/audit_acceptance.py
```

Sólo controles pequeños; no llaman main ni modifican resultados científicos. Los mocks son contraejemplos de lógica de aceptación, no observaciones sobre el error real de los modos guardados. `evidence.json` conserva los números del Green, proyección, MFS y filas históricas examinadas.
