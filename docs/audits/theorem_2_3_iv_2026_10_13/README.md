# Revisión inicial del teorema 2.3(iv): simetría Y

Fecha: 2026-10-13. Alcance: revisión del código actual, documentos de referencia y resultados históricos; controles pequeños de auditoría. No se modificó el algoritmo científico ni se ejecutó el barrido largo. Se conserva `M=24`.

Esta revisión describe el estado **legacy/pre-fix** identificado por los hashes de `evidence.json`. La implementación posterior y sus verificaciones se documentan en [validación corregida](../../../src/beyn_theorems_validations/docs/theorem_2_3_iv/CORRECTED_VALIDATION.md). El control histórico rechaza ejecutarse sobre una versión nueva para no sobrescribir silenciosamente su evidencia.

**Conclusión:** hay una base matemática y numérica útil, pero Y todavía no está listo para presentar una validación completa. Los resultados guardados respaldan un candidato que pasa los checks físicos existentes, no toda la rama ni su unicidad. Hay fallos reproducibles de certificación y una discrepancia aritmética interna en el ejemplo del artículo que debe documentarse sin forzar el código a reproducirla.

## 1. Correspondencia con el teorema

Fuentes: `docs/references/WG_03.03.25.pdf`, teorema 2.3(iv), ecuación (2.7), Remark 2.5 y sección 5.2, especialmente (5.47) y ejemplo 5.4 (página PDF 24). El Green debe contrastarse con `docs/references/linton1998greens.pdf`, como en las auditorías anteriores.

La geometría implementada es

\[
X(t)=\sin t-\frac{\beta}{2}\sin 2t,\qquad
Y(t)=-\cos t+\frac{\beta}{2}\cos 2t,
\]

con obstáculo físico `(epsilon*X, a+epsilon*Y)`. Tiene `X` impar y `Y` par. La búsqueda conjunta de frecuencia y desplazamiento `a` corresponde al caso Y: la simetría respecto del eje vertical no elimina por sí sola el primer canal propagante. Por eso **no se debe trasladar directamente la certificación del sector impar usada en X**.

Partes correctas o bien encaminadas en el código actual:

- Usa la matriz BEM completa y ajusta conjuntamente `a` y `sigma`, con continuación y varias semillas.
- Calcula `mu`, `nu` y `a1` mediante un problema exterior MFS; implementa directamente (2.7).
- La reducción periódica vertical ya aparece en Green ordinario y regularizado del ensamblaje, y en la reconstrucción del campo.
- Exige mínimo interior, caída del valle y dos medidas de singularidad, y después verifica canal abierto, paredes, simetría, ecuación integral fuera de nodos y decaimiento.
- El barrido adicional sólo se activa cuando el candidato principal pasa sus checks físicos.
- Ya usa `execution_harness.py`, directorios versionados, manifiestos, caché con fingerprint y una versión explícita de ensamblaje. Esta infraestructura debe conservarse.

Estos puntos no reemplazan una comparación independiente Green/modal ni la convergencia de los diagnósticos físicos específicos de Y.

## 2. Discrepancia del ejemplo 5.4: evidencia a favor de (2.7)

El artículo imprime

\[
f_0=-\cos t,\qquad f_1=\tfrac12\cos 2t-\tfrac14,
\qquad 2S+\pi\mu=3\pi+O(\beta),
\]

y (5.47) equivale a integrar `(Y*X' - 3*X*Y')*f` y dividir por `2S+pi*mu`. La integral de orden cero se anula. La integral del coeficiente de `beta`, al sustituir **esas mismas expresiones impresas**, es

\[
\int_{-\pi}^{\pi}\partial_\beta[(YX'-3XY')f]_{\beta=0}\,dt
=-\frac{3\pi}{2}.
\]

Por tanto resulta `a1=-beta/2+O(beta^2)`, mientras el ejemplo concluye `-beta/12+O(beta^2)`. La integración simbólica y una cuadratura trigonométrica independiente coinciden. La corrección de primer orden del denominador no puede resolver la discrepancia porque el numerador de orden cero es cero.

El MFS del script también produce `a1/beta=-0.4999998333` para `beta=.001` y `-0.4999333356` para `.02`. En puntos distintos de los de colocación, satisface la condición de Neumann con error relativo cercano a `3e-14` para esas deformaciones pequeñas; la identidad de borde para `mu` coincide cerca de precisión numérica. Además reproduce `mean(f1)=-.25` y el coeficiente de `cos(2t)` tendiendo a `.5`.

**Interpretación:** existe una incompatibilidad entre el último valor declarado en el ejemplo y sus fórmulas previas. Esto no identifica por sí solo cuál sería la corrección autoral, pero sí impide tratar `-beta/12` como referencia exacta obligatoria. Se recomienda conservar (2.7), añadir esta regresión y presentar el valor del ejemplo como referencia discrepante, no como criterio de aprobación. La equivalencia entre (2.7) y (5.47), evaluada con la misma `Psi`, es una comprobación algebraica, no una validación independiente del MFS.

## 3. Lo que realmente muestran las salidas guardadas

La carpeta `results/theorem_2_3_iv/embedded_validation/` contiene una ejecución histórica completa. No tiene un manifiesto que vincule inequívocamente todos sus resultados con el código actual; debe tratarse como evidencia histórica.

| epsilon | Candidato espectral resuelto | Checks físicos guardados | Canal abierto relativo |
| --- | --- | --- | --- |
| .02 | Sí | No | 1.42245e-5 |
| .04 | Sí | Sí | 7.78332e-7 |
| .06 | Sí | No | 1.09986e-5 |
| .08 | Sí | No | 6.96090e-5 |
| .10 | Sí | No | 1.16683e-4 |

El umbral del canal abierto es `1e-5`: sólo **1/5 casos pequeños** pasa. En los diez casos guardados sólo `.04` tiene `physical_bic_verified=True`; `.12`, `.14` y `.16` ni siquiera quedan espectralmente resueltos. No hay que interpretar los demás como inexistencia: son candidatos insuficientemente certificados o búsquedas no concluyentes.

Para los cinco casos pequeños, los residuos BIE y los errores de decaimiento guardados están dentro de sus umbrales; el fallo del canal abierto basta para impedir la certificación. No se estableció aquí si su amplitud proviene principalmente de desplazamiento mal ajustado, truncación del Green, derivación numérica o reconstrucción: eso requiere contrastes independientes.

Sólo `.04` ejecuta el barrido adicional. Su resultado significa que no encontró otro candidato no radiante resuelto en la malla muestreada, **no una exclusión certificada de toda la banda**.

`embedded_validation_v2/` es una salida parcial: contiene diagnósticos de forma/MFS y algunas imágenes de barrido, pero no el resumen final ni la tabla física completa. No sustituye a la ejecución anterior ni permite afirmar que Y ya pasó. No se borró ni alteró ninguna salida.

## 4. Hallazgos y cambios necesarios

### P1 — Unicidad puede aprobar con diagnóstico incompleto

En `whole_band_uniqueness_screen` y `validate_epsilon` (líneas 2031 y 2119 del código auditado), la aprobación se reduce a `additional_resolved == 0`.

Control reproducible: nueve valores `NaN` generan cero brackets/candidatos y la fórmula de aprobación devuelve `True`. También se descartan brackets cerca del candidato principal mediante un radio fijo `.002`, se limita a diez candidatos sin registrar los descartados, y no se exige éxito del refinamiento. Un candidato adicional que no cumple todos los umbrales no equivale automáticamente a exclusión de ese intervalo.

**Cambio necesario:** separar `mode_exists_supported` y `uniqueness_supported`; registrar cobertura, extremos omitidos, brackets evaluados/descartados, fallos y candidatos ambiguos. Cualquier intervalo pendiente, valor no finito o refinamiento fallido impide aprobar unicidad. Desduplicar raíces después del refinamiento con incertidumbre local, sin borrar anticipadamente vecindarios amplios. El barrido real necesita refinamiento/convergencia y una evidencia de conteo o exclusión compatible con el operador; no basta importar un conteo Beyn sin justificar su dominio de analiticidad y umbrales.

### P1 — El éxito del optimizador y la certificación se confunden

`_solve_complex_eigenvalue_root` (1219) acepta el punto terminal por el tamaño del autovalor sin exigir `result.success`. `root_success` (1515) expresa sólo el cumplimiento del umbral del residual. Un control con `least_squares.success=False` y autovalor pequeño reproduce esta aceptación.

**Cambio necesario:** almacenar por separado éxito/estado del optimizador, residual y certificado final, incluyendo todos los intentos. Un punto terminal fallido puede recuperarse mediante una certificación independiente explícita; no debe presentarse como optimización exitosa. Las penalizaciones `1.0` o `[1,1]` ante excepciones también necesitan un estado de fallo para no confundirse con evaluaciones válidas.

La selección del autovalor por mínimo módulo en cada evaluación (`near_zero_eigenpair`, 1141) puede cambiar de rama; la conversión de variables recorta `sigma` a la banda (1295), creando zonas planas. Deben controlarse continuidad de la rama, límites activos y condicionamiento del Jacobiano. Estas son limitaciones del algoritmo, no pruebas de que todas las raíces históricas sean falsas.

### P1 — Singularidad espectral pequeña no garantiza cancelación de radiación

El código ya usa `smin/smax` y `smin/median(row_norms)`, lo cual es mejor que un valor absoluto sin escala. Todavía falta una certificación explícita del vector singular en una matriz con escala controlada, con residual `||Av||/(||A||*||v||)`, escala de filas y controles de finitud. El residual y `smin/smax` son esencialmente equivalentes cuando ambos se calculan por la misma SVD en norma 2; su registro explícito debe complementarse con escala y verificaciones independientes del campo.

**Cambio necesario:** la aceptación de la búsqueda debe exigir también la supresión del canal abierto. Hoy puede detenerse en una raíz espectral y fallar después en física, como ocurre en cuatro de los cinco casos pequeños. Investigar si las dos ecuaciones basadas en el autovalor complejo resultan mal condicionadas cerca del BIC; contrastarlas con el coeficiente firmado del canal propagante y una fase fija. Si se reformula el residual de búsqueda, debe derivarse y comprobarse su equivalencia física, conservando la formulación del problema y sus tolerancias.

### P2 — Falta independencia y convergencia de los checks del Green/campo

`first_open_channel_diagnostics` (1700) prueba ambos signos de fase y elige el que menos varía. La convención saliente debe deducirse del Green y quedar fijada; elegir el signo a posteriori debilita la verificación. La amplitud casi nula requiere una prueba de fase adaptada a su incertidumbre, no una división inestable.

**Cambio necesario:** reproducir Linton, comparar Green lattice/modal en geometrías desplazadas y cerca de coincidencia, y certificar las matrices con `harmonic_order=20/40`, truncación y paso de derivación. Comparar la proyección del campo con un coeficiente modal obtenido directamente de la densidad. Estudiar convergencia de cuadratura transversal, distancias de observación y sobremuestreo del residual BIE, manteniendo `M=24`. Evitar pérdidas de precisión al calcular diferencias de coordenadas tras sumar el desplazamiento del strip; usar la separación directa antes de envolverla periódicamente. Considerar la representación modal para los puntos lejanos, previa comprobación de equivalencia.

### P2 — La sensibilidad actual mantiene fijo el parámetro físico que cancela radiación

`run_internal_convergence_study` (2307) cambia parámetros internos y reajusta sólo `sigma`, manteniendo `a` fijo. En Y la solución depende de ambos.

**Cambio necesario:** reajustar conjuntamente `(a,sigma)` para cada variante y registrar desplazamientos, canal abierto firmado y checks físicos. No confundir estabilidad de frecuencia con estabilidad del BIC. Convertir las advertencias de convergencia MFS/robustez de `a1` en criterios explícitos para las afirmaciones que dependen de ellas, incluyendo residual exterior fuera de nodos e identidad de `mu`.

### P2 — Asintótica y alcance de las afirmaciones

El ajuste de `a/epsilon` requiere varios puntos físicamente certificados y una exploración que no imponga de antemano la pendiente que se pretende validar. La parametrización `a=epsilon*a1+c*epsilon^2` es útil para buscar, pero la recuperación independiente de `a1` necesita verificar estabilidad frente a ampliación de ventanas/semillas y continuidad de rama. El punto histórico `.02` merece revisión por su corrección de desplazamiento claramente distinta de las siguientes.

Una lista finita de restos escalados acotados aporta consistencia asintótica; no demuestra un término `O(...)`. La falta de aprobación para epsilon moderado no contradice un teorema para epsilon suficientemente pequeño. Si sólo hay un BIC físicamente certificado, no se puede presentar como establecida la ley de toda la rama.

## 5. Orden propuesto de implementación y verificación

1. Resolver estados de fallo, aceptación del optimizador, unicidad ambigua y contabilidad completa de candidatos; añadir regresiones que reproduzcan los controles de esta auditoría.
2. Documentar la discrepancia del ejemplo 5.4 y certificar el MFS con controles independientes, sin reemplazar (2.7) por `-beta/12`.
3. Reforzar certificado de matriz/vector y validar Green/modal, fase saliente y coeficiente de canal independiente. Mantener `M=24` y las tolerancias físicas.
4. Mejorar seguimiento y ajuste conjunto de la raíz, con diagnósticos de condicionamiento y cancelación de radiación; comprobar convergencia interna reajustando ambos parámetros.
5. Ejecutar smoke/regresiones y casos `.02`, `.10`, `.20`. Los casos moderados pueden quedar ambiguos: el criterio es que el diagnóstico sea fiable, no que se fuerce un PASS.
6. Sólo después realizar una ejecución completa nueva con manifiesto y resultados versionados. Si la validación asintótica requiere epsilon menor, añadir una secuencia pequeña tras resolver la fiabilidad de búsqueda y campo.

## 6. Evidencia reproducible y límites

`audit_checks.py` realiza los controles pequeños descritos y escribe `evidence.json`, con hashes SHA-256 de código, documentos y tablas históricas. Contiene inyección de fallos para demostrar las decisiones lógicas; no prueba que aquellos fallos ocurrieron en la ejecución histórica. Se ejecutó en el entorno disponible con NumPy, SciPy y SymPy, sin iniciar el punto de entrada del experimento.

Reproducción, desde la raíz y con el entorno científico del proyecto:

```sh
PYTHONDONTWRITEBYTECODE=1 python docs/audits/theorem_2_3_iv_2026_10_13/audit_checks.py
```

Se revisó código fuente y tablas, y se contrastó simbólicamente el ejemplo del artículo. No se ejecutó una nueva reconstrucción completa de campo ni un conteo espectral independiente para Y; por ello esta revisión identifica lo necesario para cerrar la validación, pero no declara el teorema numéricamente certificado.
