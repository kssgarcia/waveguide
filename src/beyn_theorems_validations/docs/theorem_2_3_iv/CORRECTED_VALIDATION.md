# Validación corregida de Y, teorema 2.3(iv)

Versión `2026-10-13-Y-signed-channel-certified-v3`. Se conservan M=24, geometría, paredes Dirichlet, obstáculo Neumann, Green y fórmula (2.7). `embedded_validation/` y `embedded_validation_v2/` son resultados históricos/pre-fix; cada ejecución nueva se guarda en `runs/` automáticamente.

## Búsqueda y certificado

Y requiere ajustar frecuencia y desplazamiento; la simetría en x no elimina el canal propagante. Se sigue el autovector de cada semilla por overlap y se fija su fase. Las dos ecuaciones de búsqueda son el cruce espectral real y el coeficiente firmado del primer canal para la parte real de esa densidad. Se evita usar la parte imaginaria diminuta del autovalor como sustituto de cancelación de radiación.

El coeficiente modal del potencial de doble capa se deriva directamente del Green del strip:

\[
C_R=-\frac{i\,\Delta t}{2b p_1}
\sum_j u_j e^{-ip_1\xi_j}
\bigl(\xi'_j\phi'_1(\eta_j)+ip_1\eta'_j\phi_1(\eta_j)\bigr),
\]

con `phi1=cos(pi*y/(2b))` y `p1=sqrt(k^2-Lambda1)`; a la izquierda cambian ambos signos de `ip1`. El campo saliente es `C_side*exp(+i*p1*distance)*phi1(y)`: el signo de fase se fija, no se elige por mejor ajuste.

La búsqueda no sustituye la ecuación física por estas dos condiciones: el **certificado final independiente exige el operador complejo completo**, residual normalizado del vector singular, residual equilibrado por filas y mínimo singular en escala controlada (umbrales nuevos 1e-8), mínimo interior y caída, finitud, éxito del optimizador y Jacobiano condicionado. Después exige todos los checks del campo, incluidos ambos coeficientes COMPLEJOS de radiación. Un punto que sólo satisface las ecuaciones de búsqueda no se certifica.

La expansión del canal firmado es consistente con (2.7): para la densidad principal `L0(a1+Y)=a1/2+(Y-Psi)/2`, su término de orden epsilon² se anula precisamente al usar (5.47). La regresión comprueba la reducción del resto al disminuir epsilon, además de la comparación de campo y coeficiente directo.

## Estados de resultados

- `root_success` significa éxito real del optimizador; residual, condicionamiento e intentos se guardan por separado.
- `mode_exists_supported` exige certificado espectral y físico.
- `uniqueness_screen_passed` sólo describe un screen muestreado completo y estable entre dos niveles, con todos sus valles contabilizados.
- `uniqueness_supported` permanece False/`ambiguous/not-established`: este screen no aporta un conteo certificado de toda la banda ni cubre los márgenes de los umbrales. Esto no impide certificar el modo localizado.

NaN, fallos de refinamiento, candidatos pendientes o diagnósticos ausentes bloquean aprobación. La desduplicación ocurre después del refinamiento y exige frecuencia compatible con su resolución local más acuerdo del vector. Los niveles no se cuentan como modos distintos. No se declara inexistencia por falta de candidatos.

## Green, eficiencia y controles

Ensamblaje y residual BIE usan arrays para evaluar el mismo kernel y diferencias centrales que la versión escalar; regresiones contrastan ambos. Ambas separaciones verticales se reducen al período 4b, preservando la singularidad directa y evitando restar coordenadas previamente trasladadas. El campo lejano usa el Green modal con rama saliente y cola controlada, y derivadas analíticas del mismo potencial de doble capa.

El contraste independiente reutiliza el checker modal acelerado de X, ahora admitiendo una geometría explícita y matriz completa; su comportamiento anterior de X se conserva. Se compara la matriz Y y su vector con 12000/24000 términos y con harmonic_order=20/40. Se varían cuadratura del campo, secciones y sobremuestreo BIE. El orden de este diagnóstico no cambia M del problema espectral.

La sensibilidad interna reajusta **a y sigma** y repite los certificados físicos para cada configuración; un caso fallido queda explícitamente ambiguo.

## MFS y discrepancia del artículo

Además de calibración circular, convergencia de órdenes y robustez de fuentes/cuadratura, se comprueba el Neumann exterior fuera de colocación y la identidad de borde de mu. Falta de convergencia o certificados insuficientes detienen la validación.

El ejemplo 5.4 imprime `f0=-cos(t)` y `f1=.5*cos(2t)-.25`. Sustituidos en (5.47), dan una integral lineal `-3*pi/2` y `a1=-beta/2+O(beta²)`, incompatible con su última conclusión `-beta/12`. El MFS coincide con la sustitución. Se conserva (2.7); el valor impreso se presenta como discrepancia, sin forzar su reproducción. Los ajustes finitos de epsilon aportan consistencia con el término principal, no prueba de Big-O ni analiticidad.

## Ejecución

```sh
python -m unittest discover -s src/beyn_theorems_validations/tests
python src/beyn_theorems_validations/runners/validate_second_theorem_Y.py --stage smoke
python src/beyn_theorems_validations/runners/validate_second_theorem_Y.py --stage full
```

El runner smoke verifica `.02/.10/.20` e incluye contraste modal independiente y variantes físicas. El runner full añade el barrido y la sensibilidad completos, y también verifica los representantes. Ambos usan el harness existente y crean una carpeta nueva. Se puede seguir ejecutando directamente `second_theorem_Y.py` para el barrido estándar versionado; el runner añade los contrastes independientes posteriores.

El verificador admite `--check-run RUTA` para revisar una ejecución compatible sin repetir el sweep. Rechaza silencios de procedencia: el fingerprint debe coincidir con el código/configuración actuales.

## Resultado final verificado

La ejecución `2026-10-13_16-21-58_full_c622152bc891` terminó con **10/10 modos físicos apoyados**, 5/5 screens muestreados completos y estables, 86 tests aprobados y fingerprint contrastado con el código actual. Residual normalizado máximo 1.30e-11, canal abierto máximo 3.54e-11 y residual BIE máximo .002117 (umbral .005). Los tres representantes pasan el operador modal independiente con diferencia máxima 4.03e-10. El ajuste físico de a/epsilon recupera a1 con diferencia .0173% respecto de (2.7).

La sensibilidad conjunta certifica pasos de derivación 1e-8/1e-7/1e-6, lattice 25/50/100/200 y armónicos 10/15/20/30/40. Los pasos 1e-5/1e-4 fallan la verificación física posterior pese a su pequeño residual espectral, y harmonic_order=5 no cumple el certificado de matriz; quedan ambiguos. Los siete contrastes de calidad de referencia pasan el gate de estabilidad. Se mantienen los valores por defecto.

Ver [informe final](../../results/theorem_2_3_iv/runs/2026-10-13_16-21-58_full_c622152bc891/results_audit.md), certificados y logs en esa ejecución. **No se declara unicidad global ni prueba de Big-O/analiticidad.** Los resultados anteriores permanecen intactos.
