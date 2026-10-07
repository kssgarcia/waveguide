# Auditoría numérica del primer teorema de Beyn — 2026-10-07

## Dictamen y alcance

Se revisó `src/beyn_theorems_validations/first_theorem_beyn.py` y su dependencia local `lattice_sums.py`, el teorema 2.1 de WG_03.03.25.pdf, las fórmulas de Linton (1998), los CSV de `theorem_2_1_publication_ready` y la figura conjunta del barrido. No se ejecutó `main()`, el barrido largo ni los otros dos teoremas. No se modificó el código científico ni sus salidas originales. Se hicieron pruebas pequeñas de funciones, matrices de tamaño 24, series independientes y fallos controlados. Los resultados anteriores corresponden a una versión con refinamiento en M; el historial incluye el commit 2640df4 con niveles 16,24,32,40,48 y el posterior 7ca3b9b que los retira. Esto explica la diferencia de configuración, sin demostrar por sí solo la procedencia exacta de cada archivo.

El núcleo matemático es coherente. No puede aprobarse el programa como validador completo: existen defectos reproducidos en la clasificación y una evaluación de Green que pierde precisión innecesariamente por truncamiento armónico. Mantener M=24 es compatible con el estudio previo de convergencia; la auditoría no exige restablecer un barrido de M en cada ejecución.

## Correspondencia matemática

- Geometría: círculo de radio físico epsilon y centro (0,a); el contorno inflado tiene S=pi y mu=1. La normal ponderada (−Y',X') apunta al interior del obstáculo, exterior del dominio del fluido.
- Fórmulas (2.4) y (2.5): coeficiente, potencias de b, argumento pi*a/(2*b) y relación k²=Lambda_1−sigma² correctos. Para b=1,a=0.6: a0*=0.3918265520306073 y C=3.734417247 aproximadamente.
- Green: Linton (2.4), (2.46), (2.49), (2.53), (2.54) y combinación Dirichlet (4.5) concuerdan con la implementación. El período es 4b y la anchura del canal 2b. La regularización elimina Y0/4 y conserva la parte regular J0.
- Diagonal BEM: para el círculo la contribución geométrica es −1/(4pi); con el salto 1/2 y delta_theta=2pi/M resulta A=I−4pi*K/M. Signos y factores concuerdan con la convención (Delta+k²)G=delta.
- La representación de doble capa y la ecuación de contorno requieren excluir resonancias interiores Dirichlet espurias. Para los círculos y frecuencias discretas considerados, k*epsilon<pi*epsilon/(2b) permanece muy por debajo de la primera resonancia del disco; no se identificó ese mecanismo como explicación de las ramas observadas.
- Ramas complejas: Re(k)>0 y Re(k)<pi/(2b) fijan los canales abiertos/cerrados. Las raíces principales y log(k) son compatibles con la continuación en ese dominio; no se identificó un corte atravesando la elipse. Un smoke test de finitud por sí solo no verifica esta propiedad.
- Beyn: orientación antihoraria, peso dz/(i*Nq), momentos S0,S1 y B=U* S1 W Sigma^-1 correctos. Sistemas diagonales de control recuperaron 0,1,2 autovalores correctamente.
- Localización en sigma y minimización del valor singular absoluto son elecciones adecuadas cerca del umbral, aunque la matriz sigue recibiendo kb en doble precisión. No se elimina el límite de resolución de kb al aproximarse al cutoff.
- Reconstrucción: paredes y decaimiento son consistentes con el Green utilizado. El residuo fuera de nodos comparte representación y cuadratura con el operador; es una comprobación de consistencia, no una verificación completamente independiente de Neumann. Las paredes se satisfacen por construcción en el Green exacto.

## Evidencia independiente de Green

Referencia modal Dirichlet, derivada por separación de variables y equivalente a Linton (4.9):

    G_D(X,y,eta) = − sum_{n>=1} exp(−q_n |X|) sin(n*pi*y/d) sin(n*pi*eta/d)/(d*q_n)
    q_n = sqrt((n*pi/d)^2 − k^2), d=2b.

Para comparar matrices también se sustrajo y sumó analíticamente la cola de k=0, evitando la convergencia lenta del desarrollo modal en X=0. La referencia matricial usa diferencias finitas y truncamiento modal propios: sus residuos tampoco son certificados exactos.

1. Tabla 2 de Linton: para d=1,k=2,beta=sqrt(2),X=0,Y=0.01, el código obtiene −0.45952987947683693−0.35091308693822604i. Error absoluto respecto del valor impreso: 4.47e−11.
2. Quince combinaciones de frecuencia y pares de puntos, con frecuencias reales y complejas a ambos lados del eje real: error relativo máximo 7.30e−6 con H=20 y 3.17e−10 con H=40. Se probaron H=10,20,40.
3. Matrices en epsilon=0.02,0.10,0.20, M=24: diferencia relativa frente a la referencia modal acelerada con H=20 de 1.09e−7, 8.35e−6 y 2.31e−4; con H=40 de 4.17e−10, 8.72e−9 y 2.48e−6. Son errores de matriz, no errores de sigma.
4. En las 68 geometrías con status=one del barrido guardado se probaron tres pares distintos de puntos del contorno, incluyendo el arco próximo a la pared. Con H=20 sin reducción periódica, el error relativo de Green llegó a 0.19791084 en a=0.65,epsilon=0.297895. 31/68 geometrías superaron 1e−4 en alguno de los tres pares. Con H=40 todavía 17/68; con H=60, 7/68. No es una medición del error del autovalor ni una cota uniforme de todo el contorno.
5. Con reducción de ambas separaciones verticales a [−2b,2b), justificada por beta=0 y período 4b, el máximo fue 5.87e−12 con H=20 en los mismos 204 pares. H=40 y H=60 no mejoraron ese máximo. El código ya aplica esta idea en `make_extended_field_green`, pero no en el Green de ensamblaje.

**Cambio recomendado principal:** reducir periódicamente las separaciones en los Green de ensamblaje, tanto ordinario como regularizado, preservando el tratamiento de la singularidad directa. Comparar el operador y sus derivadas con la representación modal antes de cambiar resultados de publicación. El límite r<0.99*d solo evita salir del disco de convergencia; no garantiza precisión con H finito. Aumentar H indiscriminadamente tampoco es una solución universal: los coeficientes crecen y los Bessel decrecen, con riesgo de pérdida de precisión/overflow cerca de k=0.

## Efecto de la reducción periódica en sigma

Se hicieron nueve refinamientos locales de matrices M=24 en intervalos [0.8*sigma_guardada,1.2*sigma_guardada], sin Beyn ni barrido global. Las funciones originales se sustituyeron temporalmente en memoria para la prueba; no se editaron los fuentes.

| epsilon | sigma original H=20 | sigma reducida H=20 | cambio relativo | sigma reducida H=40 |
|---:|---:|---:|---:|---:|
| 0.02 | 0.001492325612938181 | 0.0014923235495743047 | 1.38265e−6 | 0.0014923235495743047 |
| 0.10 | 0.035269145975382936 | 0.03526797451049538 | 3.32161e−5 | 0.03526797451049538 |
| 0.20 | 0.11424529790946855 | 0.1142310236854055 | 1.24959e−4 | 0.1142310236854055 |

Todos los mínimos fueron interiores. Los cambios pequeños en estos tres casos respaldan su rama local, pero no autorizan trasladar esa precisión a los tamaños mayores del barrido.

## Defectos y límites de aceptación

### Crítico: no existencia después de fallos

En `independent_whole_band_svd_scan`, un barrido completamente NaN produce status=no-resolved-interior-mode y nonexistence_supported=True. Un barrido con un mínimo interior cuyo refinamiento lanza RuntimeError también produce True. Ambos casos se reprodujeron. Además, truncar la lista de candidatos al máximo configurado omite valles sin incorporar esa omisión al criterio de ausencia.

Recomendación: toda evaluación no finita, refinamiento fallido o candidato omitido debe impedir apoyar no existencia y producir un estado ambiguo con causa registrada. La ausencia de mínimos en una malla finita solo apoya no detección por encima del piso de resolución; no prueba ausencia continua.

### Importante: singularidad local aceptada solo por razón relativa

Se probó una familia no singular con s_min(sigma)=1+1e8*(sigma−0.02)^2 y s_max=1e10 en [0.01,0.03]. El criterio devuelve resolved=True: s_min=1, razón=1e−10, caída=10001. Esto demuestra que la combinación actual de razón relativa, mínimo interior y caída no basta para excluir un valle no singular con fuerte desbalance de escalas. No demuestra que los modos guardados sean falsos; sus valores singulares absolutos son mucho menores.

Recomendación: añadir un criterio absoluto coherente con la normalización natural A=I−2K y comprobar estabilidad en H, errores de ensamblaje y diagnósticos físicos. `near_singular_tolerance` existe pero actualmente no participa en `resolved`.

### Importante: rango repetido no equivale a conteo convergido

`run_beyn_convergence_study` puede detenerse con candidato estricto y rango repetido aunque S0 no converja. `unique_mode_verified` tampoco exige convergencia de S0 ni aceptación de los diagnósticos físicos. La comprobación del margen compara rangos, no posiciones ni residuos.

En las salidas, epsilon=0.02 tiene cambio final de S0 de 31.18%, ningún candidato estricto y semilla Aitken; epsilon=0.08 tiene cambio de 13.15%. Ambos figuran con unicidad verificada. Aitken es válido como semilla, no como conteo independiente.

Además, la regla del mayor salto puede subestimar rangos: con singular values [1,1e−5,1e−7,1e−9], threshold_rank=3 pero estimated_rank=1. Se reprodujo. Esto ilustra que una separación de residuos no demuestra que todos los valores pequeños sean ruido.

Recomendación: separar existencia local y unicidad apoyada; condicionar la segunda a convergencia de momentos y sensibilidad del conteo al umbral/rango, con estado pendiente cuando no se logra. Ninguna elipse con márgenes positivos cuenta literalmente los extremos omitidos de [0,Lambda_1).

### Importante: caché sin identidad suficiente

La caché verifica la rejilla (a,epsilon) y utiliza M en el nombre; no verifica b, H, términos de lattice, paso de diferencias, tolerancias ni versión del ensamblaje. Cambiar harmonic_order puede reutilizar resultados antiguos. Incorporar manifiesto de configuración y versión/hash; rechazar caché incompatible. No sobrescribir la evidencia anterior al generar una nueva versión.

### Interpretación de altura crítica y reporte

La división nominal en a<a0* y a>a0* no coincide necesariamente con la transición exacta a*(epsilon). El código trata candidatos encontrados bajo a0* como ambiguos, lo cual evita una falsa contradicción pero limita la completitud de la clasificación. El barrido crítico se detiene ante ambigüedad y lo registra correctamente. Un punto medio de un intervalo parcial amplio no es una estimación precisa de la transición.

La regla de apoyo asintótico solo exige datos finitos y que el error del primer extremo no supere al último. Cualquier Q finito en un conjunto finito es acotado: eso no establece por sí solo el Big-O. La ventana epsilon<=0.1 y la tolerancia 5% son decisiones de reporte, no constantes del teorema. La analiticidad no se verifica con esta muestra y se etiqueta correctamente como teórica.

## Auditoría de resultados guardados

- Los diez registros baseline satisfacen la dispersión a 4.44e−16 como máximo. Errores relativos y Q recalculados coinciden con los CSV.
- En epsilon=0.02,0.04,0.06,0.08,0.10, errores de sigma frente al término principal: 0.09649%,0.61835%,1.68397%,3.33124%,5.55670%. Son consistentes con mejora asintótica al disminuir epsilon. El último no pasa el umbral configurado 5%, sin que eso contradiga el teorema.
- Los diez diagnósticos físicos guardados usan M=48 y pasan sus tolerancias. No pueden atribuirse directamente a una ejecución del actual M=24.
- Barrido M32: 88 filas, 68 one,19 ambiguous,1 zero-supported. Entre los cuatro puntos subcríticos epsilon<=0.1: uno zero-supported y tres ambiguos. No existencia subcrítica no está establecida de forma completa.
- Altura crítica: un intervalo parcial de anchura 0.12 en epsilon=0.04; otros tres casos sin acotar. La ley a*(epsilon)=a0*+O(epsilon) no queda validada por estos datos.
- La gráfica conjunta muestra valores de kb muy próximos al umbral. Esa escala puede ocultar cambios relevantes en sigma; conviene conservar la comparación en sigma. La etiqueta inferior se superpone parcialmente a la leyenda; es un problema menor de presentación.

## Reproducción y límites

Scripts adjuntos: comprobaciones de Green, matrices, barrido de pares de puntos, reducción periódica, nueve refinamientos locales y fallos controlados. Requieren Python con numpy,scipy,sympy,matplotlib. Se usó un entorno temporal `/tmp/wg-audit-venv` con paquetes del sistema y sympy/matplotlib añadidos solo allí. Los scripts NO llaman `main()`. Las salidas JSON y TXT adjuntas conservan la evidencia. Los scripts de matrices y refinamiento local tardaron aproximadamente 16 y 63 segundos, respectivamente; no son el experimento completo.

No se afirma precisión uniforme de todos los autovalores ni prueba matemática del teorema. El resultado de la verificación es negativo para la aprobación completa del código actual, con causas concretas reproducidas. Las correcciones recomendadas requieren una nueva validación posterior; esta auditoría no las implementa ni sustituye los resultados anteriores.
