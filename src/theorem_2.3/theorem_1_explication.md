# Modo atrapado en una guía de ondas: del teorema al BEM

La mejor manera de entender el código es verlo como una cadena matemática completa:

$$
\text{PDE en una guía infinita} \to \text{problema espectral} \to \text{ecuación integral} \to \text{BEM} \to A(k)u=0 \to \sigma_{\min}(A(k)) \to k_{\text{BEM}} \to \sigma_{\text{BEM}} \to \text{comparación con el teorema}
$$

Y hay una distinción que quiero hacer desde el comienzo porque es fácil confundirse:

$$
k \neq K \neq \sigma \neq \sigma_{\min}(A)
$$

Son cuatro objetos completamente diferentes.

---

## 1. El problema físico original

Tenemos una guía de ondas bidimensional infinita. Si $b$ es el semi-ancho:

$$
\Pi_b = \left\{ (x,y)\in\mathbb{R}^2 : -b<y<b \right\}
$$

La guía se extiende infinitamente en $x$:

$$
-\infty<x<\infty
$$

y sus paredes están en $y=-b$ y $y=b$.

Dentro de la guía ponemos un obstáculo circular de radio $r=\varepsilon$, centrado en $(0,a)$.

En tu experimento: $b=1$, $a=0.6$.

La geometría se ve conceptualmente así:

```
        y = b
   ------------------
          ○  (0, a)
   ------------------
        y = -b
```

y la guía continúa infinitamente hacia izquierda y derecha.

El parámetro $\varepsilon$ controla el tamaño del obstáculo. En este código:

$$
\varepsilon = 0.05,\ 0.07,\ 0.09,\ 0.11
$$

Como $R_0=1$:

$$
r = R_0\,\varepsilon = \varepsilon
$$

---

## 2. ¿Qué estamos buscando?

Buscamos una función $u(x,y)$ que resuelva la ecuación de Helmholtz:

$$
-\Delta u = k^2 u
$$

Aquí aparece nuestro primer objeto importante:

$$
\boxed{k = \text{número de onda}}
$$

y $k^2$ es el valor propio espectral.

**No** estamos buscando solamente $u$. Estamos buscando pares $(k^2, u)$ para los cuales existe una solución no trivial.

Las condiciones son:

- Dirichlet en las paredes: $u=0$
- Neumann sobre el obstáculo: $\dfrac{\partial u}{\partial n}=0$
- Decaimiento en el infinito: $u(x,y)\to 0$ cuando $|x|\to\infty$

Esta última condición es crucial: significa que buscamos una onda cuya energía permanezca localizada alrededor del obstáculo, es decir:

$$
\boxed{\text{un modo atrapado}}
$$

---

## 3. Antes del obstáculo: ¿qué modos admite una guía recta?

Aquí aparece la razón matemática de $\Lambda_1$.

Sin obstáculo, separamos variables: $u(x,y)=X(x)Y(y)$. Sustituyendo en $-\Delta u = k^2 u$:

$$
-X''(x)Y(y) - X(x)Y''(y) = k^2 X(x)Y(y)
$$

Dividiendo por $XY$:

$$
-\frac{X''}{X} - \frac{Y''}{Y} = k^2
$$

La parte transversal $Y(y)$ debe satisfacer Dirichlet: $Y(-b)=Y(b)=0$. Ese problema tiene autovalores:

$$
\boxed{\Lambda_n = \left(\frac{n\pi}{2b}\right)^2}, \qquad n=1,2,\ldots
$$

En particular, el primer cutoff de la guía:

$$
\boxed{\Lambda_1 = \left(\frac{\pi}{2b}\right)^2}
$$

Para $b=1$:

$$
\Lambda_1 = \left(\frac{\pi}{2}\right)^2 \approx 2.467401100272\ldots
$$

$$
\sqrt{\Lambda_1} = \frac{\pi}{2} \approx 1.570796326795\ldots
$$

que es exactamente lo que imprime el código.

---

## 4. ¿Por qué $\Lambda_1$ separa propagación y decaimiento?

Para cada modo transversal obtenemos en $x$:

$$
X_n'' + (k^2-\Lambda_n)X_n = 0
$$

**Caso propagante:** si $k^2 > \Lambda_n$, entonces $k^2-\Lambda_n>0$ y

$$
X_n(x) = e^{\pm i\sqrt{k^2-\Lambda_n}\,x}
$$

Eso oscila: es un modo propagante.

**Caso evanescente:** si $k^2 < \Lambda_n$, entonces $k^2-\Lambda_n<0$ y

$$
X_n'' - (\Lambda_n-k^2)X_n = 0
$$

cuya solución tiene la forma

$$
X_n(x) = C_1 e^{-\sqrt{\Lambda_n-k^2}\,|x|} + C_2 e^{+\sqrt{\Lambda_n-k^2}\,|x|}
$$

Como queremos una solución acotada y localizada, conservamos la parte que decae:

$$
X_n(x) \sim e^{-\sqrt{\Lambda_n-k^2}\,|x|}
$$

Por eso, cuando

$$
\boxed{k^2<\Lambda_1}
$$

ni siquiera el primer canal puede propagarse: todos son evanescentes. Eso explica por qué el intervalo del Teorema 2.1 es:

$$
\boxed{0<k^2<\Lambda_1}
$$

Un valor propio en ese intervalo pertenece al espectro discreto.

---

## 5. Aquí aparece $\sigma$

Definimos:

$$
\boxed{\sigma = \sqrt{\Lambda_1-k^2}}
$$

de donde

$$
\boxed{k^2 = \Lambda_1-\sigma^2}
$$

que es exactamente la forma utilizada por el Teorema 2.1.

$\sigma$ tiene además una interpretación física muy bonita. Recuerda que el primer modo decae como $e^{-\sqrt{\Lambda_1-k^2}|x|}$. Como $\sigma=\sqrt{\Lambda_1-k^2}$:

$$
\boxed{u \sim e^{-\sigma|x|}}
$$

lejos del obstáculo, para la contribución dominante. Por tanto:

$$
\boxed{\sigma = \text{tasa de decaimiento longitudinal}}
$$

y aproximadamente $1/\sigma$ es una escala de longitud de localización.

- Si $\sigma$ es grande → el modo decae rápidamente.
- Si $\sigma$ es pequeño → el modo se extiende mucho más por la guía.
- Cuando $\sigma\to 0$ → $k^2\to\Lambda_1$ (el modo se acerca al cutoff).

---

## 6. ¿Qué dice realmente el Teorema 2.1?

Para un obstáculo pequeño, el teorema dice aproximadamente que existe una altura crítica

$$
a^* = a_0^* + O(\varepsilon)
$$

tal que, cuando $a>a^*$, existe **un único modo atrapado discreto**.

A orden dominante:

$$
a_0^* = \frac{2b}{\pi}\arctan\sqrt{\frac{S}{2\pi\mu}}
$$

En el círculo de referencia, $R_0=1$, de modo que $S=\pi$ y $\mu=1$. Entonces:

$$
a_0^* = \frac{2b}{\pi}\arctan\frac{1}{\sqrt{2}}
$$

Con $b=1$ resulta:

$$
\boxed{a_0^* \approx 0.391826552}
$$

y escogimos $a=0.6$. Por tanto, $0.6 > 0.391826552$. Eso explica la salida:

```text
geometric condition a > a0*: PASS
```

Antes de hacer ninguna simulación, **la geometría ya se encuentra en el régimen donde el teorema predice el modo**.

---

## 7. ¿Qué predice el teorema para $\sigma$?

Define $\alpha = \dfrac{\pi a}{b}$. Para el círculo, $\mu=1$, $S=\pi$. La expansión usada en el código es:

$$
\sigma_{\text{asym}} = \varepsilon^2 \frac{\pi^2}{4b^3}\left[ \pi\mu\sin^2\!\left(\frac{\alpha}{2}\right) - \frac{S}{2}\cos^2\!\left(\frac{\alpha}{2}\right) \right]
$$

Más precisamente, el teorema tiene un resto:

$$
\sigma = \varepsilon^2\, C(a,b,\Gamma) + O(\varepsilon^3\ln\varepsilon)
$$

Eso significa que la fórmula que usamos en el código es **la parte principal de una expansión asintótica**, no el valor exacto. Por eso la llamamos $\boxed{\sigma_{\text{asym}}}$.

Una vez que la tenemos:

$$
k_{\text{asym}}^2 = \Lambda_1-\sigma_{\text{asym}}^2 \quad\Rightarrow\quad k_{\text{asym}} = \sqrt{\Lambda_1-\sigma_{\text{asym}}^2}
$$

---

## 8. Un detalle extremadamente importante: el gap es de orden $\varepsilon^4$

El teorema dice $\sigma=O(\varepsilon^2)$. Pero $\Lambda_1-k^2=\sigma^2$. Por tanto:

$$
\boxed{\Lambda_1-k^2 = O(\varepsilon^4)}
$$

y esto explica casi todos los problemas numéricos que tuvimos. Aunque $\varepsilon$ no parezca extremadamente pequeño, el autovalor queda **muchísimo más cerca** de $\Lambda_1$.

Para $\varepsilon=0.05$, por ejemplo:

$$
\sigma_{\text{asym}} \approx 0.009336 \quad\Rightarrow\quad \sigma_{\text{asym}}^2 \approx 8.72\times10^{-5}
$$

Por tanto:

$$
k^2 = \Lambda_1 - 8.72\times10^{-5}
$$

Y en $k$, la diferencia con el cutoff es solamente:

$$
\sqrt{\Lambda_1}-k \approx 2.77\times10^{-5}
$$

Por eso en tus gráficas el modo parece **pegado a $\pi/2$**. No es un error del gráfico. Es exactamente lo que predice la teoría.

---

## 9. Tenemos una predicción analítica. ¿Cómo resolvemos el problema real?

Aquí entra **BEM** (Método de Elementos de Frontera).

El dominio es infinito ($-\infty<x<\infty$). Resolverlo directamente con un método de dominio implicaría truncar la guía en algún punto.

BEM permite evitar eso utilizando una función de Green apropiada. La idea es transformar la PDE $-\Delta u = k^2 u$ en una ecuación integral definida solamente sobre la frontera del obstáculo.

---

## 10. ¿Qué es la función de Green?

Podemos pensar en $G(P,q)$ como la respuesta de la guía en el punto $P$ ante una fuente puntual colocada en $q$. Satisface una ecuación del tipo:

$$
-\Delta G - k^2 G = \delta(P-q)
$$

con las condiciones apropiadas en las paredes de la guía.

El código utiliza:

```python
lattice.greens_dirichlet(...)
```

por lo que se está usando una función de Green que ya incorpora las condiciones de Dirichlet de las paredes. Eso es muy conveniente porque significa que las paredes no necesitan discretizarse explícitamente.

> El archivo que estamos estudiando no contiene la fórmula interna de esta función de Green; esa matemática está dentro de `lattice_sums.py`. Aquí vemos exactamente **cómo se usa**, pero para explicar la serie de Green término por término habría que mirar ese archivo también.

---

## 11. ¿Por qué aparece `y + b + a`?

El círculo se parametriza en coordenadas centradas en el obstáculo:

$$
X(t)=\varepsilon\cos t, \qquad Y(t)=\varepsilon\sin t
$$

Pero físicamente el centro está en $y=a$. Entonces la coordenada global es:

$$
y_{\text{global}} = Y+a
$$

Además, la implementación de Green trabaja con una coordenada desplazada donde la pared inferior corresponde a cero. Como la pared inferior física es $y=-b$:

$$
y_{\text{Green}} = y_{\text{global}} + b = Y + a + b
$$

Por eso aparece:

```python
y + b + config.a
```

y análogamente para la coordenada fuente. No es un término físico nuevo, es solamente un cambio de coordenadas.

---

## 12. De la PDE a la integral de frontera

Usando la identidad de Green se llega a una representación de la forma:

$$
u(P) = \int_\Gamma \left[ u(q)\frac{\partial G(P,q)}{\partial n_q} - G(P,q)\frac{\partial u(q)}{\partial n_q} \right] ds_q
$$

Pero sobre el obstáculo tenemos Neumann ($\partial u/\partial n=0$), entonces desaparece el segundo término:

$$
u(P) = \int_\Gamma u(q)\frac{\partial G(P,q)}{\partial n_q}\, ds_q
$$

Ahora hacemos que $P$ se acerque a la frontera. Aparece la conocida relación de salto del potencial de doble capa:

$$
\boxed{\frac12 u(p) = \int_\Gamma u(q)\frac{\partial G(p,q)}{\partial n_q}\, ds_q}
$$

con la convención de normales usada en tu implementación. Esta es **la ecuación fundamental que resuelve tu BEM**.

---

## 13. Aquí aparece $K$

Tenemos $k$, minúscula:

$$
\boxed{k = \text{número de onda}}
$$

Pero también podemos definir un operador integral:

$$
(Ku)(p) = \int_\Gamma u(q)\frac{\partial G(p,q)}{\partial n_q}\, ds_q
$$

Entonces $\tfrac12 u = Ku$. Este:

$$
\boxed{K = \text{operador integral de frontera}}
$$

no es el número de onda. En el código aparece como `K_weighted`, que es la versión discretizada del núcleo de ese operador.

---

## 14. ¿Qué es el núcleo?

Parametrizamos el círculo: $q(\theta) = (X(\theta), Y(\theta))$, con

$$
X(\theta)=\varepsilon\cos\theta, \qquad Y(\theta)=\varepsilon\sin\theta
$$

Sus derivadas son $X'=-\varepsilon\sin\theta$, $Y'=\varepsilon\cos\theta$. El elemento de arco es $ds = \sqrt{X'^2+Y'^2}\,d\theta$. Definimos $w(\theta)=\sqrt{X'^2+Y'^2}$; para un círculo, $w=\varepsilon$.

La normal multiplicada por $ds$ puede escribirse, dependiendo de la orientación, como $(-Y', X')\,d\theta$. Entonces:

$$
\frac{\partial G}{\partial n_q}\,ds = \nabla_q G\cdot n\, ds = -Y'G_\xi + X'G_\eta
$$

Y eso es exactamente lo que ves:

```python
return xi_p * G_eta - eta_p * G_xi
```

Es decir:

$$
\boxed{K^w(p,q) = X'(\theta)G_\eta - Y'(\theta)G_\xi}
$$

es el núcleo normal **incluyendo el Jacobiano de la parametrización**. Por eso se llama `weighted_normal_kernel`.

---

## 15. ¿De dónde salen $G_\xi$ y $G_\eta$?

Necesitamos derivar Green con respecto al punto fuente: $\partial G/\partial\xi$, $\partial G/\partial\eta$.

El código no tiene fórmulas analíticas para estas derivadas; las aproxima por diferencias centradas:

$$
G_\xi \approx \frac{G(\xi+h,\eta)-G(\xi-h,\eta)}{2h}, \qquad G_\eta \approx \frac{G(\xi,\eta+h)-G(\xi,\eta-h)}{2h}
$$

con $h=10^{-6}$. En el código:

```python
finite_difference_step = 1.0e-6
```

Este es otro parámetro numérico.

---

## 16. ¿Qué pasa cuando $p=q$?

Cuando el punto de observación y el punto fuente coinciden, el núcleo de Green presenta una singularidad. No podemos evaluar ingenuamente $\partial G/\partial n$ ahí.

Por eso el código tiene dos funciones, `green` y `green_regularized`, y cuando $p=q$ usa `G_regularized` más un término geométrico:

$$
\frac{X''Y'-Y''X'}{4\pi w^2}
$$

Su función es incorporar correctamente el límite de la parte singular del núcleo. Por eso en la matriz BEM los elementos diagonales requieren un tratamiento distinto de los elementos fuera de la diagonal.

---

## 17. Ahora discretizamos la frontera

Aquí aparece $\boxed{M}$: el número de puntos utilizados para discretizar el contorno del obstáculo.

El código toma:

$$
\theta_j = \left(j+\frac12\right)\frac{2\pi}{M}
$$

Son puntos medios uniformemente distribuidos alrededor del círculo. Con ellos aproximamos:

$$
\int_0^{2\pi} f(\theta)\,d\theta \approx \frac{2\pi}{M}\sum_{j=1}^{M} f(\theta_j)
$$

Entonces:

$$
\frac12 u_i = \frac{2\pi}{M}\sum_j K^w_{ij}u_j
$$

Multiplicamos por 2:

$$
u_i = \frac{4\pi}{M}\sum_j K^w_{ij}u_j
$$

Pasamos todo al lado izquierdo:

$$
u_i - \frac{4\pi}{M}\sum_j K^w_{ij}u_j = 0
$$

En forma matricial:

$$
\boxed{A(k)\,\mathbf{u}=0}
$$

donde

$$
\boxed{A(k) = I - \frac{4\pi}{M}K^w(k)}
$$

Eso es exactamente:

```python
return np.eye(M) - (4.0 * PI / M) * K_weighted
```

---

## 18. ¿Por qué escribimos $A(k)$ y no solamente $A$?

Porque Green depende de $k$. Por tanto $G=G(k)$, luego $K=K(k)$, y finalmente $A=A(k)$.

Cada vez que cambias $k$, **construyes un problema BEM diferente**. Esta es la razón por la que el código es costoso. Para cada valor candidato de $k$:

$$
k \to G(k) \to K(k) \to A(k) \to \text{SVD}
$$

---

## 19. ¿Cómo aparece el valor propio?

Tenemos $A(k)\mathbf{u}=0$. Para un $k$ cualquiera, normalmente $A(k)$ es invertible, y entonces $\mathbf{u}=0$ es la única solución.

Pero queremos $\mathbf{u}\neq 0$. Eso solo puede ocurrir cuando $A(k)$ pierde rango, es decir, cuando es singular:

$$
\boxed{\det A(k)=0}
$$

Esa es la condición original que escribiste en el paper.

---

## 20. ¿Por qué el código ya no busca directamente $\det A(k)=0$?

Porque numéricamente el determinante puede ser una medida bastante mala: una matriz puede estar muy cerca de ser singular mientras su determinante sea difícil de interpretar debido a escalas numéricas.

Por eso usamos SVD. Descomponemos:

$$
A = U\Sigma V^*, \qquad \Sigma = \operatorname{diag}(s_1,s_2,\ldots,s_M), \qquad s_1\geq s_2\geq\cdots\geq s_M\geq 0
$$

(Llamamos $s_j$ para evitar confusión con el $\sigma$ del teorema.)

Si $A$ es singular, $s_M=0$. Entonces definimos:

$$
\boxed{\sigma_{\min}(A) = s_M}
$$

y buscamos $\sigma_{\min}(A(k)) \approx 0$.

---

## 21. Hay DOS sigmas completamente diferentes

Esto es importantísimo.

**El primero:**

$$
\boxed{\sigma = \sqrt{\Lambda_1-k^2}}
$$

pertenece al Teorema 2.1. Es una cantidad espectral/física: mide la distancia al cutoff y la tasa de decaimiento.

**El segundo:**

$$
\boxed{\sigma_{\min}(A)}
$$

pertenece al álgebra lineal numérica. Es el menor valor singular de la matriz BEM. Sirve para detectar si $A(k)$ está cerca de ser singular.

No tienen el mismo significado. Podrías incluso cambiar el nombre en el código para evitar confusión — usar `sv_min` en lugar de `sigma_min` — y probablemente sería lo mejor. Entonces tendríamos claramente:

$$
\sigma_{\text{BEM}} = \sqrt{\Lambda_1-k_{\text{BEM}}^2} \qquad\text{frente a}\qquad sv_{\min}(A(k))
$$

Mucho más claro.

---

## 22. ¿Qué es `kb`?

El código no trabaja directamente con $k$, sino frecuentemente con:

$$
\boxed{kb = k\,b}
$$

Es simplemente el número de onda adimensionalizado por el semi-ancho de la guía. Como $[k]=L^{-1}$ y $[b]=L$, entonces $kb$ no tiene dimensiones. Por eso es muy conveniente.

El primer cutoff en esta variable es:

$$
\sqrt{\Lambda_1}\,b = \frac{\pi}{2b}\,b = \boxed{\frac{\pi}{2}}
$$

independientemente del valor de $b$.

En tu experimento, como $b=1$, numéricamente $kb=k$. Pero conceptualmente no son la misma variable: si pusieras $b=2$, entonces $kb=2k$.

---

## 23. ¿Qué es `delta`?

Para hacer la búsqueda numérica introdujimos:

$$
\boxed{\delta = \Lambda_1-k^2}
$$

Pero recuerda que $\sigma^2=\Lambda_1-k^2$. Entonces:

$$
\boxed{\delta = \sigma^2}
$$

Esto fue clave para mejorar el código. En lugar de buscar ingenuamente en $k$, podemos pensar en la distancia espectral al cutoff. El teorema nos da $\sigma_{\text{asym}}$, así que sabemos aproximadamente que:

$$
\delta_{\text{asym}} = \sigma_{\text{asym}}^2
$$

---

## 24. ¿Cómo buscamos el modo esperado?

El código toma la predicción $\delta_{\text{asym}}=\sigma_{\text{asym}}^2$ y construye una ventana alrededor:

$$
0.2\,\delta_{\text{asym}} \leq \delta \leq 5\,\delta_{\text{asym}}
$$

Es decir, no asumimos que la fórmula asintótica sea exacta: permitimos que el verdadero modo esté bastante a un lado u otro.

Después convertimos otra vez $\delta\to k$ mediante:

$$
k = \sqrt{\Lambda_1-\delta} \qquad\Rightarrow\qquad kb = b\sqrt{\Lambda_1-\delta}
$$

---

## 25. ¿Qué minimizamos exactamente?

Dentro de ese intervalo buscamos:

$$
\min_k\ \sigma_{\min}(A(k))
$$

Pero el código realmente minimiza:

$$
\log_{10}\left[\sigma_{\min}(A(k))\right]
$$

Eso no cambia dónde está el mínimo porque el logaritmo es monótono ($x_1<x_2 \Rightarrow \log x_1 < \log x_2$), pero numéricamente ayuda cuando los valores varían muchos órdenes de magnitud.

Por eso usamos `minimize_scalar(...)` con el método `"bounded"`. El resultado es $\boxed{k_{\text{BEM}}}$.

---

## 26. A partir de $k_{\text{BEM}}$, obtenemos $\sigma_{\text{BEM}}$

Una vez encontrado numéricamente $k_{\text{BEM}}$, calculamos:

$$
\boxed{\sigma_{\text{BEM}} = \sqrt{\Lambda_1-k_{\text{BEM}}^2}}
$$

Esta es la cantidad numérica que tiene sentido comparar directamente con el teorema:

$$
\sigma_{\text{asym}} \quad\text{vs}\quad \sigma_{\text{BEM}}
$$

Esto es mejor que comparar solamente $k_{\text{asym}}$ vs $k_{\text{BEM}}$, porque ambos $k$ están extremadamente cerca de $\sqrt{\Lambda_1}$. Por ejemplo: $1.5707686$ frente a $1.5707692$. El error relativo de $k$ parece diminuto, pero al transformar a $\sigma$ vemos realmente la diferencia en la distancia al cutoff.

---

## 27. ¿Por qué variamos $M$?

Porque una solución BEM depende de cuán finamente discretizamos el obstáculo. Si usamos $M=16$, tenemos una matriz $16\times16$; si usamos $M=48$, tenemos $48\times48$.

Queremos comprobar que el resultado no sea un artefacto de una discretización concreta. Por eso hacemos $M = 16,\ 24,\ 32,\ 40,\ 48$.

Por ejemplo, para $\varepsilon=0.09$:

| $M$ | $kb$           |
| --- | -------------- |
| 16  | 1.570529956631 |
| 24  | 1.570529958280 |
| 32  | 1.570529959331 |
| 40  | 1.570529959704 |
| 48  | 1.570529959870 |

Eso es una evidencia muy fuerte de convergencia: $k_M \to k_{\text{BEM}}$.

---

## 28. ¿Qué significa `mesh_change`?

El código compara $\sigma_{\text{BEM}}^{(M)}$ entre dos discretizaciones consecutivas. Calcula aproximadamente:

$$
E_M = \frac{|\sigma_M-\sigma_{M-1}|}{|\sigma_M|}
$$

Si esto es pequeño, $\sigma_M$ ya prácticamente no cambia al aumentar $M$. El código exige $E_M\leq 3\%$.

Tus resultados están muchísimo mejor que eso:

```text
mesh_change=0.000%
```

redondeado a tres decimales.

---

## 29. ¿Qué es `drop_factor`?

Esto fue introducido para no aceptar cualquier mínimo numérico insignificante.

Supongamos que buscamos en $[k_L,k_R]$ y encontramos un mínimo en $k_*$. Calculamos $s_L=\sigma_{\min}(A(k_L))$, $s_R=\sigma_{\min}(A(k_R))$ y $s_*=\sigma_{\min}(A(k_*))$. Luego:

$$
\boxed{D = \frac{\min(s_L,s_R)}{s_*}}
$$

Ese es `drop_factor`.

- Si $D=1$, prácticamente no existe valle.
- Si $D=10$, el mínimo es diez veces más bajo.
- Si $D=100$, dos órdenes de magnitud.

Tus resultados tienen, por ejemplo, $D=6.05\times10^4$. Eso significa que la medida de singularidad cae alrededor de $60{,}500$ veces respecto a los bordes de la región. Es un valle espectral extremadamente marcado.

---

## 30. ¿Por qué no exigimos ya $\sigma_{\min}(A)<10^{-6}$?

Porque descubrimos que eso era demasiado ingenuo. Para $\varepsilon=0.05$, ahora tenemos:

$$
\sigma_{\min}(A) = 1.881\times10^{-5}
$$

Según el criterio viejo, $1.881\times10^{-5} > 10^{-6}$ y habríamos dicho erróneamente: _"no hay modo"_.

Pero ahora observamos simultáneamente: $D\approx 3290$, el mínimo está en el interior, $k_M$ converge perfectamente con $M$, y no aparecen otros candidatos. Es evidencia numérica mucho más convincente. Por eso ahora un modo se considera resuelto cuando se combinan varios diagnósticos.

---

## 31. ¿Qué condiciones usa el código para decir "expected mode resolved"?

La lógica es:

$$
\text{resolved} = \text{interior} \ \land\ \text{near singular} \ \land\ \text{deep minimum} \ \land\ \text{mesh converged}
$$

Concretamente:

- $\sigma_{\min}(A) \leq 10^{-4}$
- $D \geq 100$
- $E_M \leq 3\%$
- el mínimo no puede estar simplemente en uno de los extremos de la ventana

Es importante entender que esos números ($10^{-4}$, $100$, $3\%$) son **criterios numéricos de validación**, no constantes matemáticas del Teorema 2.1.

---

## 32. Verificamos existencia. ¿Pero cómo verificamos que sea único?

Buscar solamente alrededor de $k_{\text{asym}}$ puede demostrar que encontramos **el modo esperado**. Pero el teorema dice que hay **uno solo**. Entonces preguntamos: ¿hay otro mínimo de $A(k)$ en alguna otra parte de $0<k^2<\Lambda_1$?

Por eso hacemos el _whole-band uniqueness screen_.

---

## 33. ¿Por qué no usamos solamente una malla uniforme en $k$?

Porque el modo está pegado al cutoff. Recuerda: $\Lambda_1-k^2=O(\varepsilon^4)$.

Una malla como $0,\ 0.05,\ 0.10,\ldots,1.55$ podría saltarse completamente un mínimo situado en $1.57069$.

Por eso el código combina dos muestreos: uno aproximadamente uniforme en $kb$, para cubrir la banda completa, y otro logarítmico en $\delta=\Lambda_1-k^2$. El logarítmico puede producir puntos con distancias espectrales al cutoff como $10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}, 10^{-5},\ldots$ Eso da mucha resolución cerca de $\Lambda_1$.

---

## 34. ¿Qué significa la gráfica que me mostraste?

En el eje horizontal tenemos $kb$; en el vertical, $\sigma_{\min}(A(k))$, con escala vertical logarítmica.

Durante casi toda la banda, $\sigma_{\min}(A(k)) \sim O(1)$: eso significa que $A(k)$ está bastante lejos de ser singular.

Pero cerca de $kb\approx 1.5707$ la curva cae violentamente. Ese valle significa que $A(k)$ se aproxima a una matriz singular; por eso ahí encontramos el valor propio.

- La línea vertical discontinua es $kb_{\text{asym}}$.
- La línea vertical punteada es $kb_{\text{BEM}}$.

Cuando están casi superpuestas, significa que la predicción asintótica y el cálculo numérico están muy cerca.

---

## 35. ¿Cómo buscamos modos adicionales?

Después de calcular $\sigma_{\min}(A(k))$ sobre toda la banda, buscamos puntos que sean menores que sus vecinos:

$$
s_i \leq s_{i-1}, \qquad s_i \leq s_{i+1}
$$

Esos son candidatos a mínimos locales. Sabemos que uno corresponde al modo esperado, así que ignoramos el intervalo que se superpone con su ventana. Si hubiera otro valle en otra parte de la banda, se refina de nuevo.

---

## 36. ¿Por qué usamos $M=16$ y después $M=32$ para esos candidatos?

Porque no queremos hacer BEM con $M=48$ en cientos de puntos sin necesidad. Primero hacemos un _screening_ barato con $M=16$. Si aparece algo sospechoso, repetimos con $M=32$ y preguntamos si su posición permanece aproximadamente igual:

$$
|kb_{32}-kb_{16}| \leq 0.002
$$

También tiene que ser un mínimo profundo y casi singular. Solo entonces lo consideramos un posible modo adicional persistente.

En tus cuatro casos obtuvimos:

```text
sampled local minima: 1 total
0 outside expected-mode window
```

Eso es exactamente lo que queríamos ver.

---

## 37. ¿Qué significa "exactly one resolved discrete mode: PASS"?

Significa numéricamente: primero encontramos el modo predicho por el teorema; luego exploramos $0<k^2<\Lambda_1$ y no encontramos otros candidatos persistentes. Por tanto:

$$
\boxed{N_{\text{resolved}}=1}
$$

No es una nueva demostración matemática de unicidad — esa viene del teorema — pero es una validación numérica consistente con esa afirmación.

---

## 38. Finalmente comparamos teoría y BEM

Tenemos $\sigma_{\text{asym}}$ de la fórmula teórica y $\sigma_{\text{BEM}}$ de $k_{\text{BEM}}$. Calculamos:

$$
\boxed{E_\sigma = \frac{|\sigma_{\text{BEM}}-\sigma_{\text{asym}}|}{|\sigma_{\text{asym}}|}}
$$

y decidimos usar como criterio de precisión $E_\sigma\leq 5\%$. De nuevo: $\boxed{5\%}$ no viene del teorema, es una tolerancia que nosotros elegimos para definir que la aproximación sea suficientemente precisa.

---

## 39. Tus resultados tienen una historia matemática muy bonita

| $\varepsilon$ | Único modo | $E_\sigma$ |
| ------------- | ---------- | ---------- |
| 0.05          | ✔️         | 1.096%     |
| 0.07          | ✔️         | 2.433%     |
| 0.09          | ✔️         | 4.371%     |
| 0.11          | ✔️         | 6.875%     |

Lo que esto dice es:

$$
\boxed{\text{el modo existe en los cuatro experimentos}}
$$

pero

$$
\boxed{\text{la fórmula asintótica pierde precisión al aumentar } \varepsilon}
$$

Y eso es exactamente lo que esperaríamos: la teoría fue derivada suponiendo $\varepsilon\ll1$, así que cuando $\varepsilon\to0$ esperamos $\sigma_{\text{asym}}\to\sigma_{\text{verdadero}}$.

Efectivamente observas:

$$
1.096\% < 2.433\% < 4.371\% < 6.875\%
$$

La discrepancia crece al alejarnos del régimen asintótico.

---

## 40. Toda la historia del código en una sola ecuación conceptual

El teorema dice:

$$
\boxed{\sigma_{\text{asym}}(\varepsilon,a) \longrightarrow k_{\text{asym}} = \sqrt{\Lambda_1-\sigma_{\text{asym}}^2}}
$$

Mientras que BEM hace:

$$
\boxed{k \longrightarrow G(k) \longrightarrow K(k) \longrightarrow A(k) \longrightarrow \sigma_{\min}(A(k))}
$$

y busca

$$
\boxed{k_{\text{BEM}} = \arg\min_k\ \sigma_{\min}(A(k))}
$$

cerca de una singularidad. Luego transforma:

$$
\boxed{\sigma_{\text{BEM}} = \sqrt{\Lambda_1-k_{\text{BEM}}^2}}
$$

y finalmente compara:

$$
\boxed{\sigma_{\text{asym}} \quad\text{vs}\quad \sigma_{\text{BEM}}}
$$

Mientras tanto, un barrido independiente verifica que no aparezcan otros mínimos resueltos en $0<k^2<\Lambda_1$.

---

## El mapa mental que quiero que conserves

Cuando mires este código, piensa en esta cadena:

$$
\boxed{\varepsilon,\, a,\, b} \ \text{determinan la geometría.}
$$

De ahí el teorema produce $\boxed{\sigma_{\text{asym}}}$ y después $\boxed{k_{\text{asym}}}$.

Independientemente, para cada $k$ de prueba BEM construye $\boxed{A(k)}$. El verdadero autovalor hace que $\boxed{A(k)\ \text{sea singular}}$. Lo detectamos mediante $\boxed{\sigma_{\min}(A(k))\approx 0}$. Eso produce $\boxed{k_{\text{BEM}}}$.

Luego calculamos:

$$
\boxed{\sigma_{\text{BEM}} = \sqrt{\Lambda_1-k_{\text{BEM}}^2}}
$$

Y finalmente:

$$
\boxed{\text{Teoría} \quad \sigma_{\text{asym}} \ \overset{?}{\approx}\ \sigma_{\text{BEM}} \quad \text{Numérica}}
$$

> **Idea clave para tu presentación:** BEM no usa la fórmula asintótica para producir el valor propio. La fórmula asintótica solamente nos dice aproximadamente dónde debería estar. El valor $k_{\text{BEM}}$ se obtiene independientemente, buscando dónde el sistema integral discretizado pierde invertibilidad. Esa independencia es precisamente lo que hace que la comparación sea una validación numérica real del resultado teórico.

---

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python first_theorem_2.py
=== Focused validation of Theorem 2.1 ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1) b = 1.570796326795
epsilons = (0.05, 0.07, 0.09, 0.11)
refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
predicted kb = 1.570768582119
predicted sigma = 9.33604312e-03
predicted kb cutoff gap = 2.774e-05
expected-mode mesh refinement:
M=16: kb=1.570769187049, sigma_BEM=9.23370398e-03, sigma_min=1.898e-05, drop=3.26e+03, mesh_change=--, interior=yes
M=24: kb=1.570769187077, sigma_BEM=9.23369931e-03, sigma_min=1.886e-05, drop=3.28e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570769187082, sigma_BEM=9.23369832e-03, sigma_min=1.883e-05, drop=3.29e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570769187084, sigma_BEM=9.23369802e-03, sigma_min=1.882e-05, drop=3.29e+03, mesh_change=0.000%, interior=yes
M=48: kb=1.570769187085, sigma_BEM=9.23369785e-03, sigma_min=1.881e-05, drop=3.29e+03, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570768582119
kb BEM = 1.570769187085
sigma asym = 9.33604312e-03
sigma BEM = 9.23369785e-03
final sigma_min(A) = 1.881e-05
final minimum drop = 3.29e+03
final mesh change in sigma = 0.000%
relative sigma error = 1.096%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.070 ===
predicted kb = 1.570689740172
predicted sigma = 1.82986445e-02
predicted kb cutoff gap = 1.066e-04
expected-mode mesh refinement:
M=16: kb=1.570694863808, sigma_BEM=1.78534333e-02, sigma_min=1.415e-06, drop=6.05e+04, mesh_change=--, interior=yes
M=24: kb=1.570694864415, sigma_BEM=1.78533799e-02, sigma_min=1.415e-06, drop=6.05e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570694864560, sigma_BEM=1.78533671e-02, sigma_min=1.415e-06, drop=6.05e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570694864610, sigma_BEM=1.78533627e-02, sigma_min=1.415e-06, drop=6.05e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570694864633, sigma_BEM=1.78533607e-02, sigma_min=1.415e-06, drop=6.05e+04, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570689740172
kb BEM = 1.570694864633
sigma asym = 1.82986445e-02
sigma BEM = 1.78533607e-02
final sigma_min(A) = 1.415e-06
final minimum drop = 6.05e+04
final mesh change in sigma = 0.000%
relative sigma error = 2.433%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.090 ===
predicted kb = 1.570505049848
predicted sigma = 3.02487797e-02
predicted kb cutoff gap = 2.913e-04
expected-mode mesh refinement:
M=16: kb=1.570529956631, sigma_BEM=2.89267281e-02, sigma_min=4.073e-06, drop=2.65e+04, mesh_change=--, interior=yes
M=24: kb=1.570529958280, sigma_BEM=2.89266386e-02, sigma_min=4.428e-06, drop=2.44e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570529959331, sigma_BEM=2.89265815e-02, sigma_min=4.267e-06, drop=2.53e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570529959704, sigma_BEM=2.89265612e-02, sigma_min=4.210e-06, drop=2.56e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570529959870, sigma_BEM=2.89265522e-02, sigma_min=4.184e-06, drop=2.58e+04, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570505049848
kb BEM = 1.570529959870
sigma asym = 3.02487797e-02
sigma BEM = 2.89265522e-02
final sigma_min(A) = 4.184e-06
final minimum drop = 2.58e+04
final mesh change in sigma = 0.000%
relative sigma error = 4.371%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.110 ===
predicted kb = 1.570146262336
predicted sigma = 4.51864487e-02
predicted kb cutoff gap = 6.501e-04
expected-mode mesh refinement:
M=16: kb=1.570232583665, sigma_BEM=4.20800840e-02, sigma_min=3.197e-06, drop=4.00e+04, mesh_change=--, interior=yes
M=24: kb=1.570232591832, sigma_BEM=4.20797792e-02, sigma_min=3.176e-06, drop=4.03e+04, mesh_change=0.001%, interior=yes
M=32: kb=1.570232593785, sigma_BEM=4.20797064e-02, sigma_min=3.171e-06, drop=4.04e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570232594476, sigma_BEM=4.20796806e-02, sigma_min=3.170e-06, drop=4.04e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570232594781, sigma_BEM=4.20796692e-02, sigma_min=3.169e-06, drop=4.04e+04, mesh_change=0.000%, interior=yes
whole-band uniqueness screen: 72 points (M=16)
16/72
32/72
48/72
64/72
72/72
sampled local minima: 1 total; 0 outside expected-mode window

--- validation result ---
expected mode resolved: PASS
additional resolved modes: 0
uniqueness screen: PASS
exactly one resolved discrete mode: PASS
kb asymptotic = 1.570146262336
kb BEM = 1.570232594781
sigma asym = 4.51864487e-02
sigma BEM = 4.20796692e-02
final sigma_min(A) = 3.169e-06
final minimum drop = 4.04e+04
final mesh change in sigma = 0.000%
relative sigma error = 6.875%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE/UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: 0 < k^2 < Lambda_1
Existence + no additional resolved modes: 4/4
Leading asymptotic sigma within 5.0%: 3/4
epsilon=0.050: unique=PASS, asymptotic=PASS, sigma_min=1.881e-05, mesh_change=0.000%, sigma_error=1.096%
epsilon=0.070: unique=PASS, asymptotic=PASS, sigma_min=1.415e-06, mesh_change=0.000%, sigma_error=2.433%
epsilon=0.090: unique=PASS, asymptotic=PASS, sigma_min=4.184e-06, mesh_change=0.000%, sigma_error=4.371%
epsilon=0.110: unique=PASS, asymptotic=FAIL, sigma_min=3.169e-06, mesh_change=0.000%, sigma_error=6.875%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_1_refined_validation
