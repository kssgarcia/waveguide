Sí. Estos resultados del Teorema 2.3(iii) se ven **muy consistentes numéricamente**, y además cuentan una historia incluso más clara que en el Teorema 2.1.

La conclusión principal es:

[
\boxed{\text{Para los cuatro }\varepsilon\text{ se detecta un unico BIC no radiante.}}
]

pero la aproximación asintótica

[
\sigma\_{\mathrm{asym}}
======================

\varepsilon^2\frac{\pi^3}{b^3}\mu
]

solo permanece dentro de nuestra tolerancia arbitraria del (5%) para los valores más pequeños de (\varepsilon).

## 1. Primero: las hipótesis geométricas se cumplen perfectamente

El código obtiene

[
a=0,
]

[
\mu=1.002199318791,
]

[
\nu=-1.69\times10^{-17}\approx0,
]

y además

[
R_{X,\mathrm{even}}=0,
\qquad
R_{Y,\mathrm{odd}}=0.
]

Es exactamente la situación del punto (iii):

[
X(-t)=X(t),
\qquad
Y(-t)=-Y(t),
\qquad
a=0,
\qquad
\nu=0.
]

Por eso el resultado

```text
Theorem symmetry conditions: PASS
```

está totalmente justificado.

---

## 2. El BIC de (\varepsilon=0.05) está muy bien resuelto

Para

[
\varepsilon=0.05
]

la teoría predice

[
kb\_{\mathrm{asym}}
==================

3.140631984100,
]

mientras BEM encuentra

[
kb\_{\mathrm{BEM}}
=================

3.140675361277.
]

Están extremadamente próximos a

[
\sqrt{\Lambda_2}b=\pi,
]

tal como predice

[
k^2=\Lambda_2-\sigma^2.
]

Más importante todavía, el resultado converge prácticamente de inmediato cuando aumentas (M):

[
\begin{aligned}
M=16 &: 3.140675347570,\
M=24 &: 3.140675357638,\
M=32 &: 3.140675360048,\
M=40 &: 3.140675360901,\
M=48 &: 3.140675361277.
\end{aligned}
]

Eso es **excelente estabilidad de malla**.

---

## 3. El mínimo de la matriz es realmente profundo

Para (\varepsilon=0.05):

[
s\_{\min}(A)
===========

5.611\times10^{-7},
]

y el `drop factor` es

[
1.02\times10^6.
]

Eso significa que el menor valor singular cae aproximadamente **un millón de veces** respecto a los alrededores.

Entonces no estamos viendo una pequeña fluctuación numérica. Tenemos un valle espectral extremadamente pronunciado.

Y pasa lo mismo para los demás:

[
\begin{array}{c|c|c}
\varepsilon&s_{\min}(A)&\text{drop}\
\hline
0.05&5.61\times10^{-7}&1.02\times10^6\
0.07&1.81\times10^{-6}&3.28\times10^5\
0.09&6.66\times10^{-7}&9.46\times10^5\
0.11&6.30\times10^{-7}&1.08\times10^6
\end{array}
]

Es muy buena evidencia de singularidades de (A(k)).

---

# 4. Y aqui viene lo más importante: no parece una onda propagante ordinaria

Esta es la parte que realmente distingue el Teorema 2.3 del Teorema 2.1.

Estamos en

[
\Lambda_1<k^2<\Lambda_2.
]

Por tanto, el primer canal de la guía **sí está abierto**.

En principio la energía podría escapar.

Pero obtuviste:

[
R\_{\mathrm{prop}}
=================

5.806\times10^{-7}
]

para (\varepsilon=0.05),

y luego

[
1.273\times10^{-6},
\quad
2.392\times10^{-6},
\quad
4.083\times10^{-6}.
]

Son valores extremadamente pequeños.

En otras palabras,

[
\boxed{
\text{el coeficiente del primer canal propagante es practicamente cero}.
}
]

Esa es justamente la propiedad física fundamental del BIC: **su frecuencia pertenece al continuo, pero no acopla de manera apreciable al canal que podría transportar energía al infinito**.

Para (\varepsilon=0.05), el log reporta simultáneamente singularidad, cancelación propagante y unicidad.

---

# 5. La paridad también sale espectacularmente bien

Obtienes aproximadamente

[
R_{\mathrm{odd}}
\sim10^{-10}.
]

Por ejemplo:

[
1.255\times10^{-10},
\quad
1.306\times10^{-10},
\quad
1.342\times10^{-10},
\quad
1.597\times10^{-10}.
]

Es prácticamente precisión numérica.

Esto significa que la solución que está encontrando BEM tiene el carácter de simetría esperado.

Esto también explica intuitivamente por qué el canal propagante se cancela.

El primer modo transversal es par:

[
\phi_1(-y)=\phi_1(y),
]

mientras que nuestro estado pertenece al sector impar:

[
u(x,-y)\approx-u(x,y).
]

Entonces su proyección es

[
\int_{-b}^b
u(x,y)\phi_1(y),dy.
]

El producto

[
u\phi_1
]

es impar, así que idealmente

[
\boxed{
\int_{-b}^b
u(x,y)\phi_1(y),dy=0.
}
]

Esto conecta de forma preciosa:

[
\boxed{
\text{simetria}
\rightarrow
\text{ortogonalidad}
\rightarrow
\text{cancelacion del canal propagante}
\rightarrow
\text{BIC}.
}
]

Y tus números están mostrando exactamente eso.

---

# 6. La unicidad también parece correcta

Para cada (\varepsilon), el código encuentra el BIC esperado y después explora el resto de

[
\Lambda_1<k^2<\Lambda_2.
]

Los cuatro "mínimos adicionales" que aparecen cerca de

[
kb=\pi
]

**no son modos**.

Mira sus propiedades:

[
s_{\min}\approx0.96-0.98,
]

[
\text{drop}\approx1,
]

[
R_{\mathrm{odd}}\approx2,
]

[
R_{\mathrm{prop}}\approx1.
]

Es prácticamente lo opuesto al BIC.

Por ejemplo, para (\varepsilon=0.05), esos candidatos tienen (s\_{\min}\approx0.978), drop (=1), residual de paridad (=2) y componente propagante (=1).

Por tanto, están siendo correctamente descartados.

De hecho, yo los llamaría más bien **artefactos del muestreo cerca del cutoff** que "local minima" en la descripción del paper.

El buscador detecta pequeñas fluctuaciones porque tienes muchos puntos extremadamente agrupados cerca de

[
kb=\pi.
]

Pero al refinarlos queda inmediatamente claro que

[
s_{\min}\sim1,
]

por lo que no son candidatos espectrales reales.

---

# 7. La gráfica que muestras es coherente

La segunda figura muestra precisamente lo esperado.

Durante prácticamente toda la banda,

[
s_{\min}(A(k))\sim O(1).
]

Pero cuando llegas al valor

[
kb_{\mathrm{BEM}}\approx3.140675,
]

aparece una caída abrupta de varios órdenes de magnitud.

Luego está

[
\sqrt{\Lambda_2}b
=================

\pi
\approx3.141593.
]

Entonces visualmente tienes:

[
kb_{\mathrm{BEM}}
<
\pi,
]

pero extremadamente cerca del cutoff superior.

Eso concuerda exactamente con

[
k^2
===

\Lambda_2-\sigma^2.
]

---

# 8. Ahora viene el resultado científico más interesante: el error asintótico

Tienes:

[
\begin{array}{c|c}
\varepsilon&E_\sigma\
\hline
0.05&2.283%\
0.07&5.188%\
0.09&9.229%\
0.11&14.216%
\end{array}
]

Esto es una secuencia muy limpia.

A medida que

[
\varepsilon
]

crece, el término líder

[
\sigma\_{\mathrm{asym}}
======================

\varepsilon^2
\frac{\pi^3}{b^3}\mu
]

se aleja progresivamente del BEM.

Eso es justamente lo que esperaríamos porque la fórmula completa realmente es

[
\sigma
======

\varepsilon^2
\frac{\pi^3}{b^3}\mu

- O(\varepsilon^3\ln\varepsilon).
  ]

Nosotros estamos ignorando

[
O(\varepsilon^3\ln\varepsilon).
]

Cuando (\varepsilon) aumenta, ese resto deja de ser despreciable.

---

# 9. Muy importante: el BIC NO desaparece en 0.07

Esto debe quedar muy claro en el paper.

Para

[
\varepsilon=0.07
]

obtienes:

[
\boxed{\text{BIC: PASS}}
]

pero

[
E_\sigma=5.188%>5%.
]

Entonces:

[
\boxed{
\text{el BIC sigue existiendo}
}
]

pero

[
\boxed{
\text{la aproximacion lider ya supera nuestra tolerancia del 5%}.
}
]

Exactamente igual para

[
0.09
]

y

[
0.11.
]

El log refleja correctamente esta separación.

---

# 10. Entonces, ¿cuál es (\varepsilon\_{\max})?

Con tus **puntos muestreados**, solo podemos afirmar:

[
\varepsilon=0.05
]

satisface el (5%),

mientras

[
\varepsilon=0.07
]

ya no.

Por tanto:

[
\boxed{
0.05<\varepsilon_{\max}<0.07
}
]

para esta geometría y este criterio de error.

No diría todavía

[
\varepsilon_{\max}=0.05.
]

Eso solamente sería:

> el mayor valor **muestreado** que pasa.

La primera gráfica muestra muy bien que el cruce está justo antes de (0.07).

Una interpolación lineal meramente orientativa con

[
E(0.05)=0.02283
]

y

[
E(0.07)=0.05188
]

coloca el (5%) aproximadamente en

[
\boxed{
\varepsilon\approx0.0687.
}
]

Pero no usaría ese número como resultado definitivo.

Yo haría ahora únicamente:

[
\boxed{
\varepsilon=
0.065,;
0.067,;
0.068,;
0.069,;
0.070
}
]

o incluso una pequeña búsqueda binaria.

No tienes que repetir valores (0.09) y (0.11): ya sabemos que están muy por encima de la tolerancia.

---

# 11. Hay una pequeña mejora que sí haría en el código

Los cuatro supuestos candidatos adicionales cerca de

[
\pi
]

solo generan ruido en el output.

Antes de refinarlos, pondría un filtro coarse. Por ejemplo, **ni siquiera refinar un candidato si**

[
s_{\min}^{\mathrm{coarse}}>10^{-2}
]

o si su `coarse drop` es menor que, digamos,

[
10.
]

Tus candidatos falsos tienen

[
s_{\min}\approx0.96
]

y

[
D\approx1.
]

Así que serían descartados inmediatamente.

Eso haría que en vez de:

```text
sampled local minima: 5 total
4 outside expected-BIC window
checking additional minimum 1/4
...
```

tuvieras algo conceptualmente más limpio:

```text
sampled numerical extrema: 5
spectrally significant minima: 1
expected BIC: 1
additional resolved BIC candidates: 0
```

Eso es más fiel matemáticamente a lo que realmente está pasando.

---

## En resumen

Tus resultados actuales soportan bastante bien las cuatro piezas esenciales del Teorema 2.3(iii):

[
\boxed{
\begin{array}{ll}
1.&a=0,\ X\text{ par},\ Y\text{ impar},\ \nu\approx0,[1mm]
2.&\exists,k^2\in(\Lambda_1,\Lambda_2),[1mm]
3.&\text{el estado es no radiante y tiene la paridad esperada},[1mm]
4.&\text{no se detectan otros BICs resueltos en la banda}.
\end{array}}
]

Y además observas experimentalmente:

[
\boxed{
E_\sigma\uparrow
\quad\text{cuando}\quad
\varepsilon\uparrow,
}
]

con el umbral del (5%) situado aproximadamente alrededor de

[
\boxed{\varepsilon\simeq0.069}
]

para esta geometría, pendiente de refinar numéricamente.

## **Este resultado es bastante más fuerte que simplemente comparar (k*{\rm BEM}) con (k*{\rm asym})**: tienes convergencia de malla, singularidad de la matriz, paridad correcta, cancelación del canal propagante y screening de unicidad, todos apuntando al mismo BIC.

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python second_theorem_1.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iii): x-symmetric BIC ===
b = 1.0
a = 0.0
beta = 0.1
mu = 1.002199318791e+00
nu = -1.687599000907e-17
X-even residual = 0.000e+00
Y-odd residual = 0.000e+00
Theorem symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590
epsilons = (0.05, 0.07, 0.09, 0.11)
refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
predicted kb = 3.140631984100
predicted sigma = 7.76861734e-02
predicted kb gap to sqrt(Lambda_2) = 9.607e-04
expected-BIC mesh refinement:
M=16: kb=3.140675347570, sigma_BEM=7.59128596e-02, sv_min=4.455e-07, drop=1.28e+06, mesh_change=--, odd_res=1.035e-10, Rprop=5.806e-07, interior=yes
M=24: kb=3.140675357638, sigma_BEM=7.59124431e-02, sv_min=5.307e-07, drop=1.07e+06, mesh_change=0.001%, odd_res=1.150e-10, Rprop=5.806e-07, interior=yes
M=32: kb=3.140675360048, sigma_BEM=7.59123434e-02, sv_min=5.509e-07, drop=1.04e+06, mesh_change=0.000%, odd_res=1.062e-10, Rprop=5.806e-07, interior=yes
M=40: kb=3.140675360901, sigma_BEM=7.59123081e-02, sv_min=5.580e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.104e-10, Rprop=5.806e-07, interior=yes
M=48: kb=3.140675361277, sigma_BEM=7.59122925e-02, sv_min=5.611e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.255e-10, Rprop=5.806e-07, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592488189, sv_min=9.782e-01
kb(M=32)=3.141592441525, sv_min=9.782e-01, drop=1.00e+00, shift=4.666e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592580829, sv_min=9.782e-01
kb(M=32)=3.141592580829, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592642766, sv_min=9.782e-01
kb(M=32)=3.141592642766, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652046, sv_min=9.782e-01
kb(M=32)=3.141592652046, sv_min=9.782e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.140631984100
kb BEM = 3.140675361277
sigma asym = 7.76861734e-02
sigma BEM = 7.59122925e-02
final sv_min(A) = 5.611e-07
final minimum drop = 1.02e+06
final mesh change in sigma = 0.000%
odd-parity residual = 1.255e-10
max first propagating-mode fraction = 5.806e-07
relative sigma error = 2.283%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.070 ===
predicted kb = 3.137900540385
predicted sigma = 1.52264900e-01
predicted kb gap to sqrt(Lambda_2) = 3.692e-03
expected-BIC mesh refinement:
M=16: kb=3.138273756691, sigma_BEM=1.44368380e-01, sv_min=1.959e-06, drop=3.03e+05, mesh_change=--, odd_res=1.651e-10, Rprop=1.273e-06, interior=yes
M=24: kb=3.138273862981, sigma_BEM=1.44366069e-01, sv_min=3.835e-06, drop=1.55e+05, mesh_change=0.002%, odd_res=1.160e-10, Rprop=1.273e-06, interior=yes
M=32: kb=3.138273870562, sigma_BEM=1.44365905e-01, sv_min=2.448e-06, drop=2.43e+05, mesh_change=0.000%, odd_res=1.716e-10, Rprop=1.273e-06, interior=yes
M=40: kb=3.138273873527, sigma_BEM=1.44365840e-01, sv_min=2.001e-06, drop=2.97e+05, mesh_change=0.000%, odd_res=1.247e-10, Rprop=1.273e-06, interior=yes
M=48: kb=3.138273874883, sigma_BEM=1.44365811e-01, sv_min=1.811e-06, drop=3.28e+05, mesh_change=0.000%, odd_res=1.306e-10, Rprop=1.273e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592167792, sv_min=9.696e-01
kb(M=32)=3.141592247660, sv_min=9.696e-01, drop=1.00e+00, shift=7.987e-08, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592491047, sv_min=9.696e-01
kb(M=32)=3.141592491047, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592580829, sv_min=9.696e-01
kb(M=32)=3.141592580829, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652046, sv_min=9.696e-01
kb(M=32)=3.141592652046, sv_min=9.696e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.137900540385
kb BEM = 3.138273874883
sigma asym = 1.52264900e-01
sigma BEM = 1.44365811e-01
final sv_min(A) = 1.811e-06
final minimum drop = 3.28e+05
final mesh change in sigma = 0.000%
odd-parity residual = 1.306e-10
max first propagating-mode fraction = 1.273e-06
relative sigma error = 5.188%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.090 ===
predicted kb = 3.131493237939
predicted sigma = 2.51703202e-01
predicted kb gap to sqrt(Lambda_2) = 1.010e-02
expected-BIC mesh refinement:
M=16: kb=3.133273351258, sigma_BEM=2.28478680e-01, sv_min=7.463e-07, drop=8.45e+05, mesh_change=--, odd_res=1.403e-10, Rprop=2.392e-06, interior=yes
M=24: kb=3.133273637521, sigma_BEM=2.28474755e-01, sv_min=1.827e-07, drop=3.45e+06, mesh_change=0.002%, odd_res=1.518e-10, Rprop=2.392e-06, interior=yes
M=32: kb=3.133273683641, sigma_BEM=2.28474122e-01, sv_min=1.474e-06, drop=4.28e+05, mesh_change=0.000%, odd_res=1.435e-10, Rprop=2.392e-06, interior=yes
M=40: kb=3.133273744923, sigma_BEM=2.28473282e-01, sv_min=9.288e-07, drop=6.79e+05, mesh_change=0.000%, odd_res=1.286e-10, Rprop=2.392e-06, interior=yes
M=48: kb=3.133273751158, sigma_BEM=2.28473196e-01, sv_min=6.664e-07, drop=9.46e+05, mesh_change=0.000%, odd_res=1.342e-10, Rprop=2.392e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141590702613, sv_min=9.636e-01
kb(M=32)=3.141591166236, sv_min=9.637e-01, drop=1.00e+00, shift=4.636e-07, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592567655, sv_min=9.636e-01
kb(M=32)=3.141592567655, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592613929, sv_min=9.636e-01
kb(M=32)=3.141592613929, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592650636, sv_min=9.636e-01
kb(M=32)=3.141592650636, sv_min=9.637e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.131493237939
kb BEM = 3.133273751158
sigma asym = 2.51703202e-01
sigma BEM = 2.28473196e-01
final sv_min(A) = 6.664e-07
final minimum drop = 9.46e+05
final mesh change in sigma = 0.000%
odd-parity residual = 1.342e-10
max first propagating-mode fraction = 2.392e-06
relative sigma error = 9.229%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.110 ===
predicted kb = 3.119010674785
predicted sigma = 3.76001079e-01
predicted kb gap to sqrt(Lambda_2) = 2.258e-02
expected-BIC mesh refinement:
M=16: kb=3.124989777385, sigma_BEM=3.22557425e-01, sv_min=7.528e-07, drop=9.02e+05, mesh_change=--, odd_res=1.852e-10, Rprop=4.084e-06, interior=yes
M=24: kb=3.124990573688, sigma_BEM=3.22549710e-01, sv_min=8.084e-07, drop=8.40e+05, mesh_change=0.002%, odd_res=1.642e-10, Rprop=4.083e-06, interior=yes
M=32: kb=3.124990761648, sigma_BEM=3.22547889e-01, sv_min=7.336e-07, drop=9.26e+05, mesh_change=0.001%, odd_res=1.571e-10, Rprop=4.083e-06, interior=yes
M=40: kb=3.124990826953, sigma_BEM=3.22547256e-01, sv_min=6.671e-07, drop=1.02e+06, mesh_change=0.000%, odd_res=1.640e-10, Rprop=4.083e-06, interior=yes
M=48: kb=3.124990855494, sigma_BEM=3.22546980e-01, sv_min=6.298e-07, drop=1.08e+06, mesh_change=0.000%, odd_res=1.597e-10, Rprop=4.083e-06, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592567655, sv_min=9.621e-01
kb(M=32)=3.141592567655, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592642766, sv_min=9.621e-01
kb(M=32)=3.141592642766, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592650636, sv_min=9.621e-01
kb(M=32)=3.141592650636, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652784, sv_min=9.621e-01
kb(M=32)=3.141592652784, sv_min=9.621e-01, drop=1.00e+00, shift=0.000e+00, odd_res=2.000e+00, Rprop=1.000e+00, persistent BIC=no

--- validation result ---
symmetry conditions (a=0, X even, Y odd, nu~0): PASS
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.119010674785
kb BEM = 3.124990855494
sigma asym = 3.76001079e-01
sigma BEM = 3.22546980e-01
final sv_min(A) = 6.298e-07
final minimum drop = 1.08e+06
final mesh change in sigma = 0.000%
odd-parity residual = 1.597e-10
max first propagating-mode fraction = 4.083e-06
relative sigma error = 14.216%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 4/4
Leading asymptotic sigma within 5.0%: 1/4
epsilon=0.050: BIC=PASS, asymptotic=PASS, sv_min=5.611e-07, mesh_change=0.000%, odd_res=1.255e-10, Rprop=5.806e-07, sigma_error=2.283%
epsilon=0.070: BIC=PASS, asymptotic=FAIL, sv_min=1.811e-06, mesh_change=0.000%, odd_res=1.306e-10, Rprop=1.273e-06, sigma_error=5.188%
epsilon=0.090: BIC=PASS, asymptotic=FAIL, sv_min=6.664e-07, mesh_change=0.000%, odd_res=1.342e-10, Rprop=2.392e-06, sigma_error=9.229%
epsilon=0.110: BIC=PASS, asymptotic=FAIL, sv_min=6.298e-07, mesh_change=0.000%, odd_res=1.597e-10, Rprop=4.083e-06, sigma_error=14.216%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_3_x_symmetry_validation
