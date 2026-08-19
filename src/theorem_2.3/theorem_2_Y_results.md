Sí. El enunciado (iv) también se ve **numéricamente consistente**, pero aquí aparece una diferencia muy interesante respecto al caso de simetría en (x). De hecho, estos resultados nos están diciendo algo que vale la pena investigar más.

Para empezar, las hipótesis geométricas que estás imponiendo son las correctas:

[
X(-t)=-X(t),\qquad Y(-t)=Y(t),
]

y numéricamente obtienes

[
\nu\simeq 1.04\times10^{-17}\approx0.
]

Además estás tomando

[
a_1\simeq-\frac{\beta}{12},
\qquad
a=\varepsilon a_1,
]

de manera que para cada (\varepsilon) el obstáculo se desplaza ligeramente de la línea central. El código confirma exactamente esas condiciones.

## 1. El BIC esperado aparece para los cuatro valores de (\varepsilon)

Los resultados principales son:

[
\begin{array}{c|c|c|c}
\varepsilon & kb_{\mathrm{BEM}} & s_{\min}(A) & \text{drop}\
\hline
0.05 & 3.1406661555 & 6.92\times10^{-7} & 8.19\times10^5\
0.07 & 3.1382404138 & 5.11\times10^{-6} & 1.16\times10^5\
0.09 & 3.1331895640 & 6.61\times10^{-6} & 9.48\times10^4\
0.11 & 3.1248226204 & 1.45\times10^{-5} & 4.66\times10^4
\end{array}
]

Son mínimos muy profundos y la posición del autovalor es extremadamente estable cuando pasas de

[
M=16\rightarrow24\rightarrow32\rightarrow40\rightarrow48.
]

Por ejemplo, para (\varepsilon=0.05), el resultado converge hacia

[
kb_{\mathrm{BEM}}\approx3.140666155494,
]

con prácticamente ningún cambio al aumentar (M).

Por tanto, la parte espectral es muy convincente:

[
\boxed{\text{hay una rama espectral estable dentro de }(\Lambda_1,\Lambda_2).}
]

---

## 2. La unicidad también se mantiene

Al igual que en el caso (iii), aparecen algunas falsas fluctuaciones extremadamente cerca de

[
kb=\pi=\sqrt{\Lambda_2}b,
]

pero mira cómo se comportan:

[
s_{\min}\approx0.97-0.98,
]

[
\text{drop}\approx1,
]

y

[
R_{\mathrm{prop}}\approx1.
]

Es decir, **no son singularidades de (A(k))** y no son BICs.

Para (\varepsilon=0.05), por ejemplo, el código descarta correctamente los cuatro mínimos numéricos adicionales y encuentra solamente el BIC esperado.

Así que, dentro de la resolución numérica estudiada,

[
\boxed{N_{\mathrm{BIC}}=1}
]

para los cuatro valores.

Eso concuerda con el enunciado (iv).

---

# 3. La aproximación asintótica se comporta casi igual que en el caso (iii)

Obtienes:

[
\begin{array}{c|c}
\varepsilon&E_\sigma\
\hline
0.05&2.276%\
0.07&5.179%\
0.09&9.219%\
0.11&14.207%
\end{array}
]

Es casi el mismo patrón que obtuviste para la simetría respecto al eje (x).

Por tanto:

[
\varepsilon=0.05
]

está claramente dentro de nuestra tolerancia del (5%),

mientras que

[
\varepsilon=0.07
]

ya queda ligeramente fuera:

[
5.179%.
]

El resumen del código refleja precisamente esto: existencia/unicidad (4/4), pero precisión asintótica bajo (5%) solo (1/4).

Entonces nuevamente tenemos:

[
\boxed{
\text{la aproximación asintótica pierde precisión antes que el BIC desaparezca.}
}
]

Esto está empezando a ser un patrón muy consistente en tus experimentos.

---

# 4. Pero aquí hay una diferencia MUY interesante respecto al caso (x)

Mira los valores del canal propagante.

Para simetría respecto al eje (x), teníamos aproximadamente:

[
R_{\mathrm{prop}}\sim10^{-6}.
]

Ahora, con simetría respecto al eje (y):

[
\begin{array}{c|c}
\varepsilon&R_{\mathrm{prop}}\
\hline
0.05&1.56\times10^{-4}\
0.07&4.76\times10^{-4}\
0.09&1.15\times10^{-3}\
0.11&2.41\times10^{-3}
\end{array}
]

Sigue siendo pequeño.

Para (\varepsilon=0.05),

[
R_{\mathrm{prop}}\approx0.000156
]

es apenas un

[
0.0156%.
]

Y para (\varepsilon=0.11),

[
R_{\mathrm{prop}}\approx0.002407
]

es alrededor de

[
0.24%.
]

Pero observa que **crece sistemáticamente con (\varepsilon)**.

Esto me parece muy interesante.

---

# 5. Creo que sabemos por qué ocurre

En el caso (iii), la simetría respecto al eje (x) produce de manera exacta la desacoplación del modo propagante.

Ahí tienes una simetría exacta:

[
u(x,-y)=-u(x,y),
]

mientras el primer modo propagante es par. Entonces la ortogonalidad se produce estructuralmente.

Por eso encontrábamos cosas del orden de

[
R_{\mathrm{prop}}\sim10^{-6}.
]

Pero el caso (iv) es diferente.

El teorema no dice simplemente:

> toma cualquier obstáculo simétrico respecto a (y).

Dice que además tiene que cumplirse

[
\boxed{
a=\varepsilon a_1+\mathcal O(\varepsilon^2).
}
]

Y nosotros estamos utilizando solamente

[
\boxed{
a=\varepsilon\left(-\frac{\beta}{12}\right).
}
]

Pero incluso

[
a_1=-\frac{\beta}{12}
]

es en sí una aproximación:

[
a_1
===

-\frac{\beta}{12}

- \mathcal O(\beta^2).
  ]

Entonces tenemos **dos aproximaciones acumuladas**:

[
a_1
===

-\frac{\beta}{12}

- O(\beta^2),
  ]

y posteriormente

[
a
=

\varepsilon a_1

- O(\varepsilon^2).
  ]

Por eso no me sorprende en absoluto que encontremos

[
R_{\mathrm{prop}}>0
]

y que crezca al aumentar (\varepsilon).

---

# 6. Hay que ser cuidadosos con este `PASS`

El código imprime:

```text
placement condition: PASS (residual=0)
```

pero conceptualmente eso no significa que hayamos **validado** la condición

[
a=\varepsilon a_1+O(\varepsilon^2).
]

El residual es cero porque nosotros mismos definimos en el código

[
a:=\varepsilon a_1.
]

Es decir, estamos verificando

[
a-\varepsilon a_1=0
]

porque así lo construimos.

Eso está bien como configuración del experimento, pero en un paper no lo presentaría como una comprobación independiente.

Diría:

> Para la simulación se tomó la aproximación de orden dominante (a=\varepsilon a_1), con (a_1\simeq-\beta/12).

Eso es más exacto.

---

# 7. Y esto abre un experimento MUY bueno

Aquí creo que podemos mejorar significativamente la validación del punto (iv).

En vez de imponer solamente

[
a\_{\mathrm{asym}}
=================

-\frac{\varepsilon\beta}{12},
]

podemos utilizar ese valor como **punto inicial** y buscar numericamente el valor exacto de (a) que minimiza el acoplamiento propagante.

Es decir:

[
a\_{\mathrm{asym}}
=================

-\frac{\varepsilon\beta}{12}.
]

Después buscamos

[
\boxed{
a\_{\mathrm{BIC}}
================

\arg\min*a
R*{\mathrm{prop}}(a)
}
]

en una pequeña vecindad:

[
a
\in
[
a_{\mathrm{asym}}-C\varepsilon^2,,
a_{\mathrm{asym}}+C\varepsilon^2
].
]

Y para cada (a), simultáneamente localizamos el mínimo espectral en (k).

Entonces podríamos obtener algo como:

[
a\_{\mathrm{BIC}}
================

\varepsilon a_1

- C\varepsilon^2+\cdots.
  ]

Ese sería un experimento **mucho más directamente conectado con la afirmación (iv)**.

---

# 8. Y podríamos comprobar la escala del error en (a)

El teorema dice

[
a-\varepsilon a_1
=================

O(\varepsilon^2).
]

Entonces, después de encontrar numéricamente (a\_{\mathrm{BIC}}), podríamos calcular

[
\Delta a
========

a\_{\mathrm{BIC}}-\varepsilon a_1.
]

Si realmente

[
\Delta a=O(\varepsilon^2),
]

entonces

[
\boxed{
\frac{\Delta a}{\varepsilon^2}
}
]

debería tender aproximadamente a una constante cuando

[
\varepsilon\to0.
]

Esa sería una validación numérica muy elegante:

[
\begin{array}{c|c|c}
\varepsilon
&
a_{\mathrm{BIC}}-\varepsilon a_1
&
\dfrac{a_{\mathrm{BIC}}-\varepsilon a_1}{\varepsilon^2}
\
\hline
0.05 & \cdots & C_1\
0.07 & \cdots & C_2\
0.09 & \cdots & C_3
\end{array}
]

y esperaríamos

[
C_1\approx C_2\approx C_3
]

en el régimen pequeño.

Eso sería, desde mi perspectiva, **más fuerte científicamente que simplemente marcar `placement condition: PASS`**.

---

# 9. La paridad también nos dice algo interesante

Tu solución tiene:

[
\text{y-axis even residual}\sim10^{-11}
]

y

[
\text{y-axis odd residual}=2.
]

Eso significa que el estado encontrado pertenece claramente al sector par respecto al eje (y).

Para (\varepsilon=0.05), por ejemplo,

[
R_{\mathrm{even}}\approx2.58\times10^{-11}.
]

Es una confirmación excelente de que estamos siguiendo la rama de simetría correcta.

---

# 10. Mi lectura global de los tres experimentos hasta ahora

Ahora tenemos una historia bastante limpia.

### Teorema 2.1 — modo discreto

[
0<k^2<\Lambda_1
]

Un único modo detectado.

La fórmula pierde precisión al aumentar (\varepsilon).

---

### Teorema 2.3(iii) — BIC, simetría (x)

[
\Lambda_1<k^2<\Lambda_2.
]

Un único BIC.

Simetría exacta produce una cancelación espectacular:

[
R_{\mathrm{prop}}\sim10^{-6}.
]

---

### Teorema 2.3(iv) — BIC, simetría (y)

También un único BIC.

Pero el mecanismo depende del tuning geométrico

[
a=\varepsilon a_1+O(\varepsilon^2).
]

Al utilizar solamente el término líder

[
a=\varepsilon\left(-\frac{\beta}{12}\right),
]

el acoplamiento propagante permanece pequeño pero ya no prácticamente nulo:

[
R_{\mathrm{prop}}
:
1.56\times10^{-4}
\rightarrow
2.41\times10^{-3}.
]

Eso, para mí, no es un problema de los resultados.

**Es probablemente información sobre el término (O(\varepsilon^2)) que estamos ignorando.**

Y eso nos da un siguiente paso bastante natural: en el caso (iv), **optimizar (a) numericamente para encontrar el BIC exacto y comprobar que la corrección respecto a (\varepsilon a_1) escala como (O(\varepsilon^2))**.

---

# Resultados de theorem_2_X_results.md

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python second_theorem_Y.py
0.13058424044339537
0.13058424044372274
=== Focused validation of Theorem 2.3(iv): y-symmetric BIC ===
b = 1.0
beta = 0.1
leading a1 = -beta/12 = -8.333333333333e-03
mu = 1.007140755654e+00
nu = 1.042239055775e-17
X-odd residual = 0.000e+00
Y-even residual = 0.000e+00
Theorem shape symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
sqrt(Lambda_1) b = 1.570796326795
sqrt(Lambda_2) b = 3.141592653590
epsilons = (0.05, 0.07, 0.09, 0.11)
refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
a = epsilon\*a1 = -4.166666666667e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.140622485938
predicted sigma = 7.80692123e-02
predicted kb gap to sqrt(Lambda_2) = 9.702e-04
expected-BIC mesh refinement:
M=16: kb=3.140666142067, sigma_BEM=7.62927596e-02, sv_min=7.244e-07, drop=7.83e+05, mesh_change=--, y-axis parity=even(4.197e-11), Rprop=1.564e-04, interior=yes
M=24: kb=3.140666151929, sigma_BEM=7.62923536e-02, sv_min=7.003e-07, drop=8.10e+05, mesh_change=0.001%, y-axis parity=even(3.322e-11), Rprop=1.564e-04, interior=yes
M=32: kb=3.140666154290, sigma_BEM=7.62922564e-02, sv_min=6.950e-07, drop=8.16e+05, mesh_change=0.000%, y-axis parity=even(5.053e-11), Rprop=1.564e-04, interior=yes
M=40: kb=3.140666155126, sigma_BEM=7.62922220e-02, sv_min=6.932e-07, drop=8.18e+05, mesh_change=0.000%, y-axis parity=even(2.793e-11), Rprop=1.564e-04, interior=yes
M=48: kb=3.140666155494, sigma_BEM=7.62922068e-02, sv_min=6.924e-07, drop=8.19e+05, mesh_change=0.000%, y-axis parity=even(2.584e-11), Rprop=1.564e-04, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592441525, sv_min=9.841e-01
kb(M=32)=3.141592375422, sv_min=9.841e-01, drop=1.00e+00, shift=6.610e-08, y-axis parity=odd(8.789e-08), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592580829, sv_min=9.841e-01
kb(M=32)=3.141592580829, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.734e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592642766, sv_min=9.841e-01
kb(M=32)=3.141592642766, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.469e-07), Rprop=9.999e-01, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652046, sv_min=9.841e-01
kb(M=32)=3.141592652046, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.059e-06), Rprop=9.999e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-4.16666667e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.140622485938
kb BEM = 3.140666155494
sigma asym = 7.80692123e-02
sigma BEM = 7.62922068e-02
final sv_min(A) = 6.924e-07
final minimum drop = 8.19e+05
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 2.584e-11
max first propagating-mode fraction = 1.564e-04
relative sigma error = 2.276%
asymptotic accuracy <= 5.0%: PASS
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: PASS

=== epsilon=0.070 ===
a = epsilon\*a1 = -5.833333333333e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.137864020328
predicted sigma = 1.53015656e-01
predicted kb gap to sqrt(Lambda_2) = 3.729e-03
expected-BIC mesh refinement:
M=16: kb=3.138240305549, sigma_BEM=1.45093714e-01, sv_min=3.417e-06, drop=1.73e+05, mesh_change=--, y-axis parity=even(3.234e-11), Rprop=4.759e-04, interior=yes
M=24: kb=3.138240338969, sigma_BEM=1.45092991e-01, sv_min=3.947e-06, drop=1.50e+05, mesh_change=0.000%, y-axis parity=even(5.922e-11), Rprop=4.759e-04, interior=yes
M=32: kb=3.138240344147, sigma_BEM=1.45092879e-01, sv_min=5.428e-06, drop=1.09e+05, mesh_change=0.000%, y-axis parity=even(1.713e-11), Rprop=4.759e-04, interior=yes
M=40: kb=3.138240413163, sigma_BEM=1.45091386e-01, sv_min=5.385e-06, drop=1.10e+05, mesh_change=0.001%, y-axis parity=even(3.287e-11), Rprop=4.759e-04, interior=yes
M=48: kb=3.138240413756, sigma_BEM=1.45091373e-01, sv_min=5.115e-06, drop=1.16e+05, mesh_change=0.000%, y-axis parity=even(3.469e-11), Rprop=4.759e-04, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 6 total; 5 outside expected-BIC window
checking additional minimum 1/5
kb(M=16)=3.141592247660, sv_min=9.755e-01
kb(M=32)=3.141592325862, sv_min=9.755e-01, drop=1.00e+00, shift=7.820e-08, y-axis parity=odd(6.729e-08), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 2/5
kb(M=16)=3.141592491047, sv_min=9.755e-01
kb(M=32)=3.141592491047, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.307e-07), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 3/5
kb(M=16)=3.141592632870, sv_min=9.755e-01
kb(M=32)=3.141592632870, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(4.822e-07), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 4/5
kb(M=16)=3.141592647935, sv_min=9.755e-01
kb(M=32)=3.141592647935, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.950e-07), Rprop=9.998e-01, persistent BIC=no
checking additional minimum 5/5
kb(M=16)=3.141592652046, sv_min=9.755e-01
kb(M=32)=3.141592652046, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.133e-06), Rprop=9.998e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-5.83333333e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.137864020328
kb BEM = 3.138240413756
sigma asym = 1.53015656e-01
sigma BEM = 1.45091373e-01
final sv_min(A) = 5.115e-06
final minimum drop = 1.16e+05
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 3.469e-11
max first propagating-mode fraction = 4.759e-04
relative sigma error = 5.179%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.090 ===
a = epsilon\*a1 = -7.500000000000e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.131393237609
predicted sigma = 2.52944248e-01
predicted kb gap to sqrt(Lambda_2) = 1.020e-02
expected-BIC mesh refinement:
M=16: kb=3.133189193309, sigma_BEM=2.29629876e-01, sv_min=6.697e-06, drop=9.36e+04, mesh_change=--, y-axis parity=even(7.267e-11), Rprop=1.147e-03, interior=yes
M=24: kb=3.133189445989, sigma_BEM=2.29626428e-01, sv_min=6.561e-06, drop=9.56e+04, mesh_change=0.002%, y-axis parity=even(2.897e-11), Rprop=1.147e-03, interior=yes
M=32: kb=3.133189499024, sigma_BEM=2.29625705e-01, sv_min=6.651e-06, drop=9.43e+04, mesh_change=0.000%, y-axis parity=even(2.945e-11), Rprop=1.147e-03, interior=yes
M=40: kb=3.133189558673, sigma_BEM=2.29624891e-01, sv_min=6.657e-06, drop=9.42e+04, mesh_change=0.000%, y-axis parity=even(3.635e-11), Rprop=1.147e-03, interior=yes
M=48: kb=3.133189564015, sigma_BEM=2.29624818e-01, sv_min=6.610e-06, drop=9.48e+04, mesh_change=0.000%, y-axis parity=even(2.801e-11), Rprop=1.147e-03, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 4 total; 3 outside expected-BIC window
checking additional minimum 1/3
kb(M=16)=3.141592375422, sv_min=9.697e-01
kb(M=32)=3.141592375422, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.257e-08), Rprop=9.995e-01, persistent BIC=no
checking additional minimum 2/3
kb(M=16)=3.141592642766, sv_min=9.697e-01
kb(M=32)=3.141592642766, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(6.442e-07), Rprop=9.995e-01, persistent BIC=no
checking additional minimum 3/3
kb(M=16)=3.141592652784, sv_min=9.697e-01
kb(M=32)=3.141592652784, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.896e-06), Rprop=9.995e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-7.50000000e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.131393237609
kb BEM = 3.133189564015
sigma asym = 2.52944248e-01
sigma BEM = 2.29624818e-01
final sv_min(A) = 6.610e-06
final minimum drop = 9.48e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 2.801e-11
max first propagating-mode fraction = 1.147e-03
relative sigma error = 9.219%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== epsilon=0.110 ===
a = epsilon\*a1 = -9.166666666667e-04
placement condition: PASS (residual=0.000e+00)
predicted kb = 3.118786624544
predicted sigma = 3.77854988e-01
predicted kb gap to sqrt(Lambda_2) = 2.281e-02
expected-BIC mesh refinement:
M=16: kb=3.124821514981, sigma_BEM=3.24183437e-01, sv_min=1.446e-05, drop=4.67e+04, mesh_change=--, y-axis parity=even(3.585e-11), Rprop=2.408e-03, interior=yes
M=24: kb=3.124822313310, sigma_BEM=3.24175741e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.002%, y-axis parity=even(3.300e-11), Rprop=2.407e-03, interior=yes
M=32: kb=3.124822516708, sigma_BEM=3.24173781e-01, sv_min=1.447e-05, drop=4.67e+04, mesh_change=0.001%, y-axis parity=even(3.151e-11), Rprop=2.407e-03, interior=yes
M=40: kb=3.124822588659, sigma_BEM=3.24173087e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.000%, y-axis parity=even(3.662e-11), Rprop=2.407e-03, interior=yes
M=48: kb=3.124822620361, sigma_BEM=3.24172782e-01, sv_min=1.448e-05, drop=4.66e+04, mesh_change=0.000%, y-axis parity=even(2.617e-11), Rprop=2.407e-03, interior=yes
whole-band uniqueness screen: 79 points (M=16)
16/79
32/79
48/79
64/79
79/79
sampled local minima: 5 total; 4 outside expected-BIC window
checking additional minimum 1/4
kb(M=16)=3.141592567655, sv_min=9.682e-01
kb(M=32)=3.141592567655, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(1.752e-07), Rprop=9.983e-01, persistent BIC=no
checking additional minimum 2/4
kb(M=16)=3.141592613929, sv_min=9.682e-01
kb(M=32)=3.141592613929, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(3.686e-07), Rprop=9.983e-01, persistent BIC=no
checking additional minimum 3/4
kb(M=16)=3.141592642766, sv_min=9.682e-01
kb(M=32)=3.141592642766, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(7.500e-07), Rprop=9.983e-01, persistent BIC=no
checking additional minimum 4/4
kb(M=16)=3.141592652046, sv_min=9.682e-01
kb(M=32)=3.141592652046, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, y-axis parity=odd(2.278e-06), Rprop=9.983e-01, persistent BIC=no

--- validation result ---
shape symmetry (X odd, Y even, nu~0): PASS
placement a=epsilon\*a1 (leading order): PASS [a1=-8.33333333e-03, a=-9.16666667e-04]
expected spectral minimum resolved: PASS
non-radiating BIC diagnostics: PASS
additional resolved BICs: 0
uniqueness screen: PASS
exactly one resolved BIC in [Lambda_1,Lambda_2): PASS
kb asymptotic = 3.118786624544
kb BEM = 3.124822620361
sigma asym = 3.77854988e-01
sigma BEM = 3.24172782e-01
final sv_min(A) = 1.448e-05
final minimum drop = 4.66e+04
final mesh change in sigma = 0.000%
y-axis parity diagnostic = even residual 2.617e-11
max first propagating-mode fraction = 2.407e-03
relative sigma error = 14.207%
asymptotic accuracy <= 5.0%: FAIL
EXISTENCE / BIC / UNIQUENESS CHECK: PASS
ASYMPTOTIC APPROXIMATION CHECK: FAIL

=== FINAL SUMMARY ===
Interval screened: Lambda_1 < k^2 < Lambda_2 (excluding tiny numerical layers at the cutoffs).
Unique non-radiating BIC checks passed: 4/4
Leading asymptotic sigma within 5.0%: 1/4
epsilon=0.050: a=-4.166667e-04, BIC=PASS, asymptotic=PASS, sv_min=6.924e-07, mesh_change=0.000%, yparity=even(2.584e-11), Rprop=1.564e-04, sigma_error=2.276%
epsilon=0.070: a=-5.833333e-04, BIC=PASS, asymptotic=FAIL, sv_min=5.115e-06, mesh_change=0.000%, yparity=even(3.469e-11), Rprop=4.759e-04, sigma_error=5.179%
epsilon=0.090: a=-7.500000e-04, BIC=PASS, asymptotic=FAIL, sv_min=6.610e-06, mesh_change=0.000%, yparity=even(2.801e-11), Rprop=1.147e-03, sigma_error=9.219%
epsilon=0.110: a=-9.166667e-04, BIC=PASS, asymptotic=FAIL, sv_min=1.448e-05, mesh_change=0.000%, yparity=even(2.617e-11), Rprop=2.407e-03, sigma_error=14.207%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_3_y_symmetry_validation

# Resultados de theorem_2_Y_results.md

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python second_theorem_Y_2.py
0.13058424044339537
0.13058424044372274
=== Theorem 2.3(iv), version 2: tune a numerically ===
b = 1.0
beta = 0.1
a1 used = -beta/12 = -8.333333333333e-03 [small-beta approximation]
mu = 1.007140755654e+00
nu = 1.042239055775e-17
X-odd residual = 0.000e+00
Y-even residual = 0.000e+00
symmetry conditions: PASS
Lambda_1 = 2.467401100272
Lambda_2 = 9.869604401089
epsilons = (0.05, 0.07, 0.09, 0.11)
a search half-width = 1.0 \* epsilon^2
placement optimization M = (16, 24, 32)
final spectral refinement M = (16, 24, 32, 40, 48)

=== epsilon=0.050 ===
leading a_asym = -4.1666666667e-04
O(eps^2) search window = [-2.9166666667e-03, +2.0833333333e-03]
predicted kb = 3.140622485938
predicted sigma = 7.80692123e-02
tuning a_BIC by minimizing the open-channel projection:
M=16: a_BIC=-2.4860261635e-03, Delta a=-2.069e-03, Delta a/eps^2=-0.827744, kb=3.140666046884, Rprop=3.754e-07, sv_min=9.485e-06, drop=6.38e+04, a-interior=yes
M=24: a_BIC=-2.4860448331e-03, Delta a=-2.069e-03, Delta a/eps^2=-0.827751, kb=3.140666063491, Rprop=3.754e-07, sv_min=1.311e-05, drop=4.62e+04, a-interior=yes
M=32: a_BIC=-2.4860411337e-03, Delta a=-2.069e-03, Delta a/eps^2=-0.827750, kb=3.140666067533, Rprop=3.754e-07, sv_min=1.401e-05, drop=4.32e+04, a-interior=yes
fixed-a_BIC spectral mesh refinement:

    M=16: kb=3.140666046884, sigma_BEM=7.62966778e-02, sv_min=9.485e-06, drop=6.38e+04, mesh_change=--, Rprop=3.755e-07
    M=24: kb=3.140666063491, sigma_BEM=7.62959942e-02, sv_min=1.311e-05, drop=4.62e+04, mesh_change=0.001%, Rprop=3.754e-07
    M=32: kb=3.140666067533, sigma_BEM=7.62958278e-02, sv_min=1.401e-05, drop=4.32e+04, mesh_change=0.000%, Rprop=3.754e-07
    M=40: kb=3.140666068971, sigma_BEM=7.62957686e-02, sv_min=1.434e-05, drop=4.22e+04, mesh_change=0.000%, Rprop=3.754e-07
    M=48: kb=3.140666069606, sigma_BEM=7.62957424e-02, sv_min=1.448e-05, drop=4.18e+04, mesh_change=0.000%, Rprop=3.754e-07

whole-band uniqueness screen at tuned a_BIC: 71 points (M=16)
16/71
32/71
48/71
64/71
71/71
sampled local minima: 5 total; refining all of them
checking minimum 1/5
kb(M=32)=3.140666029622, sv_min=6.538e-06, drop=1.47e+04, shift=1.938e-08, Rprop=3.754e-07, persistent=YES, class=expected BIC
checking minimum 2/5
kb(M=32)=3.141592219547, sv_min=9.841e-01, drop=1.00e+00, shift=7.300e-08, Rprop=9.999e-01, persistent=no, class=other
checking minimum 3/5
kb(M=32)=3.141592581497, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other
checking minimum 4/5
kb(M=32)=3.141592634565, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other
checking minimum 5/5
kb(M=32)=3.141592652576, sv_min=9.841e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other

--- tuned-placement validation result ---
a_asym = epsilon\*a1 = -4.1666666667e-04
a_BIC (numerical) = -2.4860411337e-03
Delta a = -2.069374e-03
Delta a / epsilon^2 = -0.82774979
O(epsilon^2) search consistency: PASS
Rprop at leading a = 1.564e-04
Rprop at tuned a_BIC = 3.754e-07
radiation reduction = 4.17e+02x
spectral BIC resolved = PASS
non-radiating diagnostic = PASS
additional resolved BICs = 0
uniqueness screen = PASS
unique tuned BIC = PASS
kb asymptotic = 3.140622485938
kb BEM = 3.140666069606
sigma asymptotic = 7.80692123e-02
sigma BEM = 7.62957424e-02
final sv_min(A) = 1.448e-05
final minimum drop = 4.18e+04
relative sigma error = 2.272%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.070 ===
leading a_asym = -5.8333333333e-04
O(eps^2) search window = [-5.4833333333e-03, +4.3166666667e-03]
predicted kb = 3.137864020328
predicted sigma = 1.53015656e-01
tuning a_BIC by minimizing the open-channel projection:
M=16: a_BIC=-3.4827953701e-03, Delta a=-2.899e-03, Delta a/eps^2=-0.591727, kb=3.138239541965, Rprop=8.004e-07, sv_min=3.021e-06, drop=2.08e+05, a-interior=yes
M=24: a_BIC=-3.4827580060e-03, Delta a=-2.899e-03, Delta a/eps^2=-0.591719, kb=3.138239604597, Rprop=8.004e-07, sv_min=2.047e-06, drop=3.07e+05, a-interior=yes
M=32: a_BIC=-3.4827543329e-03, Delta a=-2.899e-03, Delta a/eps^2=-0.591719, kb=3.138239619485, Rprop=8.004e-07, sv_min=1.797e-06, drop=3.50e+05, a-interior=yes
fixed-a_BIC spectral mesh refinement:
M=16: kb=3.138239541964, sigma_BEM=1.45110228e-01, sv_min=3.021e-06, drop=2.08e+05, mesh_change=--, Rprop=8.004e-07
M=24: kb=3.138239604597, sigma_BEM=1.45108874e-01, sv_min=2.047e-06, drop=3.07e+05, mesh_change=0.001%, Rprop=8.004e-07
M=32: kb=3.138239619485, sigma_BEM=1.45108552e-01, sv_min=1.797e-06, drop=3.50e+05, mesh_change=0.000%, Rprop=8.004e-07
M=40: kb=3.138239624743, sigma_BEM=1.45108438e-01, sv_min=1.707e-06, drop=3.68e+05, mesh_change=0.000%, Rprop=8.004e-07
M=48: kb=3.138239627059, sigma_BEM=1.45108388e-01, sv_min=1.667e-06, drop=3.77e+05, mesh_change=0.000%, Rprop=8.004e-07
whole-band uniqueness screen at tuned a_BIC: 71 points (M=16)
16/71
32/71
48/71
64/71
71/71
sampled local minima: 6 total; refining all of them
checking minimum 1/6
kb(M=32)=3.138239616102, sv_min=1.281e-06, drop=4.13e+04, shift=8.976e-08, Rprop=8.004e-07, persistent=YES, class=expected BIC
checking minimum 2/6
kb(M=32)=3.141592445061, sv_min=9.755e-01, drop=1.00e+00, shift=5.675e-08, Rprop=9.999e-01, persistent=no, class=other
checking minimum 3/6
kb(M=32)=3.141592560039, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other
checking minimum 4/6
kb(M=32)=3.141592613991, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other
checking minimum 5/6
kb(M=32)=3.141592649199, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other
checking minimum 6/6
kb(M=32)=3.141592653103, sv_min=9.755e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.999e-01, persistent=no, class=other

--- tuned-placement validation result ---
a_asym = epsilon\*a1 = -5.8333333333e-04
a_BIC (numerical) = -3.4827543329e-03
Delta a = -2.899421e-03
Delta a / epsilon^2 = -0.59171857
O(epsilon^2) search consistency: PASS
Rprop at leading a = 4.759e-04
Rprop at tuned a_BIC = 8.004e-07
radiation reduction = 5.95e+02x
spectral BIC resolved = PASS
non-radiating diagnostic = PASS
additional resolved BICs = 0
uniqueness screen = PASS
unique tuned BIC = PASS
kb asymptotic = 3.137864020328
kb BEM = 3.138239627059
sigma asymptotic = 1.53015656e-01
sigma BEM = 1.45108388e-01
final sv_min(A) = 1.667e-06
final minimum drop = 3.77e+05
relative sigma error = 5.168%
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.090 ===
leading a_asym = -7.5000000000e-04
O(eps^2) search window = [-8.8500000000e-03, +7.3500000000e-03]
predicted kb = 3.131393237609
predicted sigma = 2.52944248e-01
tuning a_BIC by minimizing the open-channel projection:
M=16: a_BIC=-4.4795598709e-03, Delta a=-3.730e-03, Delta a/eps^2=-0.460439, kb=3.133186111136, Rprop=1.444e-06, sv_min=4.626e-07, drop=1.44e+06, a-interior=yes
M=24: a_BIC=-4.4795036922e-03, Delta a=-3.730e-03, Delta a/eps^2=-0.460433, kb=3.133186356620, Rprop=1.444e-06, sv_min=2.473e-06, drop=2.69e+05, a-interior=yes
M=32: a_BIC=-4.4794852473e-03, Delta a=-3.729e-03, Delta a/eps^2=-0.460430, kb=3.133186440658, Rprop=1.444e-06, sv_min=1.375e-06, drop=4.83e+05, a-interior=yes
fixed-a_BIC spectral mesh refinement:
M=16: kb=3.133186111142, sigma_BEM=2.29671927e-01, sv_min=4.623e-07, drop=1.44e+06, mesh_change=--, Rprop=1.444e-06
M=24: kb=3.133186356619, sigma_BEM=2.29668578e-01, sv_min=2.473e-06, drop=2.69e+05, mesh_change=0.001%, Rprop=1.444e-06
M=32: kb=3.133186440658, sigma_BEM=2.29667432e-01, sv_min=1.375e-06, drop=4.83e+05, mesh_change=0.000%, Rprop=1.444e-06
M=40: kb=3.133186470602, sigma_BEM=2.29667023e-01, sv_min=9.745e-07, drop=6.81e+05, mesh_change=0.000%, Rprop=1.444e-06
M=48: kb=3.133186483837, sigma_BEM=2.29666843e-01, sv_min=7.954e-07, drop=8.35e+05, mesh_change=0.000%, Rprop=1.444e-06
whole-band uniqueness screen at tuned a_BIC: 71 points (M=16)
16/71
32/71
48/71
64/71
71/71
sampled local minima: 4 total; refining all of them
checking minimum 1/4
kb(M=32)=3.133186447468, sv_min=9.499e-07, drop=1.02e+05, shift=3.392e-07, Rprop=1.444e-06, persistent=YES, class=expected BIC
checking minimum 2/4
kb(M=32)=3.141592613991, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.997e-01, persistent=no, class=other
checking minimum 3/4
kb(M=32)=3.141592649199, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.997e-01, persistent=no, class=other
checking minimum 4/4
kb(M=32)=3.141592653103, sv_min=9.697e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.997e-01, persistent=no, class=other

--- tuned-placement validation result ---
a_asym = epsilon\*a1 = -7.5000000000e-04
a_BIC (numerical) = -4.4794852473e-03
Delta a = -3.729485e-03
Delta a / epsilon^2 = -0.46043028
O(epsilon^2) search consistency: PASS
Rprop at leading a = 1.147e-03
Rprop at tuned a_BIC = 1.444e-06
radiation reduction = 7.95e+02x
spectral BIC resolved = PASS
non-radiating diagnostic = PASS
additional resolved BICs = 0
uniqueness screen = PASS
unique tuned BIC = PASS
kb asymptotic = 3.131393237609
kb BEM = 3.133186483837
sigma asymptotic = 2.52944248e-01
sigma BEM = 2.29666843e-01
final sv_min(A) = 7.954e-07
final minimum drop = 8.35e+05
relative sigma error = 9.203%
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.110 ===
leading a_asym = -9.1666666667e-04
O(eps^2) search window = [-1.3016666667e-02, +1.1183333333e-02]
predicted kb = 3.118786624544
predicted sigma = 3.77854988e-01
tuning a_BIC by minimizing the open-channel projection:
M=16: a_BIC=-5.4763707555e-03, Delta a=-4.560e-03, Delta a/eps^2=-0.376835, kb=3.124812802143, Rprop=2.326e-06, sv_min=8.352e-07, drop=8.52e+05, a-interior=yes
M=24: a_BIC=-5.4762477343e-03, Delta a=-4.560e-03, Delta a/eps^2=-0.376825, kb=3.124813625745, Rprop=2.326e-06, sv_min=7.914e-08, drop=8.99e+06, a-interior=yes
M=32: a_BIC=-5.4762136181e-03, Delta a=-4.560e-03, Delta a/eps^2=-0.376822, kb=3.124813809990, Rprop=2.326e-06, sv_min=1.204e-07, drop=5.91e+06, a-interior=yes
fixed-a_BIC spectral mesh refinement:
M=16: kb=3.124812802138, sigma_BEM=3.24267409e-01, sv_min=8.352e-07, drop=8.52e+05, mesh_change=--, Rprop=2.328e-06
M=24: kb=3.124813625742, sigma_BEM=3.24259472e-01, sv_min=7.910e-08, drop=9.00e+06, mesh_change=0.002%, Rprop=2.326e-06
M=32: kb=3.124813809990, sigma_BEM=3.24257697e-01, sv_min=1.204e-07, drop=5.91e+06, mesh_change=0.001%, Rprop=2.326e-06
M=40: kb=3.124813875166, sigma_BEM=3.24257069e-01, sv_min=1.920e-07, drop=3.71e+06, mesh_change=0.000%, Rprop=2.326e-06
M=48: kb=3.124813903881, sigma_BEM=3.24256792e-01, sv_min=2.237e-07, drop=3.18e+06, mesh_change=0.000%, Rprop=2.326e-06
whole-band uniqueness screen at tuned a_BIC: 71 points (M=16)
16/71
32/71
48/71
64/71
71/71
sampled local minima: 3 total; refining all of them
checking minimum 1/3
kb(M=32)=3.124813804981, sv_min=2.828e-07, drop=2.88e+05, shift=9.684e-07, Rprop=2.326e-06, persistent=YES, class=expected BIC
checking minimum 2/3
kb(M=32)=3.141592613991, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.987e-01, persistent=no, class=other
checking minimum 3/3
kb(M=32)=3.141592653103, sv_min=9.682e-01, drop=1.00e+00, shift=0.000e+00, Rprop=9.987e-01, persistent=no, class=other

--- tuned-placement validation result ---
a_asym = epsilon\*a1 = -9.1666666667e-04
a_BIC (numerical) = -5.4762136181e-03
Delta a = -4.559547e-03
Delta a / epsilon^2 = -0.37682206
O(epsilon^2) search consistency: PASS
Rprop at leading a = 2.407e-03
Rprop at tuned a_BIC = 2.326e-06
radiation reduction = 1.03e+03x
spectral BIC resolved = PASS
non-radiating diagnostic = PASS
additional resolved BICs = 0
uniqueness screen = PASS
unique tuned BIC = PASS
kb asymptotic = 3.118786624544
kb BEM = 3.124813903881
sigma asymptotic = 3.77854988e-01
sigma BEM = 3.24256792e-01
final sv_min(A) = 2.237e-07
final minimum drop = 3.18e+06
relative sigma error = 14.185%
asymptotic accuracy <= 5.0%: FAIL

=== FINAL SUMMARY ===
epsilon=0.050: a_asym=-4.166667e-04, a_BIC=-2.486041e-03, Delta a/eps^2=-0.827750, Rprop 1.564e-04 -> 3.754e-07 (4.17e+02x), BIC=PASS, sigma_error=2.272%
epsilon=0.070: a_asym=-5.833333e-04, a_BIC=-3.482754e-03, Delta a/eps^2=-0.591719, Rprop 4.759e-04 -> 8.004e-07 (5.95e+02x), BIC=PASS, sigma_error=5.168%
epsilon=0.090: a_asym=-7.500000e-04, a_BIC=-4.479485e-03, Delta a/eps^2=-0.460430, Rprop 1.147e-03 -> 1.444e-06 (7.95e+02x), BIC=PASS, sigma_error=9.203%
epsilon=0.110: a_asym=-9.166667e-04, a_BIC=-5.476214e-03, Delta a/eps^2=-0.376822, Rprop 2.407e-03 -> 2.326e-06 (1.03e+03x), BIC=PASS, sigma_error=14.185%

Key quantity for Theorem 2.3(iv):
(a_BIC - epsilon\*a1) / epsilon^2
If the O(epsilon^2) placement correction is resolved, this quantity should
remain O(1) and ideally approach a finite limit as epsilon -> 0.

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_3_y_symmetry_tuned_a_validation
