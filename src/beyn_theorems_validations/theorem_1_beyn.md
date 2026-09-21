Sí. La explicación anterior seguía toda la cadena PDE \(\to\) BEM \(\to A(k)\to\) SVD. Con el código nuevo, **todo eso sigue igual hasta construir \(A(k)\)**; lo que cambia principalmente es cómo encontramos globalmente los valores de \(k\).

# Del problema físico a Beyn + BEM

La cadena completa del código nuevo puede pensarse así:

$$
\boxed{
\text{PDE}
\rightarrow
\text{BEM}
\rightarrow
A(k)u=0
\rightarrow
\text{Beyn}
\rightarrow
k_{\rm candidate}
\rightarrow
\text{SVD local}
\rightarrow
k_{\rm BEM}
\rightarrow
\sigma_{\rm BEM}
\rightarrow
\text{Teorema}
}
$$

## 1. El problema que queremos resolver

Partimos del mismo problema:

$$
-\Delta u=k^2u,
$$

con Dirichlet en las paredes, Neumann sobre el obstáculo y decaimiento cuando \(|x|\to\infty\).

Para el Teorema 2.1 buscamos un eigenvalue discreto en

$$
\boxed{0<k^2<\Lambda_1},
\qquad
\Lambda_1=\left(\frac{\pi}{2b}\right)^2.
$$

El teorema además predice

$$
k^2=\Lambda_1-\sigma^2
$$

y una aproximación asintótica

$$
\sigma_{\rm asym}=C\varepsilon^2.
$$

Esta parte es exactamente la misma que en el código anterior.

---

# 2. BEM convierte la PDE en un problema matricial

Utilizando la Green de la guía y la integral de frontera obtenemos

$$
\frac12u(p)
=
\int_\Gamma
u(q)\frac{\partial G(p,q;k)}{\partial n_q}\,ds_q.
$$

Después de discretizar:

$$
\boxed{
A(k)\mathbf u=0
}
$$

con

$$
\boxed{
A(k)=I-\frac{4\pi}{M}K^w(k).
}
$$

Esto tampoco ha cambiado.

El punto importante es que \(A\) depende de \(k\) de manera **no lineal**, principalmente a través de la Green:

$$
\boxed{A=A(k)}.
$$

Por tanto tenemos un **nonlinear eigenvalue problem**:

$$
\boxed{
A(k)\mathbf u=0,\qquad \mathbf u\neq0.
}
$$

Queremos encontrar los \(k\) para los cuales \(A(k)\) pierde invertibilidad.

---

# 3. ¿Qué hacía el código anterior?

El código anterior recorría valores reales de \(k\) y calculaba

$$
s_{\min}(A(k)).
$$

Después buscaba valles donde

$$
\boxed{s_{\min}(A(k))\approx0}.
$$

Eso funciona, pero para demostrar que no había otros modos tuvimos que hacer un **whole-band exhaustive screen**. Esa era precisamente la parte global del método anterior.

El código nuevo reemplaza esa búsqueda global por **Beyn**.

---

# 4. Entra el método de Beyn

En lugar de recorrer todos los \(k\) reales, extendemos

$$
A(k)\longrightarrow A(z),
\qquad z\in\mathbb C.
$$

Escogemos un contorno cerrado \(\Gamma\) en el plano complejo que encierra la región espectral que queremos estudiar:

$$
0<k^2<\Lambda_1.
$$

Conceptualmente:

$$
\boxed{
\Gamma
\text{ encierra todos los eigenvalues que queremos contar.}
}
$$

El contorno se mantiene separado del cutoff \(\sqrt{\Lambda_1}\), porque allí aparece un branch point de la Green.

---

# 5. Beyn construye dos integrales matriciales

Tomamos una pequeña matriz de prueba \(V\) y calculamos

$$
\boxed{
S_0=
\frac{1}{2\pi i}
\oint_\Gamma
A(z)^{-1}V\,dz
}
$$

y

$$
\boxed{
S_1=
\frac{1}{2\pi i}
\oint_\Gamma
z\,A(z)^{-1}V\,dz.
}
$$

Esta es la parte fundamental del método.

Para cada punto \(z_j\) del contorno el código:

$$
z_j
\rightarrow
A(z_j)
\rightarrow
A(z_j)X_j=V.
$$

Es decir, **no buscamos singularidades directamente sobre el contorno**: resolvemos sistemas lineales en muchos puntos alrededor de él.

---

# 6. ¿Cómo sabemos cuántos eigenvalues hay?

Hacemos SVD de \(S_0\):

$$
\boxed{
S_0=U\Sigma W^*.
}
$$

Si

$$
\Sigma=
\operatorname{diag}(s_1,s_2,\ldots),
$$

un salto grande como

$$
s_1\gg s_2
$$

indica

$$
\boxed{\operatorname{rank}(S_0)=1}.
$$

Bajo las hipótesis del método de Beyn y con una matriz de probing suficiente, ese rango representa el número de eigenvalues encerrados, contando multiplicidad.

Esto es muy importante porque ahora la evidencia de unicidad no viene simplemente de:

> “muestreamos la banda y solo vimos un valle”,

sino de una propiedad espectral global del problema.

En nuestros cuatro experimentos el código obtuvo

$$
\boxed{\operatorname{rank}(S_0)=1}.
$$

---

# 7. ¿Cómo obtiene Beyn la posición del eigenvalue?

Si el rango detectado es \(r\), conservamos

$$
U_r,\qquad
\Sigma_r,\qquad
W_r.
$$

Luego construimos el pequeño problema

$$
\boxed{
B=
U_r^*S_1W_r\Sigma_r^{-1}.
}
$$

Los eigenvalues de \(B\),

$$
B y=\lambda y,
$$

son aproximaciones a los eigenvalues no lineales originales:

$$
\boxed{
\lambda\approx k.
}
$$

Así obtenemos el `raw Beyn eigenvalue`.

---

# 8. ¿Por qué repetimos \(N_q=96,192,384,\ldots\)?

Las integrales de contorno se calculan numéricamente. Por eso usamos

$$
N_q=96,\;192,\;384
$$

puntos y observamos si:

$$
\operatorname{rank}(S_0)
$$

permanece estable y si

$$
k_{\rm Beyn}^{(N_q)}
$$

converge.

Normalmente \(384\) fue suficiente.

Pero para \(\varepsilon=0.05\), el eigenvalue está extremadamente cerca del cutoff. El código observó

$$
1.5709633
\rightarrow
1.5708350
\rightarrow
1.5707912
$$

y todavía no tenía un candidato dentro del contorno. Entonces escaló automáticamente a

$$
N_q=768,
$$

obteniendo

$$
\boxed{
k_{\rm Beyn}=1.570775857305.
}
$$

Esto es la parte **adaptive Beyn** del código.

---

# 9. Beyn descubre; SVD certifica

Aquí hay una distinción importante.

No tomamos directamente

$$
k_{\rm Beyn}
$$

como resultado final.

Lo usamos como una semilla para una búsqueda local de

$$
\boxed{
\min_k s_{\min}(A(k)).
}
$$

Por ejemplo, para \(\varepsilon=0.05\):

$$
k_{\rm Beyn}=1.570775857305
$$

pero el refinamiento SVD encuentra

$$
\boxed{
k_{\rm BEM}=1.570769183935.
}
$$

Además repetimos el cálculo para

$$
M=16,24,32,40,48
$$

y verificamos que la posición converge.

Por eso la arquitectura realmente es

$$
\boxed{
\underbrace{\text{Beyn}}_{\text{global discovery/counting}}
+
\underbrace{\text{SVD+BEM}}_{\text{local certification}}.
}
$$

---

# 10. ¿Qué papel tienen `near-contour` y Aitken?

Son mecanismos de robustez.

Si una aproximación de Beyn queda ligeramente fuera del contorno pero parece converger hacia él, se puede marcar como

$$
\texttt{near-contour-seed}.
$$

No significa “eigenvalue confirmado”.

Igualmente, Aitken puede extrapolar una secuencia

$$
k_{96},k_{192},k_{384}
$$

para estimar hacia dónde converge.

Pero ninguno de ellos reemplaza la certificación:

$$
\boxed{
\text{la posición final sigue viniendo del SVD local}.
}
$$

En la última corrida ni siquiera hizo falta ese fallback: los cuatro modos terminaron usando `strict-beyn` como semilla.

---

# 11. Finalmente volvemos al teorema

Una vez encontrado

$$
k_{\rm BEM},
$$

calculamos

$$
\boxed{
\sigma_{\rm BEM}
=
\sqrt{\Lambda_1-k_{\rm BEM}^2}.
}
$$

El teorema nos había dado independientemente

$$
\sigma_{\rm asym}.
$$

Entonces calculamos

$$
\boxed{
E_\sigma
=
\frac{
|\sigma_{\rm BEM}-\sigma_{\rm asym}|
}{
|\sigma_{\rm asym}|
}.
}
$$

Los resultados finales fueron:

$$
\begin{array}{c|c|c}
\varepsilon & \operatorname{rank}(S_0)&E_\sigma\\
\hline
0.05&1&1.090\%\\
0.07&1&2.436\%\\
0.09&1&4.375\%\\
0.11&1&6.876\%
\end{array}
$$

Por tanto:

$$
\boxed{\text{un único modo detectado en }4/4}
$$

y

$$
\boxed{\text{aproximación asintótica dentro del 5\% en }3/4.}
$$

## El mapa mental que conservaría

El código nuevo puede resumirse casi completamente con:

$$
\boxed{
\begin{aligned}
&\text{BEM:} &&
A(k)\mathbf u=0,
\\[1mm]
&\text{Beyn:} &&
S_0=\frac{1}{2\pi i}\oint_\Gamma A(z)^{-1}V\,dz,
\\
&&&
S_1=\frac{1}{2\pi i}\oint_\Gamma zA(z)^{-1}V\,dz,
\\[1mm]
&\text{Conteo:} &&
\operatorname{rank}(S_0)=N_{\rm eig},
\\[1mm]
&\text{Localización:} &&
\operatorname{eig}
\left(
U_r^*S_1W_r\Sigma_r^{-1}
\right)
\approx k,
\\[1mm]
&\text{Certificación:} &&
k_{\rm BEM}
=
\arg\min_k s_{\min}(A(k)),
\\[1mm]
&\text{Validación:} &&
\sigma_{\rm BEM}
=
\sqrt{\Lambda_1-k_{\rm BEM}^2}
\overset{?}{\approx}
\sigma_{\rm asym}.
\end{aligned}}
$$

La diferencia conceptual más importante frente al código anterior es:

$$
\boxed{
\text{antes: búsqueda global por muestreo de }s_{\min}(A(k))
}
$$

$$
\boxed{
\text{ahora: conteo/localización global por Beyn
+ certificación local por SVD}.
}
$$

Y algo que vale la pena destacar en la presentación: **la fórmula asintótica no se usa para localizar el eigenvalue con Beyn**. Beyn examina globalmente el intervalo encerrado por el contorno; solo al final comparamos el resultado obtenido con la predicción del teorema. Eso hace la validación más independiente.

---

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python first_theorem_4.py
=== Theorem 2.1 validation with Beyn contour discovery (v2) ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
epsilons = (0.05, 0.07, 0.09, 0.11)
Beyn ellipse: real endpoints=[0.00010000, 1.57078633], center=0.78544316, rx=0.78534316, ry=3.000e-02
Beyn settings: M=24, Nq-levels=(96, 192, 384), probe_dim=8
rank diagnostics: rel_tol=1.0e-08, gap_threshold=1.0e+03
local verification M = (16, 24, 32, 40, 48)

Checking complex-kb support required by Beyn...
complex-kb assembly: PASS

=== epsilon=0.050 ===
predicted kb = 1.570768582119
predicted sigma = 9.33604312e-03
predicted cutoff gap = 2.774e-05
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.719e+04
Beyn estimated enclosed rank = 1
singular values of S0: [2.825e-04, 4.939e-09, 2.411e-12, 1.946e-12, 1.675e-12, 1.460e-12, 1.204e-12, 9.691e-13]
max contour linear-solve residual = 3.779e-16
raw Beyn eigenvalues:
[0] +1.570963315921 -2.186e-09i outside-contour/band
Nq=96: rank=1, candidates=0, S0-change=--, candidate-shift=--
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.688e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.285e-04, 1.353e-09, 1.919e-12, 1.732e-12, 1.499e-12, 1.274e-12, 9.016e-13, 7.497e-13]
max contour linear-solve residual = 4.364e-16
raw Beyn eigenvalues:
[0] +1.570834970811 -6.507e-10i outside-contour/band
Nq=192: rank=1, candidates=0, S0-change=2.362e-01, candidate-shift=--
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.058e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.022e-04, 3.997e-10, 1.762e-12, 1.593e-12, 1.296e-12, 1.215e-12, 1.119e-12, 9.452e-13]
max contour linear-solve residual = 5.810e-16
raw Beyn eigenvalues:
[0] +1.570791155740 -1.601e-10i outside-contour/band
Nq=384: rank=1, candidates=0, S0-change=1.302e-01, candidate-shift=--
final clustered near-real Beyn candidates = 0

--- Beyn-v2 + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final near-real candidates = 0
locally resolved modes = 0
one-mode count supported = FAIL
kb asymptotic = 1.570768582119
kb BEM = nan
sigma asym = 9.33604312e-03
sigma BEM = nan
final sigma_min(A) = nan
final minimum drop = nan
final mesh change in sigma = nan%
relative sigma error = nan%
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.070 ===
predicted kb = 1.570689740172
predicted sigma = 1.82986445e-02
predicted cutoff gap = 1.066e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=6.119e+04
Beyn estimated enclosed rank = 1
singular values of S0: [5.862e-04, 9.580e-09, 3.160e-12, 2.506e-12, 2.434e-12, 2.079e-12, 1.607e-12, 1.439e-12]
max contour linear-solve residual = 4.050e-16
raw Beyn eigenvalues:
[0] +1.570823960741 -2.099e-09i outside-contour/band
Nq=96: rank=1, candidates=0, S0-change=--, candidate-shift=--
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.195e+05
Beyn estimated enclosed rank = 1
singular values of S0: [5.763e-04, 2.626e-09, 2.702e-12, 2.194e-12, 1.954e-12, 1.830e-12, 1.431e-12, 1.133e-12]
max contour linear-solve residual = 6.219e-16
raw Beyn eigenvalues:
[0] +1.570730866534 +5.850e-11i real-candidate
Nq=192: rank=1, candidates=1, S0-change=1.716e-02, candidate-shift=--
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.232e+05
Beyn estimated enclosed rank = 1
singular values of S0: [6.391e-04, 7.764e-10, 2.464e-12, 2.334e-12, 1.797e-12, 1.505e-12, 1.315e-12, 9.538e-13]
max contour linear-solve residual = 9.130e-16
raw Beyn eigenvalues:
[0] +1.570704461001 +2.110e-10i real-candidate
Nq=384: rank=1, candidates=1, S0-change=9.829e-02, candidate-shift=2.641e-05
final clustered near-real Beyn candidates = 1
candidate 1: Beyn kb=1.570704461001 +2.110e-10i, local bracket=[1.568704461001, 1.570786326795]
M=16: kb=1.570694867809, sigma_BEM=1.78530812e-02, sv_min=1.668e-06, drop=6.51e+04, mesh_change=--, interior=yes
M=24: kb=1.570694869053, sigma_BEM=1.78529718e-02, sv_min=2.158e-06, drop=5.03e+04, mesh_change=0.001%, interior=yes
M=32: kb=1.570694869347, sigma_BEM=1.78529460e-02, sv_min=2.272e-06, drop=4.78e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570694869451, sigma_BEM=1.78529368e-02, sv_min=2.313e-06, drop=4.69e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570694869497, sigma_BEM=1.78529328e-02, sv_min=2.331e-06, drop=4.66e+04, mesh_change=0.000%, interior=yes
resolved=YES

--- Beyn-v2 + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final near-real candidates = 1
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570689740172
kb BEM = 1.570694869497
sigma asym = 1.82986445e-02
sigma BEM = 1.78529328e-02
final sigma_min(A) = 2.331e-06
final minimum drop = 4.66e+04
final mesh change in sigma = 0.000%
relative sigma error = 2.436%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.090 ===
predicted kb = 1.570505049848
predicted sigma = 3.02487797e-02
predicted cutoff gap = 2.913e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.144e+04
Beyn estimated enclosed rank = 1
singular values of S0: [1.278e-03, 1.569e-08, 4.657e-12, 3.981e-12, 3.622e-12, 2.643e-12, 2.399e-12, 1.756e-12]
max contour linear-solve residual = 5.174e-16
raw Beyn eigenvalues:
[0] +1.570604821989 -1.567e-09i real-candidate
Nq=96: rank=1, candidates=1, S0-change=--, candidate-shift=--
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=3.382e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.454e-03, 4.300e-09, 3.737e-12, 3.154e-12, 2.598e-12, 2.342e-12, 1.972e-12, 1.692e-12]
max contour linear-solve residual = 8.628e-16
raw Beyn eigenvalues:
[0] +1.570548000754 -4.111e-10i real-candidate
Nq=192: rank=1, candidates=1, S0-change=1.214e-01, candidate-shift=5.682e-05
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.260e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.601e-03, 1.271e-09, 3.742e-12, 3.386e-12, 2.903e-12, 2.657e-12, 2.186e-12, 1.895e-12]
max contour linear-solve residual = 1.385e-15
raw Beyn eigenvalues:
[0] +1.570534810125 -3.545e-10i real-candidate
Nq=384: rank=1, candidates=1, S0-change=9.179e-02, candidate-shift=1.319e-05
final clustered near-real Beyn candidates = 1
candidate 1: Beyn kb=1.570534810125 -3.545e-10i, local bracket=[1.568534810125, 1.570786326795]
M=16: kb=1.570529958290, sigma_BEM=2.89266380e-02, sv_min=3.453e-06, drop=4.75e+04, mesh_change=--, interior=yes
M=24: kb=1.570529959047, sigma_BEM=2.89265969e-02, sv_min=4.141e-06, drop=3.96e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570529959225, sigma_BEM=2.89265873e-02, sv_min=4.307e-06, drop=3.81e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570529959287, sigma_BEM=2.89265839e-02, sv_min=4.366e-06, drop=3.75e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570529982676, sigma_BEM=2.89253140e-02, sv_min=4.345e-06, drop=3.77e+04, mesh_change=0.004%, interior=yes
resolved=YES

--- Beyn-v2 + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final near-real candidates = 1
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570505049848
kb BEM = 1.570529982676
sigma asym = 3.02487797e-02
sigma BEM = 2.89253140e-02
final sigma_min(A) = 4.345e-06
final minimum drop = 3.77e+04
final mesh change in sigma = 0.004%
relative sigma error = 4.375%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.110 ===
predicted kb = 1.570146262336
predicted sigma = 4.51864487e-02
predicted cutoff gap = 6.501e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.181e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.744e-03, 2.323e-08, 4.685e-12, 4.382e-12, 3.527e-12, 3.118e-12, 2.955e-12, 2.141e-12]
max contour linear-solve residual = 6.182e-16
raw Beyn eigenvalues:
[0] +1.570274341561 -1.955e-09i real-candidate
Nq=96: rank=1, candidates=1, S0-change=--, candidate-shift=--
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.473e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.848e-03, 6.368e-09, 3.934e-12, 3.481e-12, 3.111e-12, 2.573e-12, 2.288e-12, 1.982e-12]
max contour linear-solve residual = 9.857e-16
raw Beyn eigenvalues:
[0] +1.570243629925 -3.585e-10i real-candidate
Nq=192: rank=1, candidates=1, S0-change=3.658e-02, candidate-shift=3.071e-05
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.505e+06
Beyn estimated enclosed rank = 1
singular values of S0: [2.832e-03, 1.882e-09, 4.634e-12, 4.044e-12, 3.593e-12, 3.053e-12, 2.442e-12, 1.875e-12]
max contour linear-solve residual = 1.421e-15
raw Beyn eigenvalues:
[0] +1.570235883656 -1.913e-11i real-candidate
Nq=384: rank=1, candidates=1, S0-change=5.667e-03, candidate-shift=7.746e-06
final clustered near-real Beyn candidates = 1
candidate 1: Beyn kb=1.570235883656 -1.913e-11i, local bracket=[1.568235883656, 1.570786326795]
M=16: kb=1.570232590865, sigma_BEM=4.20798153e-02, sv_min=1.660e-06, drop=1.29e+05, mesh_change=--, interior=yes
M=24: kb=1.570232598950, sigma_BEM=4.20795136e-02, sv_min=1.658e-06, drop=1.30e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.570232600882, sigma_BEM=4.20794415e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.570232601567, sigma_BEM=4.20794160e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.570232601869, sigma_BEM=4.20794047e-02, sv_min=1.656e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
resolved=YES

--- Beyn-v2 + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final near-real candidates = 1
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570146262336
kb BEM = 1.570232601869
sigma asym = 4.51864487e-02
sigma BEM = 4.20794047e-02
final sigma_min(A) = 1.656e-06
final minimum drop = 1.30e+05
final mesh change in sigma = 0.000%
relative sigma error = 6.876%
asymptotic accuracy <= 5.0%: FAIL

=== FINAL SUMMARY ===
Global discovery: Beyn contour method with Nq convergence study
Local certification: SVD minima + BEM mesh refinement
Interval enclosed: 0 < k^2 < Lambda_1, excluding endpoint margins
Stable Beyn rank=1 + one locally resolved mode: 3/4
Leading asymptotic sigma within 5.0%: 2/4
epsilon=0.050: Nq=384, rank=1, rank_stable=yes, near-real=0, resolved=0, one-mode=FAIL, sigma_error=--
epsilon=0.070: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=2.436%
epsilon=0.090: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=4.375%
epsilon=0.110: Nq=384, rank=1, rank_stable=yes, near-real=1, resolved=1, one-mode=PASS, sigma_error=6.876%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_1_beyn_v2_validation

---

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 % python first_theorem_beyn_2.py
=== Theorem 2.1 validation with adaptive Beyn contour discovery (v3) ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
epsilons = (0.05, 0.07, 0.09, 0.11)
Beyn ellipse: real endpoints=[0.00010000, 1.57078633], center=0.78544316, rx=0.78534316, ry=3.000e-02
Beyn settings: M=24, base-Nq=(96, 192, 384), adaptive-Nq=(768, 1536), probe_dim=8
rank diagnostics: rel_tol=1.0e-08, gap_threshold=1.0e+03
local verification M = (16, 24, 32, 40, 48)

Checking complex-kb support required by Beyn...
complex-kb assembly: PASS

=== epsilon=0.050 ===
predicted kb = 1.570768582119
predicted sigma = 9.33604312e-03
predicted cutoff gap = 2.774e-05
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.719e+04
Beyn estimated enclosed rank = 1
singular values of S0: [2.825e-04, 4.939e-09, 2.411e-12, 1.946e-12, 1.675e-12, 1.460e-12, 1.204e-12, 9.691e-13]
max contour linear-solve residual = 3.779e-16
raw Beyn eigenvalues:
[0] +1.570963315921 -2.186e-09i outside-contour/band
Nq=96: rank=1, strict=0, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570963315921-2.19e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.688e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.285e-04, 1.353e-09, 1.919e-12, 1.732e-12, 1.499e-12, 1.274e-12, 9.016e-13, 7.497e-13]
max contour linear-solve residual = 4.364e-16
raw Beyn eigenvalues:
[0] +1.570834970811 -6.507e-10i outside-contour/band
Nq=192: rank=1, strict=0, near-contour=0, S0-change=2.362e-01, candidate-shift=--, raw=1.570834970811-6.51e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.058e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.022e-04, 3.997e-10, 1.762e-12, 1.593e-12, 1.296e-12, 1.215e-12, 1.119e-12, 9.452e-13]
max contour linear-solve residual = 5.810e-16
raw Beyn eigenvalues:
[0] +1.570791155740 -1.601e-10i near-contour-seed
Nq=384: rank=1, strict=0, near-contour=1, S0-change=1.302e-01, candidate-shift=--, raw=1.570791155740-1.60e-10i
continuing adaptive quadrature: no strict candidate
adaptive Beyn escalation -> Nq=768
Beyn contour discovery: M=24, quadrature=768, probe_dim=8
contour node 128/768
contour node 256/768
contour node 384/768
contour node 512/768
contour node 640/768
contour node 768/768
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.661e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.981e-04, 1.193e-10, 1.817e-12, 1.656e-12, 1.340e-12, 1.288e-12, 1.023e-12, 8.857e-13]
max contour linear-solve residual = 1.031e-15
raw Beyn eigenvalues:
[0] +1.570775857305 +1.499e-10i real-candidate
Nq=768: rank=1, strict=1, near-contour=0, S0-change=2.049e-02, candidate-shift=--, raw=1.570775857305+1.50e-10i
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570775857305 +1.499e-10i, local bracket=[1.568775857305, 1.570786326795]
M=16: kb=1.570769183776, sigma_BEM=9.23426071e-03, sv_min=1.221e-05, drop=3.63e+03, mesh_change=--, interior=yes
M=24: kb=1.570769183892, sigma_BEM=9.23424097e-03, sv_min=1.227e-05, drop=3.61e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570769183920, sigma_BEM=9.23423625e-03, sv_min=1.229e-05, drop=3.61e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570769183930, sigma_BEM=9.23423456e-03, sv_min=1.230e-05, drop=3.61e+03, mesh_change=0.000%, interior=yes
M=48: kb=1.570769183935, sigma_BEM=9.23423381e-03, sv_min=1.230e-05, drop=3.60e+03, mesh_change=0.000%, interior=yes
resolved=YES
NOTE: final strict Beyn candidate position is still moving by more than 2.0e-04 between the last two Nq levels. Local SVD refinement remains the trusted position.

--- Beyn-v3 adaptive + local-SVD validation result ---
final Beyn quadrature Nq = 768
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570767650092
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570768582119
kb BEM = 1.570769183935
sigma asym = 9.33604312e-03
sigma BEM = 9.23423381e-03
final sigma_min(A) = 1.230e-05
final minimum drop = 3.60e+03
final mesh change in sigma = 0.000%
relative sigma error = 1.090%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.070 ===
predicted kb = 1.570689740172
predicted sigma = 1.82986445e-02
predicted cutoff gap = 1.066e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=6.119e+04
Beyn estimated enclosed rank = 1
singular values of S0: [5.862e-04, 9.580e-09, 3.160e-12, 2.506e-12, 2.434e-12, 2.079e-12, 1.607e-12, 1.439e-12]
max contour linear-solve residual = 4.050e-16
raw Beyn eigenvalues:
[0] +1.570823960741 -2.099e-09i outside-contour/band
Nq=96: rank=1, strict=0, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570823960741-2.10e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.195e+05
Beyn estimated enclosed rank = 1
singular values of S0: [5.763e-04, 2.626e-09, 2.702e-12, 2.194e-12, 1.954e-12, 1.830e-12, 1.431e-12, 1.133e-12]
max contour linear-solve residual = 6.219e-16
raw Beyn eigenvalues:
[0] +1.570730866534 +5.850e-11i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=1.716e-02, candidate-shift=--, raw=1.570730866534+5.85e-11i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.232e+05
Beyn estimated enclosed rank = 1
singular values of S0: [6.391e-04, 7.764e-10, 2.464e-12, 2.334e-12, 1.797e-12, 1.505e-12, 1.315e-12, 9.538e-13]
max contour linear-solve residual = 9.130e-16
raw Beyn eigenvalues:
[0] +1.570704461001 +2.110e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=9.829e-02, candidate-shift=2.641e-05, raw=1.570704461001+2.11e-10i
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570704461001 +2.110e-10i, local bracket=[1.568704461001, 1.570786326795]
M=16: kb=1.570694867809, sigma_BEM=1.78530812e-02, sv_min=1.668e-06, drop=6.51e+04, mesh_change=--, interior=yes
M=24: kb=1.570694869053, sigma_BEM=1.78529718e-02, sv_min=2.158e-06, drop=5.03e+04, mesh_change=0.001%, interior=yes
M=32: kb=1.570694869347, sigma_BEM=1.78529460e-02, sv_min=2.272e-06, drop=4.78e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570694869451, sigma_BEM=1.78529368e-02, sv_min=2.313e-06, drop=4.69e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570694869497, sigma_BEM=1.78529328e-02, sv_min=2.331e-06, drop=4.66e+04, mesh_change=0.000%, interior=yes
resolved=YES

--- Beyn-v3 adaptive + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570694005670
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570689740172
kb BEM = 1.570694869497
sigma asym = 1.82986445e-02
sigma BEM = 1.78529328e-02
final sigma_min(A) = 2.331e-06
final minimum drop = 4.66e+04
final mesh change in sigma = 0.000%
relative sigma error = 2.436%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.090 ===
predicted kb = 1.570505049848
predicted sigma = 3.02487797e-02
predicted cutoff gap = 2.913e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.144e+04
Beyn estimated enclosed rank = 1
singular values of S0: [1.278e-03, 1.569e-08, 4.657e-12, 3.981e-12, 3.622e-12, 2.643e-12, 2.399e-12, 1.756e-12]
max contour linear-solve residual = 5.174e-16
raw Beyn eigenvalues:
[0] +1.570604821989 -1.567e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570604821989-1.57e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=3.382e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.454e-03, 4.300e-09, 3.737e-12, 3.154e-12, 2.598e-12, 2.342e-12, 1.972e-12, 1.692e-12]
max contour linear-solve residual = 8.628e-16
raw Beyn eigenvalues:
[0] +1.570548000754 -4.111e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=1.214e-01, candidate-shift=5.682e-05, raw=1.570548000754-4.11e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.260e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.601e-03, 1.271e-09, 3.742e-12, 3.386e-12, 2.903e-12, 2.657e-12, 2.186e-12, 1.895e-12]
max contour linear-solve residual = 1.385e-15
raw Beyn eigenvalues:
[0] +1.570534810125 -3.545e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=9.179e-02, candidate-shift=1.319e-05, raw=1.570534810125-3.55e-10i
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570534810125 -3.545e-10i, local bracket=[1.568534810125, 1.570786326795]
M=16: kb=1.570529958290, sigma_BEM=2.89266380e-02, sv_min=3.453e-06, drop=4.75e+04, mesh_change=--, interior=yes
M=24: kb=1.570529959047, sigma_BEM=2.89265969e-02, sv_min=4.141e-06, drop=3.96e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570529959225, sigma_BEM=2.89265873e-02, sv_min=4.307e-06, drop=3.81e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570529959287, sigma_BEM=2.89265839e-02, sv_min=4.366e-06, drop=3.75e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570529982676, sigma_BEM=2.89253140e-02, sv_min=4.345e-06, drop=3.77e+04, mesh_change=0.004%, interior=yes
resolved=YES

--- Beyn-v3 adaptive + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570530822266
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570505049848
kb BEM = 1.570529982676
sigma asym = 3.02487797e-02
sigma BEM = 2.89253140e-02
final sigma_min(A) = 4.345e-06
final minimum drop = 3.77e+04
final mesh change in sigma = 0.004%
relative sigma error = 4.375%
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.110 ===
predicted kb = 1.570146262336
predicted sigma = 4.51864487e-02
predicted cutoff gap = 6.501e-04
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.181e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.744e-03, 2.323e-08, 4.685e-12, 4.382e-12, 3.527e-12, 3.118e-12, 2.955e-12, 2.141e-12]
max contour linear-solve residual = 6.182e-16
raw Beyn eigenvalues:
[0] +1.570274341561 -1.955e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570274341561-1.96e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.473e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.848e-03, 6.368e-09, 3.934e-12, 3.481e-12, 3.111e-12, 2.573e-12, 2.288e-12, 1.982e-12]
max contour linear-solve residual = 9.857e-16
raw Beyn eigenvalues:
[0] +1.570243629925 -3.585e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=3.658e-02, candidate-shift=3.071e-05, raw=1.570243629925-3.59e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.505e+06
Beyn estimated enclosed rank = 1
singular values of S0: [2.832e-03, 1.882e-09, 4.634e-12, 4.044e-12, 3.593e-12, 3.053e-12, 2.442e-12, 1.875e-12]
max contour linear-solve residual = 1.421e-15
raw Beyn eigenvalues:
[0] +1.570235883656 -1.913e-11i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=5.667e-03, candidate-shift=7.746e-06, raw=1.570235883656-1.91e-11i
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570235883656 -1.913e-11i, local bracket=[1.568235883656, 1.570786326795]
M=16: kb=1.570232590865, sigma_BEM=4.20798153e-02, sv_min=1.660e-06, drop=1.29e+05, mesh_change=--, interior=yes
M=24: kb=1.570232598950, sigma_BEM=4.20795136e-02, sv_min=1.658e-06, drop=1.30e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.570232600882, sigma_BEM=4.20794415e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.570232601567, sigma_BEM=4.20794160e-02, sv_min=1.657e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.570232601869, sigma_BEM=4.20794047e-02, sv_min=1.656e-06, drop=1.30e+05, mesh_change=0.000%, interior=yes
resolved=YES

--- Beyn-v3 adaptive + local-SVD validation result ---
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570233270822
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570146262336
kb BEM = 1.570232601869
sigma asym = 4.51864487e-02
sigma BEM = 4.20794047e-02
final sigma_min(A) = 1.656e-06
final minimum drop = 1.30e+05
final mesh change in sigma = 0.000%
relative sigma error = 6.876%
asymptotic accuracy <= 5.0%: FAIL

=== FINAL SUMMARY ===
Global discovery: adaptive Beyn contour method with Nq convergence study
Local certification: SVD minima + BEM mesh refinement
Interval enclosed: 0 < k^2 < Lambda_1, excluding endpoint margins
Stable Beyn rank=1 + one locally resolved mode: 4/4
Leading asymptotic sigma within 5.0% among resolved modes: 3/4
epsilon=0.050: Nq=768, rank=1, rank_stable=yes, strict=1, seed=strict-beyn, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=1.090%
epsilon=0.070: Nq=384, rank=1, rank_stable=yes, strict=1, seed=strict-beyn, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=2.436%
epsilon=0.090: Nq=384, rank=1, rank_stable=yes, strict=1, seed=strict-beyn, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=4.375%
epsilon=0.110: Nq=384, rank=1, rank_stable=yes, strict=1, seed=strict-beyn, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=6.876%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/theorem_2.3/theorem_2_1_beyn_v3_adaptive_validation
(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/theorem_2.3 %

---
