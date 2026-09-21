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

(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations % python first_theorem_beyn.py
=== Theorem 2.1 numerical validation: Beyn + BEM + SVD (v4) ===
b = 1.0
a = 0.6
leading-order a0* = 0.391826552031
geometric leading condition a > a0*: PASS
Lambda_1 = 2.467401100272
sqrt(Lambda_1)b = 1.570796326795
epsilons = (0.01, 0.03111111, 0.05222222, 0.07333333, 0.09444444, 0.11555556, 0.13666667, 0.15777778, 0.17888889, 0.2)
Beyn global discovery uses fixed M=24 (no multi-M global Beyn sweep, by design).
Beyn Nq base=(96, 192, 384), adaptive=(768, 1536), probe_dim=8
rank diagnostics: empty-S0=1.0e-09, rel_tol=1.0e-08, gap_threshold=1.0e+03
local BEM refinement M = (16, 24, 32, 40, 48)

Checking complex-kb support required by Beyn...
complex-kb assembly: PASS

=== epsilon=0.01000000 ===
predicted kb = 1.570796282404
predicted sigma = 3.73441725e-04
predicted cutoff gap = 4.439e-08
Beyn ellipse: real=[0.000100000000, 1.570796315697], center=0.785448157849, rx=0.785348157849, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.929e+05
Beyn estimated enclosed rank = 1
singular values of S0: [3.974e-05, 2.060e-10, 4.511e-13, 4.116e-13, 3.397e-13, 3.108e-13, 2.583e-13, 2.302e-13]
max contour linear-solve residual = 3.769e-16
raw Beyn eigenvalues:
[0] +1.571082768559 -3.972e-10i outside-contour/band
Nq=96: rank=1, strict=0, near-contour=0, S0-change=--, candidate-shift=--, raw=1.571082768559-3.97e-10i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.650e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.690e-05, 5.784e-11, 4.131e-13, 3.328e-13, 2.990e-13, 2.579e-13, 1.834e-13, 1.569e-13]
max contour linear-solve residual = 3.935e-16
raw Beyn eigenvalues:
[0] +1.570915147056 -2.295e-10i outside-contour/band
Nq=192: rank=1, strict=0, near-contour=0, S0-change=4.776e-01, candidate-shift=--, raw=1.570915147056-2.30e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.038e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.867e-05, 1.799e-11, 3.745e-13, 3.082e-13, 2.688e-13, 2.546e-13, 2.071e-13, 1.841e-13]
max contour linear-solve residual = 3.986e-16
raw Beyn eigenvalues:
[0] +1.570849499731 -5.941e-10i outside-contour/band
Nq=384: rank=1, strict=0, near-contour=0, S0-change=4.404e-01, candidate-shift=--, raw=1.570849499731-5.94e-10i
continuing adaptive quadrature: no strict candidate
adaptive Beyn escalation -> Nq=768
Beyn contour discovery: M=24, quadrature=768, probe_dim=8
contour node 128/768
contour node 256/768
contour node 384/768
contour node 512/768
contour node 640/768
contour node 768/768
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.213e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.320e-05, 5.963e-12, 3.369e-13, 2.869e-13, 2.374e-13, 2.041e-13, 1.961e-13, 1.593e-13]
max contour linear-solve residual = 4.191e-16
raw Beyn eigenvalues:
[0] +1.570821137983 +2.022e-10i outside-contour/band
Nq=768: rank=1, strict=0, near-contour=0, S0-change=4.149e-01, candidate-shift=--, raw=1.570821137983+2.02e-10i
continuing adaptive quadrature: no strict candidate
adaptive Beyn escalation -> Nq=1536
Beyn contour discovery: M=24, quadrature=1536, probe_dim=8
contour node 256/1536
contour node 512/1536
contour node 768/1536
contour node 1024/1536
contour node 1280/1536
contour node 1536/1536
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.649e+06
Beyn estimated enclosed rank = 1
singular values of S0: [9.453e-06, 2.033e-12, 2.235e-13, 1.860e-13, 1.560e-13, 1.309e-13, 1.113e-13, 8.372e-14]
max contour linear-solve residual = 4.232e-16
raw Beyn eigenvalues:
[0] +1.570808129446 -1.301e-09i outside-contour/band
Nq=1536: rank=1, strict=0, near-contour=0, S0-change=3.961e-01, candidate-shift=--, raw=1.570808129446-1.30e-09i
effective cutoff margin = 1.110e-08
final strict near-real Beyn candidates = 0
local refinement seeds = 0
tighter-cutoff diagnostic: margin 1.110e-08 -> 2.220e-09
Beyn contour discovery: M=24, quadrature=1536, probe_dim=8
contour node 256/1536
contour node 512/1536
contour node 768/1536
contour node 1024/1536
contour node 1280/1536
contour node 1536/1536
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.655e+06
Beyn estimated enclosed rank = 1
singular values of S0: [9.456e-06, 2.032e-12, 2.235e-13, 1.891e-13, 1.520e-13, 1.416e-13, 1.087e-13, 9.110e-14]
max contour linear-solve residual = 4.180e-16
raw Beyn eigenvalues:
[0] +1.570808132284 +8.404e-10i outside-contour/band
tighter-cutoff rank consistency = PASS (1 -> 1)

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.110e-08
final Beyn quadrature Nq = 1536
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 0
local refinement seeds = 0
local seed source = none
Aitken diagnostic kb = 1.570797107514
locally resolved modes = 0
one-mode count supported = FAIL
kb asymptotic = 1.570796282404
kb BEM = --
sigma asym = 3.73441725e-04
sigma BEM = --
relative singular value = --
relative sigma error = --
asymptotic accuracy <= 5.0%: N/A (requires exactly one certified BEM mode)

=== epsilon=0.03111111 ===
predicted kb = 1.570792168087
predicted sigma = 3.61454681e-03
predicted cutoff gap = 4.159e-06
Beyn ellipse: real=[0.000100000000, 1.570795287118], center=0.785447643559, rx=0.785347643559, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=3, gap-rank=1, selected gap=7.150e+04
Beyn estimated enclosed rank = 1
singular values of S0: [1.403e-04, 1.963e-09, 1.626e-12, 1.349e-12, 1.104e-12, 1.025e-12, 7.921e-13, 6.599e-13]
max contour linear-solve residual = 4.091e-16
raw Beyn eigenvalues:
[0] +1.571041826753 -1.099e-09i outside-contour/band
Nq=96: rank=1, strict=0, near-contour=0, S0-change=--, candidate-shift=--, raw=1.571041826753-1.10e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=3, gap-rank=1, selected gap=1.847e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.015e-04, 5.497e-10, 1.081e-12, 9.541e-13, 8.382e-13, 7.573e-13, 6.806e-13, 4.714e-13]
max contour linear-solve residual = 3.969e-16
raw Beyn eigenvalues:
[0] +1.570888853351 -7.383e-10i outside-contour/band
Nq=192: rank=1, strict=0, near-contour=0, S0-change=3.824e-01, candidate-shift=--, raw=1.570888853351-7.38e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=4, gap-rank=1, selected gap=4.542e+05
Beyn estimated enclosed rank = 1
singular values of S0: [7.729e-05, 1.702e-10, 1.209e-12, 1.052e-12, 9.486e-13, 8.194e-13, 7.386e-13, 5.711e-13]
max contour linear-solve residual = 4.091e-16
raw Beyn eigenvalues:
[0] +1.570831453789 +8.098e-10i outside-contour/band
Nq=384: rank=1, strict=0, near-contour=0, S0-change=3.134e-01, candidate-shift=--, raw=1.570831453789+8.10e-10i
continuing adaptive quadrature: no strict candidate
adaptive Beyn escalation -> Nq=768
Beyn contour discovery: M=24, quadrature=768, probe_dim=8
contour node 128/768
contour node 256/768
contour node 384/768
contour node 512/768
contour node 640/768
contour node 768/768
rank diagnostics: threshold-rank=3, gap-rank=1, selected gap=1.118e+06
Beyn estimated enclosed rank = 1
singular values of S0: [6.211e-05, 5.555e-11, 1.044e-12, 9.664e-13, 7.975e-13, 6.860e-13, 6.118e-13, 4.818e-13]
max contour linear-solve residual = 5.131e-16
raw Beyn eigenvalues:
[0] +1.570808154879 -2.173e-11i outside-contour/band
Nq=768: rank=1, strict=0, near-contour=0, S0-change=2.443e-01, candidate-shift=--, raw=1.570808154879-2.17e-11i
continuing adaptive quadrature: no strict candidate
adaptive Beyn escalation -> Nq=1536
Beyn contour discovery: M=24, quadrature=1536, probe_dim=8
contour node 256/1536
contour node 512/1536
contour node 768/1536
contour node 1024/1536
contour node 1280/1536
contour node 1536/1536
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.875e+06
Beyn estimated enclosed rank = 1
singular values of S0: [5.324e-05, 1.852e-11, 7.344e-13, 6.178e-13, 5.148e-13, 4.222e-13, 3.904e-13, 3.346e-13]
max contour linear-solve residual = 7.480e-16
raw Beyn eigenvalues:
[0] +1.570798418157 -2.737e-10i outside-contour/band
Nq=1536: rank=1, strict=0, near-contour=0, S0-change=1.667e-01, candidate-shift=--, raw=1.570798418157-2.74e-10i
no strict final candidate; using Aitken Delta^2 only as local-SVD seed: kb=1.570791427858
effective cutoff margin = 1.040e-06
final strict near-real Beyn candidates = 0
local refinement seeds = 1
candidate 1 [aitken]: seed kb=1.570791427858 +0.000e+00i, local bracket=[1.568791427858, 1.570795287118]
M=16: kb=1.570792195434, sigma_BEM=3.60264275e-03, sv_min=4.137e-06, sv_min/sv_max=1.440e-07, drop=8.45e+03, mesh_change=--, interior=no
M=24: kb=1.570792195436, sigma_BEM=3.60264204e-03, sv_min=4.111e-06, sv_min/sv_max=1.430e-07, drop=8.51e+03, mesh_change=0.000%, interior=no
M=32: kb=1.570792195436, sigma_BEM=3.60264188e-03, sv_min=4.103e-06, sv_min/sv_max=1.428e-07, drop=8.52e+03, mesh_change=0.000%, interior=no
M=40: kb=1.570792195437, sigma_BEM=3.60264182e-03, sv_min=4.099e-06, sv_min/sv_max=1.426e-07, drop=8.53e+03, mesh_change=0.000%, interior=no
M=48: kb=1.570792195437, sigma_BEM=3.60264179e-03, sv_min=4.099e-06, sv_min/sv_max=1.426e-07, drop=8.53e+03, mesh_change=0.000%, interior=no
resolved=no
tighter-cutoff diagnostic: margin 1.040e-06 -> 2.079e-07
Beyn contour discovery: M=24, quadrature=1536, probe_dim=8
contour node 256/1536
contour node 512/1536
contour node 768/1536
contour node 1024/1536
contour node 1280/1536
contour node 1536/1536
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.859e+06
Beyn estimated enclosed rank = 1
singular values of S0: [5.442e-05, 1.904e-11, 7.509e-13, 6.029e-13, 5.204e-13, 4.749e-13, 3.702e-13, 3.329e-13]
max contour linear-solve residual = 7.378e-16
raw Beyn eigenvalues:
[0] +1.570798486469 -1.195e-09i outside-contour/band
tighter-cutoff rank consistency = PASS (1 -> 1)

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.040e-06
final Beyn quadrature Nq = 1536
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 0
local refinement seeds = 1
local seed source = aitken
Aitken diagnostic kb = 1.570791427858
locally resolved modes = 0
one-mode count supported = FAIL
kb asymptotic = 1.570792168087
kb BEM = --
sigma asym = 3.61454681e-03
sigma BEM = --
relative singular value = --
relative sigma error = --
asymptotic accuracy <= 5.0%: N/A (requires exactly one certified BEM mode)

=== epsilon=0.05222222 ===
predicted kb = 1.570763311005
predicted sigma = 1.01843543e-02
predicted cutoff gap = 3.302e-05
Beyn ellipse: real=[0.000100000000, 1.570788072847], center=0.785444036424, rx=0.785344036424, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.685e+04
Beyn estimated enclosed rank = 1
singular values of S0: [3.069e-04, 5.397e-09, 2.393e-12, 2.110e-12, 1.822e-12, 1.657e-12, 1.297e-12, 9.638e-13]
max contour linear-solve residual = 4.209e-16
raw Beyn eigenvalues:
[0] +1.570951022631 -1.992e-09i outside-contour/band
Nq=96: rank=1, strict=0, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570951022631-1.99e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.709e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.538e-04, 1.486e-09, 2.024e-12, 1.868e-12, 1.766e-12, 1.488e-12, 1.192e-12, 9.181e-13]
max contour linear-solve residual = 3.979e-16
raw Beyn eigenvalues:
[0] +1.570826316386 +3.048e-10i outside-contour/band
Nq=192: rank=1, strict=0, near-contour=0, S0-change=2.091e-01, candidate-shift=--, raw=1.570826316386+3.05e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.231e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.318e-04, 4.431e-10, 2.269e-12, 1.979e-12, 1.842e-12, 1.580e-12, 1.282e-12, 1.047e-12]
max contour linear-solve residual = 6.429e-16
raw Beyn eigenvalues:
[0] +1.570784419936 +9.014e-12i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=9.516e-02, candidate-shift=--, raw=1.570784419936+9.01e-12i
effective cutoff margin = 8.254e-06
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570784419936 +9.014e-12i, local bracket=[1.568784419936, 1.570788072847]
M=16: kb=1.570764092065, sigma_BEM=1.00631681e-02, sv_min=1.340e-05, sv_min/sv_max=7.720e-07, drop=4.35e+03, mesh_change=--, interior=yes
M=24: kb=1.570764092267, sigma_BEM=1.00631365e-02, sv_min=1.322e-05, sv_min/sv_max=7.621e-07, drop=4.40e+03, mesh_change=0.000%, interior=yes
M=32: kb=1.570764092315, sigma_BEM=1.00631290e-02, sv_min=1.319e-05, sv_min/sv_max=7.598e-07, drop=4.42e+03, mesh_change=0.000%, interior=yes
M=40: kb=1.570764092333, sigma_BEM=1.00631262e-02, sv_min=1.317e-05, sv_min/sv_max=7.589e-07, drop=4.42e+03, mesh_change=0.000%, interior=yes
M=48: kb=1.570764092341, sigma_BEM=1.00631250e-02, sv_min=1.316e-05, sv_min/sv_max=7.585e-07, drop=4.42e+03, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 8.254e-06 -> 1.651e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.123e+05
Beyn estimated enclosed rank = 1
singular values of S0: [2.410e-04, 4.704e-10, 2.143e-12, 1.878e-12, 1.519e-12, 1.174e-12, 1.011e-12, 7.582e-13]
max contour linear-solve residual = 6.160e-16
raw Beyn eigenvalues:
[0] +1.570784856218 -3.158e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=1.900e-10, Neumann residual=6.290e+00
decay rates left/right=1.006021e-02/1.006020e-02, expected sigma=1.006312e-02
NOTE: final strict Beyn candidate position is still moving by more than 2.0e-04 between the last two Nq levels. Local SVD refinement remains the trusted position.

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 8.254e-06
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570763223018
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570763311005
kb BEM = 1.570764092341
sigma asym = 1.01843543e-02
sigma BEM = 1.00631250e-02
sigma/epsilon^2 = 3.68996466e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 2.88329444e-01
final sigma_min(A) = 1.316e-05
final sigma_max(A) = 1.735e+01
final sigma_min/sigma_max = 7.585e-07
final minimum drop = 4.42e+03
final mesh change in sigma = 0.000%
relative sigma error = 1.190%
wall Dirichlet residual = 1.900e-10
obstacle Neumann residual = 6.290e+00
decay rate left/right = 1.006021e-02/1.006020e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.07333333 ===
predicted kb = 1.570667940347
predicted sigma = 2.00828643e-02
predicted cutoff gap = 1.284e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=6.335e+04
Beyn estimated enclosed rank = 1
singular values of S0: [6.650e-04, 1.050e-08, 3.776e-12, 2.980e-12, 2.567e-12, 2.037e-12, 1.732e-12, 1.614e-12]
max contour linear-solve residual = 4.157e-16
raw Beyn eigenvalues:
[0] +1.570793717492 -1.595e-09i near-contour-seed
Nq=96: rank=1, strict=0, near-contour=1, S0-change=--, candidate-shift=--, raw=1.570793717492-1.59e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=2.350e+05
Beyn estimated enclosed rank = 1
singular values of S0: [6.760e-04, 2.877e-09, 3.138e-12, 2.589e-12, 2.460e-12, 1.830e-12, 1.574e-12, 1.140e-12]
max contour linear-solve residual = 6.221e-16
raw Beyn eigenvalues:
[0] +1.570706888499 -5.861e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=1.625e-02, candidate-shift=--, raw=1.570706888499-5.86e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=9.013e+05
Beyn estimated enclosed rank = 1
singular values of S0: [7.665e-04, 8.504e-10, 3.035e-12, 2.373e-12, 1.992e-12, 1.654e-12, 1.589e-12, 1.224e-12]
max contour linear-solve residual = 1.178e-15
raw Beyn eigenvalues:
[0] +1.570683181864 +2.769e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=1.180e-01, candidate-shift=2.371e-05, raw=1.570683181864+2.77e-10i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570683181864 +2.769e-10i, local bracket=[1.568683181864, 1.570786326795]
M=16: kb=1.570674808348, sigma_BEM=1.95383391e-02, sv_min=8.757e-06, sv_min/sv_max=6.929e-07, drop=1.35e+04, mesh_change=--, interior=yes
M=24: kb=1.570674809063, sigma_BEM=1.95382816e-02, sv_min=8.811e-06, sv_min/sv_max=6.972e-07, drop=1.34e+04, mesh_change=0.000%, interior=yes
M=32: kb=1.570674809237, sigma_BEM=1.95382676e-02, sv_min=8.822e-06, sv_min/sv_max=6.981e-07, drop=1.34e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570674809303, sigma_BEM=1.95382623e-02, sv_min=8.824e-06, sv_min/sv_max=6.982e-07, drop=1.34e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570674809326, sigma_BEM=1.95382604e-02, sv_min=8.828e-06, sv_min/sv_max=6.986e-07, drop=1.34e+04, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.560e+05
Beyn estimated enclosed rank = 1
singular values of S0: [7.833e-04, 9.151e-10, 3.026e-12, 2.422e-12, 2.019e-12, 1.756e-12, 1.563e-12, 1.202e-12]
max contour linear-solve residual = 9.941e-16
raw Beyn eigenvalues:
[0] +1.570683623958 -2.336e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=1.353e-10, Neumann residual=6.287e+00
decay rates left/right=1.953176e-02/1.953176e-02, expected sigma=1.953826e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570674278450
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570667940347
kb BEM = 1.570674809326
sigma asym = 2.00828643e-02
sigma BEM = 1.95382604e-02
sigma/epsilon^2 = 3.63314793e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 5.28542931e-01
final sigma_min(A) = 8.828e-06
final sigma_max(A) = 1.264e+01
final sigma_min/sigma_max = 6.986e-07
final minimum drop = 1.34e+04
final mesh change in sigma = 0.000%
relative sigma error = 2.712%
wall Dirichlet residual = 1.353e-10
obstacle Neumann residual = 6.287e+00
decay rate left/right = 1.953176e-02/1.953176e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.09444444 ===
predicted kb = 1.570443102779
predicted sigma = 3.33100766e-02
predicted cutoff gap = 3.532e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=8.846e+04
Beyn estimated enclosed rank = 1
singular values of S0: [1.525e-03, 1.724e-08, 4.751e-12, 4.179e-12, 3.769e-12, 2.604e-12, 2.399e-12, 2.050e-12]
max contour linear-solve residual = 5.244e-16
raw Beyn eigenvalues:
[0] +1.570542281546 -2.198e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570542281546-2.20e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=3.681e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.739e-03, 4.725e-09, 3.702e-12, 3.252e-12, 2.958e-12, 2.451e-12, 2.117e-12, 1.804e-12]
max contour linear-solve residual = 9.234e-16
raw Beyn eigenvalues:
[0] +1.570492504463 -4.804e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=1.231e-01, candidate-shift=4.978e-05, raw=1.570492504463-4.80e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.327e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.854e-03, 1.397e-09, 3.932e-12, 3.344e-12, 3.122e-12, 2.342e-12, 2.113e-12, 1.484e-12]
max contour linear-solve residual = 1.336e-15
raw Beyn eigenvalues:
[0] +1.570481120514 -7.526e-11i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=6.180e-02, candidate-shift=1.138e-05, raw=1.570481120514-7.53e-11i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570481120514 -7.526e-11i, local bracket=[1.568481120514, 1.570786326795]
M=16: kb=1.570476733888, sigma_BEM=3.16848322e-02, sv_min=4.711e-06, sv_min/sv_max=4.654e-07, drop=3.73e+04, mesh_change=--, interior=yes
M=24: kb=1.570476738028, sigma_BEM=3.16846269e-02, sv_min=4.476e-06, sv_min/sv_max=4.422e-07, drop=3.92e+04, mesh_change=0.001%, interior=yes
M=32: kb=1.570476739019, sigma_BEM=3.16845778e-02, sv_min=4.420e-06, sv_min/sv_max=4.367e-07, drop=3.97e+04, mesh_change=0.000%, interior=yes
M=40: kb=1.570476739369, sigma_BEM=3.16845605e-02, sv_min=4.401e-06, sv_min/sv_max=4.347e-07, drop=3.99e+04, mesh_change=0.000%, interior=yes
M=48: kb=1.570476739523, sigma_BEM=3.16845529e-02, sv_min=4.392e-06, sv_min/sv_max=4.339e-07, drop=4.00e+04, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.236e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.857e-03, 1.503e-09, 3.665e-12, 3.423e-12, 2.817e-12, 2.568e-12, 2.162e-12, 1.831e-12]
max contour linear-solve residual = 1.309e-15
raw Beyn eigenvalues:
[0] +1.570481443313 -1.121e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=1.025e-10, Neumann residual=6.285e+00
decay rates left/right=3.167454e-02/3.167454e-02, expected sigma=3.168455e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570477745059
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570443102779
kb BEM = 1.570476739523
sigma asym = 3.33100766e-02
sigma BEM = 3.16845529e-02
sigma/epsilon^2 = 3.55217858e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 8.17710040e-01
final sigma_min(A) = 4.392e-06
final sigma_max(A) = 1.012e+01
final sigma_min/sigma_max = 4.339e-07
final minimum drop = 4.00e+04
final mesh change in sigma = 0.000%
relative sigma error = 4.880%
wall Dirichlet residual = 1.025e-10
obstacle Neumann residual = 6.285e+00
decay rate left/right = 3.167454e-02/3.167454e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: PASS

=== epsilon=0.11555556 ===
predicted kb = 1.570004612194
predicted sigma = 4.98660001e-02
predicted cutoff gap = 7.917e-04
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.286e+05
Beyn estimated enclosed rank = 1
singular values of S0: [3.290e-03, 2.558e-08, 5.879e-12, 5.085e-12, 4.462e-12, 3.866e-12, 3.194e-12, 2.406e-12]
max contour linear-solve residual = 7.067e-16
raw Beyn eigenvalues:
[0] +1.570157721550 -2.132e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.570157721550-2.13e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=4.647e+05
Beyn estimated enclosed rank = 1
singular values of S0: [3.258e-03, 7.011e-09, 4.700e-12, 3.644e-12, 3.549e-12, 2.771e-12, 2.433e-12, 1.972e-12]
max contour linear-solve residual = 9.592e-16
raw Beyn eigenvalues:
[0] +1.570131436145 -7.865e-11i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=9.627e-03, candidate-shift=2.629e-05, raw=1.570131436145-7.86e-11i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.560e+06
Beyn estimated enclosed rank = 1
singular values of S0: [3.232e-03, 2.071e-09, 4.455e-12, 3.154e-12, 2.778e-12, 2.554e-12, 2.223e-12, 1.781e-12]
max contour linear-solve residual = 1.230e-15
raw Beyn eigenvalues:
[0] +1.570124373581 +2.340e-11i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=8.036e-03, candidate-shift=7.063e-06, raw=1.570124373581+2.34e-11i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.570124373581 +2.340e-11i, local bracket=[1.568124373581, 1.570786326795]
M=16: kb=1.570121353802, sigma_BEM=4.60438335e-02, sv_min=2.170e-06, sv_min/sv_max=2.523e-07, drop=1.05e+05, mesh_change=--, interior=yes
M=24: kb=1.570121371033, sigma_BEM=4.60432459e-02, sv_min=9.349e-07, sv_min/sv_max=1.087e-07, drop=2.43e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.570121376211, sigma_BEM=4.60430694e-02, sv_min=4.422e-07, sv_min/sv_max=5.140e-08, drop=5.15e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.570121378104, sigma_BEM=4.60430048e-02, sv_min=2.563e-07, sv_min/sv_max=2.979e-08, drop=8.88e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.570121378949, sigma_BEM=4.60429760e-02, sv_min=1.727e-07, sv_min/sv_max=2.008e-08, drop=1.32e+06, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.451e+06
Beyn estimated enclosed rank = 1
singular values of S0: [3.233e-03, 2.229e-09, 3.959e-12, 3.382e-12, 2.944e-12, 2.763e-12, 2.562e-12, 1.969e-12]
max contour linear-solve residual = 1.296e-15
raw Beyn eigenvalues:
[0] +1.570124600270 -1.207e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=1.391e-10, Neumann residual=6.283e+00
decay rates left/right=4.602595e-02/4.602595e-02, expected sigma=4.604298e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.570121778761
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.570004612194
kb BEM = 1.570121378949
sigma asym = 4.98660001e-02
sigma BEM = 4.60429760e-02
sigma/epsilon^2 = 3.44811462e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.14810678e+00
final sigma_min(A) = 1.727e-07
final sigma_max(A) = 8.603e+00
final sigma_min/sigma_max = 2.008e-08
final minimum drop = 1.32e+06
final mesh change in sigma = 0.000%
relative sigma error = 7.667%
wall Dirichlet residual = 1.391e-10
obstacle Neumann residual = 6.283e+00
decay rate left/right = 4.602595e-02/4.602595e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.13666667 ===
predicted kb = 1.569246937686
predicted sigma = 6.97506189e-02
predicted cutoff gap = 1.549e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.481e+05
Beyn estimated enclosed rank = 1
singular values of S0: [5.262e-03, 3.552e-08, 5.815e-12, 5.175e-12, 4.570e-12, 4.262e-12, 3.453e-12, 2.405e-12]
max contour linear-solve residual = 8.229e-16
raw Beyn eigenvalues:
[0] +1.569595648850 -2.647e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.569595648850-2.65e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.163e+05
Beyn estimated enclosed rank = 1
singular values of S0: [5.027e-03, 9.736e-09, 5.322e-12, 4.944e-12, 4.395e-12, 3.503e-12, 2.945e-12, 2.359e-12]
max contour linear-solve residual = 1.044e-15
raw Beyn eigenvalues:
[0] +1.569576980368 -4.309e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=4.667e-02, candidate-shift=1.867e-05, raw=1.569576980368-4.31e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.742e+06
Beyn estimated enclosed rank = 1
singular values of S0: [5.011e-03, 2.877e-09, 5.509e-12, 4.485e-12, 4.205e-12, 3.462e-12, 2.678e-12, 2.196e-12]
max contour linear-solve residual = 1.241e-15
raw Beyn eigenvalues:
[0] +1.569571695662 -4.834e-11i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=3.158e-03, candidate-shift=5.285e-06, raw=1.569571695662-4.83e-11i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.569571695662 -4.834e-11i, local bracket=[1.567571695662, 1.570786326795]
M=16: kb=1.569569435517, sigma_BEM=6.20716309e-02, sv_min=6.993e-07, sv_min/sv_max=9.177e-08, drop=2.43e+05, mesh_change=--, interior=yes
M=24: kb=1.569569475286, sigma_BEM=6.20706253e-02, sv_min=9.040e-07, sv_min/sv_max=1.186e-07, drop=1.88e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.569569479191, sigma_BEM=6.20705265e-02, sv_min=6.181e-07, sv_min/sv_max=8.112e-08, drop=2.75e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.569569480571, sigma_BEM=6.20704916e-02, sv_min=5.169e-07, sv_min/sv_max=6.784e-08, drop=3.29e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.569569481179, sigma_BEM=6.20704763e-02, sv_min=4.723e-07, sv_min/sv_max=6.198e-08, drop=3.61e+05, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.619e+06
Beyn estimated enclosed rank = 1
singular values of S0: [5.012e-03, 3.095e-09, 5.859e-12, 4.405e-12, 3.758e-12, 2.913e-12, 2.356e-12, 2.211e-12]
max contour linear-solve residual = 1.187e-15
raw Beyn eigenvalues:
[0] +1.569571864782 +1.081e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=1.027e-10, Neumann residual=6.281e+00
decay rates left/right=6.204497e-02/6.204497e-02, expected sigma=6.207048e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.569569608948
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.569246937686
kb BEM = 1.569569481179
sigma asym = 6.97506189e-02
sigma BEM = 6.20704763e-02
sigma/epsilon^2 = 3.32322581e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.51175864e+00
final sigma_min(A) = 4.723e-07
final sigma_max(A) = 7.620e+00
final sigma_min/sigma_max = 6.198e-08
final minimum drop = 3.61e+05
final mesh change in sigma = 0.000%
relative sigma error = 11.011%
wall Dirichlet residual = 1.027e-10
obstacle Neumann residual = 6.281e+00
decay rate left/right = 6.204497e-02/6.204497e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.15777778 ===
predicted kb = 1.568042986053
predicted sigma = 9.29639401e-02
predicted cutoff gap = 2.753e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.464e+05
Beyn estimated enclosed rank = 1
singular values of S0: [6.895e-03, 4.711e-08, 7.473e-12, 6.072e-12, 4.898e-12, 4.582e-12, 3.491e-12, 2.298e-12]
max contour linear-solve residual = 7.876e-16
raw Beyn eigenvalues:
[0] +1.568822008946 -3.623e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.568822008946-3.62e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.548e+05
Beyn estimated enclosed rank = 1
singular values of S0: [7.165e-03, 1.291e-08, 6.097e-12, 5.249e-12, 4.767e-12, 4.163e-12, 2.928e-12, 2.740e-12]
max contour linear-solve residual = 9.625e-16
raw Beyn eigenvalues:
[0] +1.568805533087 -7.297e-10i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=3.760e-02, candidate-shift=1.648e-05, raw=1.568805533087-7.30e-10i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.875e+06
Beyn estimated enclosed rank = 1
singular values of S0: [7.157e-03, 3.818e-09, 5.957e-12, 5.376e-12, 5.022e-12, 4.240e-12, 3.649e-12, 2.875e-12]
max contour linear-solve residual = 1.086e-15
raw Beyn eigenvalues:
[0] +1.568801375483 -2.048e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=1.068e-03, candidate-shift=4.158e-06, raw=1.568801375483-2.05e-10i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.568801375483 -2.048e-10i, local bracket=[1.566801375483, 1.570786326795]
M=16: kb=1.568799576685, sigma_BEM=7.91769440e-02, sv_min=3.306e-07, sv_min/sv_max=4.750e-08, drop=3.95e+05, mesh_change=--, interior=yes
M=24: kb=1.568799619952, sigma_BEM=7.91760867e-02, sv_min=6.916e-07, sv_min/sv_max=9.938e-08, drop=1.89e+05, mesh_change=0.001%, interior=yes
M=32: kb=1.568799632328, sigma_BEM=7.91758415e-02, sv_min=7.674e-07, sv_min/sv_max=1.103e-07, drop=1.70e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.568799636922, sigma_BEM=7.91757505e-02, sv_min=7.761e-07, sv_min/sv_max=1.115e-07, drop=1.68e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.568799638981, sigma_BEM=7.91757097e-02, sv_min=7.771e-07, sv_min/sv_max=1.117e-07, drop=1.68e+05, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.742e+06
Beyn estimated enclosed rank = 1
singular values of S0: [7.157e-03, 4.108e-09, 5.955e-12, 5.172e-12, 4.376e-12, 4.249e-12, 3.116e-12, 2.586e-12]
max contour linear-solve residual = 1.192e-15
raw Beyn eigenvalues:
[0] +1.568801508127 -2.800e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=7.236e-11, Neumann residual=6.280e+00
decay rates left/right=7.914006e-02/7.914006e-02, expected sigma=7.917571e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.568799972227
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.568042986053
kb BEM = 1.568799638981
sigma asym = 9.29639401e-02
sigma BEM = 7.91757097e-02
sigma/epsilon^2 = 3.18053576e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 1.90110357e+00
final sigma_min(A) = 7.771e-07
final sigma_max(A) = 6.959e+00
final sigma_min/sigma_max = 1.117e-07
final minimum drop = 1.68e+05
final mesh change in sigma = 0.000%
relative sigma error = 14.832%
wall Dirichlet residual = 7.236e-11
obstacle Neumann residual = 6.280e+00
decay rate left/right = 7.914006e-02/7.914006e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.17888889 ===
predicted kb = 1.566243730998
predicted sigma = 1.19505964e-01
predicted cutoff gap = 4.553e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.562e+05
Beyn estimated enclosed rank = 1
singular values of S0: [9.451e-03, 6.049e-08, 8.592e-12, 6.863e-12, 5.989e-12, 5.781e-12, 4.545e-12, 3.821e-12]
max contour linear-solve residual = 6.808e-16
raw Beyn eigenvalues:
[0] +1.567831307108 -4.531e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.567831307108-4.53e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.793e+05
Beyn estimated enclosed rank = 1
singular values of S0: [9.604e-03, 1.658e-08, 6.497e-12, 6.166e-12, 4.962e-12, 4.600e-12, 3.289e-12, 2.743e-12]
max contour linear-solve residual = 8.682e-16
raw Beyn eigenvalues:
[0] +1.567818239959 -1.102e-09i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=1.597e-02, candidate-shift=1.307e-05, raw=1.567818239959-1.10e-09i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.956e+06
Beyn estimated enclosed rank = 1
singular values of S0: [9.585e-03, 4.900e-09, 7.551e-12, 6.699e-12, 5.298e-12, 5.037e-12, 3.841e-12, 2.832e-12]
max contour linear-solve residual = 1.091e-15
raw Beyn eigenvalues:
[0] +1.567814842485 -4.152e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=2.003e-03, candidate-shift=3.397e-06, raw=1.567814842485-4.15e-10i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.567814842485 -4.152e-10i, local bracket=[1.565814842485, 1.569814842485]
M=16: kb=1.567813315236, sigma_BEM=9.67600581e-02, sv_min=3.675e-07, sv_min/sv_max=5.648e-08, drop=2.81e+05, mesh_change=--, interior=yes
M=24: kb=1.567813425969, sigma_BEM=9.67582639e-02, sv_min=8.160e-07, sv_min/sv_max=1.254e-07, drop=1.27e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.567813452477, sigma_BEM=9.67578343e-02, sv_min=9.257e-07, sv_min/sv_max=1.423e-07, drop=1.12e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.567813461858, sigma_BEM=9.67576823e-02, sv_min=9.648e-07, sv_min/sv_max=1.483e-07, drop=1.07e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.567813465992, sigma_BEM=9.67576153e-02, sv_min=9.821e-07, sv_min/sv_max=1.509e-07, drop=1.05e+05, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.818e+06
Beyn estimated enclosed rank = 1
singular values of S0: [9.585e-03, 5.272e-09, 6.631e-12, 5.950e-12, 4.929e-12, 3.975e-12, 3.727e-12, 3.004e-12]
max contour linear-solve residual = 1.148e-15
raw Beyn eigenvalues:
[0] +1.567814951137 -4.127e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=2.152e-10, Neumann residual=6.280e+00
decay rates left/right=9.670573e-02/9.670573e-02, expected sigma=9.675762e-02

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.567813648770
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.566243730998
kb BEM = 1.567813465992
sigma asym = 1.19505964e-01
sigma BEM = 9.67576153e-02
sigma/epsilon^2 = 3.02355879e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 2.30898610e+00
final sigma_min(A) = 9.821e-07
final sigma_max(A) = 6.507e+00
final sigma_min/sigma_max = 1.509e-07
final minimum drop = 1.05e+05
final mesh change in sigma = 0.000%
relative sigma error = 19.035%
wall Dirichlet residual = 2.152e-10
obstacle Neumann residual = 6.280e+00
decay rate left/right = 9.670573e-02/9.670573e-02
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== epsilon=0.20000000 ===
predicted kb = 1.563677621758
predicted sigma = 1.49376690e-01
predicted cutoff gap = 7.119e-03
Beyn ellipse: real=[0.000100000000, 1.570786326795], center=0.785443163397, rx=0.785343163397, ry=3.000e-02
Beyn quadrature-convergence study:
Beyn contour discovery: M=24, quadrature=96, probe_dim=8
contour node 16/96
contour node 32/96
contour node 48/96
contour node 64/96
contour node 80/96
contour node 96/96
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.692e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.284e-02, 7.588e-08, 1.111e-11, 9.968e-12, 8.103e-12, 7.453e-12, 6.426e-12, 4.924e-12]
max contour linear-solve residual = 7.354e-16
raw Beyn eigenvalues:
[0] +1.566650425610 -5.564e-09i real-candidate
Nq=96: rank=1, strict=1, near-contour=0, S0-change=--, candidate-shift=--, raw=1.566650425610-5.56e-09i
Beyn contour discovery: M=24, quadrature=192, probe_dim=8
contour node 32/192
contour node 64/192
contour node 96/192
contour node 128/192
contour node 160/192
contour node 192/192
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=5.857e+05
Beyn estimated enclosed rank = 1
singular values of S0: [1.218e-02, 2.080e-08, 7.478e-12, 6.761e-12, 5.836e-12, 4.949e-12, 4.230e-12, 3.365e-12]
max contour linear-solve residual = 8.661e-16
raw Beyn eigenvalues:
[0] +1.566640338773 -1.359e-09i real-candidate
Nq=192: rank=1, strict=1, near-contour=0, S0-change=5.374e-02, candidate-shift=1.009e-05, raw=1.566640338773-1.36e-09i
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.982e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.219e-02, 6.149e-09, 7.118e-12, 5.994e-12, 5.475e-12, 4.835e-12, 4.153e-12, 3.721e-12]
max contour linear-solve residual = 1.079e-15
raw Beyn eigenvalues:
[0] +1.566637451734 -3.260e-10i real-candidate
Nq=384: rank=1, strict=1, near-contour=0, S0-change=5.051e-04, candidate-shift=2.887e-06, raw=1.566637451734-3.26e-10i
effective cutoff margin = 1.000e-05
final strict near-real Beyn candidates = 1
local refinement seeds = 1
candidate 1 [strict-beyn]: seed kb=1.566637451734 -3.260e-10i, local bracket=[1.564637451734, 1.568637451734]
M=16: kb=1.566636070811, sigma_BEM=1.14247634e-01, sv_min=1.875e-07, sv_min/sv_max=3.025e-08, drop=4.49e+05, mesh_change=--, interior=yes
M=24: kb=1.566636232516, sigma_BEM=1.14245417e-01, sv_min=4.203e-07, sv_min/sv_max=6.781e-08, drop=2.00e+05, mesh_change=0.002%, interior=yes
M=32: kb=1.566636272315, sigma_BEM=1.14244871e-01, sv_min=5.097e-07, sv_min/sv_max=8.223e-08, drop=1.65e+05, mesh_change=0.000%, interior=yes
M=40: kb=1.566636293893, sigma_BEM=1.14244575e-01, sv_min=1.791e-07, sv_min/sv_max=2.889e-08, drop=4.71e+05, mesh_change=0.000%, interior=yes
M=48: kb=1.566636300070, sigma_BEM=1.14244490e-01, sv_min=1.943e-07, sv_min/sv_max=3.134e-08, drop=4.34e+05, mesh_change=0.000%, interior=yes
resolved=YES
tighter-cutoff diagnostic: margin 1.000e-05 -> 2.000e-06
Beyn contour discovery: M=24, quadrature=384, probe_dim=8
contour node 64/384
contour node 128/384
contour node 192/384
contour node 256/384
contour node 320/384
contour node 384/384
rank diagnostics: threshold-rank=2, gap-rank=1, selected gap=1.843e+06
Beyn estimated enclosed rank = 1
singular values of S0: [1.219e-02, 6.615e-09, 7.366e-12, 5.799e-12, 5.463e-12, 4.601e-12, 3.999e-12, 3.224e-12]
max contour linear-solve residual = 1.071e-15
raw Beyn eigenvalues:
[0] +1.566637543623 -4.070e-10i real-candidate
tighter-cutoff rank consistency = PASS (1 -> 1)
reconstructing field for independent physical diagnostics...
wall residual=8.393e-11, Neumann residual=6.279e+00
decay rates left/right=1.141771e-01/1.141771e-01, expected sigma=1.142445e-01

--- Beyn-v4 + local-SVD validation result ---
effective cutoff margin = 1.000e-05
final Beyn quadrature Nq = 384
final Beyn estimated rank = 1
Beyn rank stable (last 2 Nq) = YES
tighter-cutoff rank check = PASS
final strict Beyn candidates = 1
local refinement seeds = 1
local seed source = strict-beyn
Aitken diagnostic kb = 1.566636294064
locally resolved modes = 1
one-mode count supported = PASS
kb asymptotic = 1.563677621758
kb BEM = 1.566636300070
sigma asym = 1.49376690e-01
sigma BEM = 1.14244490e-01
sigma/epsilon^2 = 2.85611226e+00
asymptotic coefficient C(a) = 3.73441725e+00
scaled asymptotic remainder = 2.72860786e+00
final sigma_min(A) = 1.943e-07
final sigma_max(A) = 6.199e+00
final sigma_min/sigma_max = 3.134e-08
final minimum drop = 4.34e+05
final mesh change in sigma = 0.000%
relative sigma error = 23.519%
wall Dirichlet residual = 8.393e-11
obstacle Neumann residual = 6.279e+00
decay rate left/right = 1.141771e-01/1.141771e-01
physical decay check = PASS
asymptotic accuracy <= 5.0%: FAIL

=== INTERNAL GREEN / DERIVATIVE CONVERGENCE at epsilon=0.094 ===
finite_difference_step=3e-06: kb=1.570476739075, sigma=3.16845751e-02, rel_sv=4.349e-07, delta_sigma=0.000%
finite_difference_step=1e-06: kb=1.570476739076, sigma=3.16845750e-02, rel_sv=4.348e-07, delta_sigma=0.000%
finite_difference_step=3e-07: kb=1.570476739075, sigma=3.16845750e-02, rel_sv=4.348e-07, delta_sigma=0.000%
lattice_terms=100: kb=1.570476739076, sigma=3.16845750e-02, rel_sv=4.349e-07, delta_sigma=0.000%
lattice_terms=200: kb=1.570476739076, sigma=3.16845750e-02, rel_sv=4.348e-07, delta_sigma=0.000%
lattice_terms=300: kb=1.570476739076, sigma=3.16845750e-02, rel_sv=4.349e-07, delta_sigma=0.000%
harmonic_order=12: kb=1.570476305063, sigma=3.17060800e-02, rel_sv=2.623e-07, delta_sigma=0.068%
harmonic_order=20: kb=1.570476739076, sigma=3.16845750e-02, rel_sv=4.348e-07, delta_sigma=0.000%
harmonic_order=28: kb=1.570476768163, sigma=3.16831332e-02, rel_sv=8.412e-08, delta_sigma=0.004%

=== FINAL SUMMARY ===
Global discovery: adaptive Beyn contour method with Nq convergence study
Global Beyn BEM order: one fixed M, as requested
Local certification: relative SVD singularity + BEM mesh refinement
Near-cutoff safeguard: epsilon-adaptive margin + tighter-margin rank check
Asymptotics: direct sigma/epsilon^2 and scaled O(epsilon^3 log epsilon) diagnostics
Stable rank=1 + exactly one locally resolved mode: 8/10
Leading asymptotic sigma within 5.0%: 3/8
Independent physical decay checks: 8/8
Critical-height a\*(epsilon) study: IMPLEMENTED but disabled. Set run_critical_height_study=True for the full existence/non-existence transition test.
epsilon=0.01000: margin=1.11e-08, Nq=1536, rank=1, rank_stable=yes, margin_check=PASS, resolved=0, one-mode=FAIL, asymptotic=N/A, sigma_error=N/A
epsilon=0.03111: margin=1.04e-06, Nq=1536, rank=1, rank_stable=yes, margin_check=PASS, resolved=0, one-mode=FAIL, asymptotic=N/A, sigma_error=N/A
epsilon=0.05222: margin=8.25e-06, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=1.190%
epsilon=0.07333: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=2.712%
epsilon=0.09444: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=PASS, sigma_error=4.880%
epsilon=0.11556: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=7.667%
epsilon=0.13667: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=11.011%
epsilon=0.15778: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=14.832%
epsilon=0.17889: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=19.035%
epsilon=0.20000: margin=1.00e-05, Nq=384, rank=1, rank_stable=yes, margin_check=PASS, resolved=1, one-mode=PASS, asymptotic=FAIL, sigma_error=23.519%

Files written to: /Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/the(base) kevinsepulveda@Kevins-Mac-mini:~/Documents/waveguide/src/beyn_theorems_validations %
