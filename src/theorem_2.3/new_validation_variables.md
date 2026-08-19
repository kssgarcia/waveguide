Sí. En los codigos de simetria agregamos varias cantidades que **no aparecen literalmente en el teorema**, sino que sirven como diagnosticos numericos para asegurarnos de que el candidato encontrado por BEM realmente corresponde al BIC que queremos y no a una singularidad espuria.

La idea general es esta:

[
\boxed{
\text{Teorema}
\rightarrow
\text{predice un BIC}
\rightarrow
\text{BEM encuentra un candidato}
\rightarrow
\text{diagnosticos confirman que realmente es un BIC}
}
]

Las variables nuevas mas importantes son estas.

### 1. `sv_min` o (s\_{\min}(A))

Esta es la misma idea que usamos en el Teorema 2.1.

Tenemos

[
A(k)\mathbf u=0.
]

Para que exista una solucion no trivial necesitamos que (A(k)) sea singular. En lugar de evaluar directamente

[
\det A(k)=0,
]

calculamos el menor valor singular:

[
\boxed{
s_{\min}(A(k))
}
]

mediante SVD.

Cuando

[
s_{\min}(A(k))\approx0,
]

la matriz esta cerca de ser singular y tenemos un candidato a valor propio.

Por ejemplo,

[
s_{\min}=10^{-6}
]

es una señal fuerte, mientras que

[
s_{\min}\approx0.97
]

indica que la matriz no esta ni remotamente cerca de ser singular.

Por eso los falsos candidatos que vimos cerca de (\Lambda_2), con

[
s_{\min}\approx0.97,
]

se descartan inmediatamente.

---

### 2. `drop_factor`

Esta variable mide **que tan profundo es el minimo de (s\_{\min}(A(k)))**.

Si el minimo esta en (k\_\*), tomamos valores a izquierda y derecha:

[
s_L=s_{\min}(A(k_L)),
]

[
s_R=s_{\min}(A(k_R)),
]

y en el minimo

[
s_*=s_{\min}(A(k_*)).
]

Entonces definimos algo del tipo

[
\boxed{
D=
\frac{\min(s_L,s_R)}{s_*}
}
]

que es el `drop_factor`.

Por ejemplo,

[
D=1
]

significa que practicamente no hay valle.

Pero

[
D=10^5
]

significa que el valor singular cae unas cien mil veces en el candidato.

En tus BICs obtuvimos cosas como

[
10^5,;10^6,
]

que son minimos espectrales muy claros.

---

# 3. `parity residual`

Esta es una de las nuevas variables mas importantes.

La usamos para verificar que la **solucion propia** encontrada por BEM tiene la simetria esperada.

No estamos comprobando solamente la simetria del obstaculo:

[
X(t),Y(t).
]

Tambien comprobamos la simetria del vector propio

[
\mathbf u.
]

## En el caso (iii), simetria respecto al eje (x)

El BIC que encontramos pertenece al sector impar respecto a la coordenada transversal:

[
u(x,-y)\approx-u(x,y).
]

Por eso construimos un residual impar.

Idealmente:

[
u(\theta)+u(2\pi-\theta)=0.
]

Entonces medimos algo parecido a

[
\boxed{
R\_{\rm odd}
===========

\frac{
|u(\theta)+u(2\pi-\theta)|
}{
|u|
}
}
]

Si

[
R_{\rm odd}=0,
]

la funcion es perfectamente impar.

Si obtenemos

[
R_{\rm odd}\sim10^{-10},
]

como en tus resultados, significa que numericamente la solucion tiene la paridad esperada con muchisima precision.

---

# 4. ¿Por que aparecia `odd_res = 2` en candidatos falsos?

Porque esos candidatos no pertenecian al sector impar.

Supongamos que la funcion fuera par:

[
u(-y)=u(y).
]

Entonces en el residual impar:

[
u(y)+u(-y)=2u(y).
]

Por tanto,

[
R_{\rm odd}\approx2.
]

Eso explica perfectamente valores como

```text
odd_res = 2.000
```

que vimos en los falsos candidatos.

No es un error.

Significa:

[
\boxed{\text{esta solucion no tiene la paridad requerida}.}
]

---

# 5. `yaxis_even_residual` y `yaxis_odd_residual`

En el caso (iv), con simetria respecto al eje (y), calculamos ambos sectores.

Uno mide cuanto se parece el vector propio a una funcion par:

[
\boxed{
R\_{\rm even}
============

\frac{
|u(t)-u(-t)|
}{
|u|
}
}
]

y el otro cuanto se parece a una funcion impar:

[
\boxed{
R\_{\rm odd}
===========

\frac{
|u(t)+u(-t)|
}{
|u|
}.
}
]

Si obtenemos por ejemplo

[
R_{\rm even}=2.5\times10^{-11},
]

[
R_{\rm odd}=2,
]

eso nos dice:

[
\boxed{
\text{el modo es claramente par}.
}
]

Que es justamente lo que estabamos viendo en el caso (iv).

---

# 6. `Rprop` o `propagating_ratio`

Esta es probablemente **la variable nueva mas importante para el Teorema 2.3**.

Recuerda que estamos en

[
\Lambda_1<k^2<\Lambda_2.
]

En ese intervalo el primer modo transversal ya puede propagarse.

El primer modo transversal para Dirichlet es esencialmente

[
\phi_1(y)
=========

\cos\left(\frac{\pi y}{2b}\right).
]

Entonces una solucion general lejos del obstaculo puede contener una componente del tipo

[
C_1e^{iqx}\phi_1(y).
]

Si

[
C_1\neq0,
]

la energia puede propagarse hacia el infinito.

Pero un BIC debe permanecer atrapado.

Por tanto necesitamos

[
\boxed{C_1=0.}
]

Numericamente no esperamos cero exacto, asi que reconstruimos el campo lejos del obstaculo y calculamos cuanto de ese campo pertenece al primer modo propagante.

Definimos aproximadamente

[
C_1
===

\frac{
\langle u,\phi_1\rangle
}{
|\phi_1|^2
}.
]

Luego comparamos el tamaño de esa proyeccion con el tamaño total del campo.

Eso produce:

[
\boxed{
R\_{\rm prop}
============

\frac{
|\text{componente propagante}|
}{
|\text{campo total}|
}.
}
]

---

# 7. ¿Como interpreto (R\_{\rm prop})?

Si

[
R_{\rm prop}\approx1,
]

significa que casi todo el campo pertenece al canal propagante.

Eso **no es un BIC**.

Si

[
R_{\rm prop}=10^{-2},
]

hay aproximadamente una fraccion del (1%).

Si

[
R_{\rm prop}=10^{-4},
]

la fraccion es aproximadamente

[
0.01%.
]

Y si

[
R_{\rm prop}=10^{-6},
]

el acoplamiento es practicamente nulo.

En el caso de simetria (x) obtuvimos cosas como

[
R_{\rm prop}\sim10^{-6},
]

que es espectacular.

En el caso de simetria (y), con solamente

[
a=\varepsilon a_1
]

de orden dominante, obtuvimos aproximadamente

[
10^{-4}-10^{-3}.
]

Sigue siendo pequeño, pero no tan pequeño.

Y justamente eso motivo el nuevo experimento de optimizar (a).

---

# 8. `propagating_ratio_left` y `propagating_ratio_right`

Calculamos (R\_{\rm prop}) en ambos lados del obstaculo.

Por ejemplo:

[
x=-2b,
]

y

[
x=+2b.
]

Por eso tenemos

```text
propagating_ratio_left
propagating_ratio_right
```

Idealmente ambos deben ser pequeños:

[
R_{\rm prop}^{(-)}
\approx0,
\qquad
R_{\rm prop}^{(+)}
\approx0.
]

Despues solemos usar

[
\boxed{
R\_{\rm prop}
============

\max
\left(
R*{\rm prop}^{(-)},
R*{\rm prop}^{(+)}
\right)
}
]

como criterio conservador.

Es decir: tomamos el peor lado.

---

# 9. ¿Por que no basta con (s\_{\min}(A))?

Esta es una distincion fundamental.

Podriamos encontrar

[
s_{\min}(A)\approx0
]

dentro de

[
\Lambda_1<k^2<\Lambda_2,
]

pero eso por si solo indica una solucion espectral del sistema.

No necesariamente nos dice que sea **no radiante**.

Por eso para el Teorema 2.3 hacemos:

[
\boxed{
s_{\min}(A)\ll1
}
]

**y ademas**

[
\boxed{
R_{\rm prop}\ll1.
}
]

La primera dice:

> hay un candidato espectral.

La segunda dice:

> ese candidato no se acopla significativamente al canal abierto.

Juntas son mucho mas convincentes como validacion de un BIC.

---

# 10. `mesh_change`

Esta mide si el resultado cambia cuando aumentamos (M).

Por ejemplo:

[
M=16,24,32,40,48.
]

Para cada uno obtenemos

[
\sigma_{\rm BEM}^{(M)}.
]

Luego calculamos algo como

[
\boxed{
E_M
===

\frac{
|\sigma*M-\sigma*{M-1}|
}{
|\sigma_M|
}.
}
]

Queremos

[
E_M\rightarrow0.
]

Cuando ves:

```text
mesh_change = 0.000%
```

significa que el resultado practicamente no cambia al refinar la frontera.

Eso es evidencia de convergencia numerica.

---

# 11. `kb_shift`

Esta aparece sobre todo cuando examinamos candidatos adicionales.

Supongamos que con

[
M=16
]

encontramos un supuesto modo en

[
kb_{16}.
]

Lo refinamos con

[
M=32
]

y obtenemos

[
kb_{32}.
]

Definimos

[
\boxed{
\Delta kb
=========

|kb*{32}-kb*{16}|.
}
]

Si el candidato es fisico, esperamos que su posicion sea relativamente estable.

Si cambia mucho cuando modificamos (M), probablemente era un artefacto numerico.

---

# 12. `persistent_nonradiating`

Esta es una variable booleana:

```python
True / False
```

que resume varias condiciones para un candidato adicional.

Conceptualmente quiere decir:

[
\boxed{
\text{¿este candidato sobrevive al refinamiento y sigue siendo no radiante?}
}
]

Para que sea `True`, normalmente queremos simultaneamente:

[
s_{\min}(A)\ll1,
]

un `drop_factor` grande,

[
kb_{16}\approx kb_{32},
]

y

[
R_{\rm prop}\ll1.
]

Si falla cualquiera de las condiciones importantes, queda:

```text
persistent_nonradiating = False
```

como todos los falsos minimos cercanos a (\Lambda_2).

---

# 13. `expected_bic_resolved`

Esta variable responde a:

> ¿Encontramos numericamente el BIC que el teorema predice cerca de (k\_{\rm asym})?

Combina principalmente:

- minimo interior;
- (s\_{\min}(A)) pequeño;
- `drop_factor` grande;
- estabilidad con (M).

Entonces:

```text
expected_bic_resolved = True
```

no significa todavia necesariamente que hayamos pasado todos los tests del BIC.

Solo significa que **la rama espectral esperada fue resuelta**.

---

# 14. `nonradiating_bic_verified`

Esta agrega precisamente la condicion fisica faltante:

[
\boxed{R_{\rm prop}\ll1.}
]

Entonces conceptualmente:

[
\texttt{expected_bic_resolved}

- R\_{\rm prop}\ll1
  ]

nos da:

```text
nonradiating_bic_verified = True
```

En el caso (iii) tambien usamos la paridad como una comprobacion fuerte adicional.

---

# 15. `additional_resolved_bics`

Esta es simplemente la cantidad de BICs adicionales encontrados fuera del esperado.

Queremos:

[
\boxed{
\texttt{additional_resolved_bics}=0.
}
]

Porque el Teorema 2.3 dice que existe un **unico** valor propio.

Si algun dia obtenemos:

```text
additional_resolved_bics = 1
```

tendriamos:

[
N_{\rm BIC}=2.
]

Y eso seria precisamente el tipo de fenomeno que queremos investigar cuando salgamos del regimen (\varepsilon\ll1).

---

# 16. `uniqueness_screen_passed`

Es:

```text
True
```

si no encontramos BICs adicionales.

Conceptualmente:

[
\boxed{
\texttt{uniqueness_screen_passed}
\iff
N_{\rm additional}=0.
}
]

---

# 17. `unique_bic_verified`

Combina todo:

[
\boxed{
\begin{aligned}
&\text{BIC esperado resuelto}\
&+\text{BIC no radiante}\
&+\text{ningun BIC adicional}
\end{aligned}
}
]

y entonces declara:

```text
unique_bic_verified = True
```

Es nuestra validacion numerica de la afirmacion:

[
\boxed{
\text{existe un unico BIC en }[\Lambda_1,\Lambda_2).
}
]

---

# 18. `relative_error_sigma`

Esta ya la conoces, pero ahora adquiere especial importancia.

Tenemos

[
\sigma\_{\rm asym}
=================

\varepsilon^2\frac{\pi^3}{b^3}\mu
]

y

[
\sigma\_{\rm BEM}
================

\sqrt{\Lambda*2-k*{\rm BEM}^2}.
]

Entonces:

[
\boxed{
E\_\sigma
========

\frac{
|\sigma*{\rm BEM}-\sigma*{\rm asym}|
}{
|\sigma\_{\rm asym}|
}.
}
]

Esta mide **la precision de la formula asintotica**.

No mide existencia.

Por eso podemos tener perfectamente:

```text
unique_bic_verified = True
asymptotic_agreement_verified = False
```

como ocurre para (\varepsilon=0.11).

Eso significa:

[
\boxed{\text{el BIC existe}}
]

pero

[
\boxed{\text{la formula lider ya no lo localiza con error menor al 5%}.}
]

---

# 19. Variables nuevas del experimento con (a\_{\rm BIC})

En la ultima version agregamos otro grupo.

Tenemos primero:

[
\boxed{
a\_{\rm asym}
============

\varepsilon a_1
}
]

que es la prediccion de primer orden.

Luego optimizamos:

[
\boxed{
a\_{\rm BIC}
===========

\arg\min*a R*{\rm prop}(a).
}
]

Esto representa la posicion numerica que produce la mejor cancelacion del canal propagante.

Despues definimos

[
\boxed{
\Delta a
========

a*{\rm BIC}-a*{\rm asym}.
}
]

Y finalmente:

[
\boxed{
C_2(\varepsilon)
================

\frac{\Delta a}{\varepsilon^2}.
}
]

¿Por que?

Porque el teorema dice

[
a\_{\rm BIC}
===========

\varepsilon a_1

- O(\varepsilon^2).
  ]

Si escribimos mas explicitamente

[
a\_{\rm BIC}
===========

\varepsilon a_1

- C_2\varepsilon^2
  +\cdots,
  ]

entonces

[
\frac{
a_{\rm BIC}-\varepsilon a_1
}{
\varepsilon^2
}
\approx C_2.
]

Ese es precisamente el coeficiente que intentamos observar numericamente.

---

## El mapa mental final

Cuando veas las variables del codigo, separalas en cuatro grupos:

[
\boxed{\text{1. Espectro}}
]

[
kb,\quad
\sigma_{\rm BEM},\quad
s_{\min}(A),\quad
\text{drop}
]

responden:

> ¿hay un valor propio?

Luego:

[
\boxed{\text{2. Simetria}}
]

[
R_{\rm even},\quad
R_{\rm odd}
]

responden:

> ¿la solucion tiene la paridad que esperamos?

Luego:

[
\boxed{\text{3. Radiacion}}
]

[
R_{\rm prop}^{L},
\quad
R_{\rm prop}^{R}
]

responden:

> ¿el estado esta acoplado al canal propagante o realmente esta atrapado?

Luego:

[
\boxed{\text{4. Convergencia y unicidad}}
]

[
\text{mesh change},
\quad
kb\text{ shift},
\quad
N_{\rm additional}
]

responden:

> ¿es estable numericamente y es realmente el unico?

Y para el caso (iv) agregamos una quinta capa:

[
\boxed{\text{5. Tuning geometrico}}
]

[
a_{\rm asym},
\quad
a_{\rm BIC},
\quad
\Delta a,
\quad
\frac{\Delta a}{\varepsilon^2}
]

que responde:

> ¿cual es realmente la posicion del obstaculo que produce el BIC, y como se relaciona con (a=\varepsilon a_1+O(\varepsilon^2))?

Esa es la logica completa de todas esas variables.
