# Computational Wave Equations and Helmholtz Solvers

Este proyecto contiene implementaciones para resolver ecuaciones de ondas, ecuaciones de Helmholtz y métodos numéricos relacionados utilizando Python y MATLAB. El código se enfoca en métodos de elementos finitos (FEM), métodos de elementos de frontera (BEM), análisis de Fourier y validación contra soluciones analíticas.

## 📋 Tabla de Contenidos

- [Descripción General](#descripción-general)
- [Estructura del Proyecto](#estructura-del-proyecto)
- [Componentes Principales](#componentes-principales)
- [Teoremas Implementados](#teoremas-implementados)
- [Instalación](#instalación)
- [Uso](#uso)
- [Dependencias](#dependencias)
- [Documentación](#documentación)
- [Contribuciones](#contribuciones)

## 📖 Descripción General

Este proyecto implementa métodos numéricos para resolver problemas de propagación de ondas en guías de ondas bidimensionales con obstáculos. El enfoque principal está en:

1. **Ecuación de Helmholtz**: -∇²u - k²u = f(x,y)
2. **Ecuación de Onda**: ∂²u/∂t² = c²∇²u
3. **Ecuación de Poisson**: -∇²u = f(x,y)
4. **Modos Atrapados**: Soluciones localizadas en guías de ondas con obstáculos

El proyecto incluye validaciones contra resultados de referencia (Linton 1992, Linton 2002) y análisis de estabilidad numérica.

## 📁 Estructura del Proyecto

```
ondas/
├── src/                              # Código fuente principal
│   ├── helmholtz/                    # Solvers de ecuación de Helmholtz
│   │   ├── DtN_base.py              # Implementación de condiciones DtN
│   │   ├── helmholtz_solution.py    # Solver principal FEM
│   │   ├── helmholtz_solution_stability.py  # Análisis de estabilidad
│   │   ├── lattice_sums.py          # Sumas de red para función de Green
│   │   ├── validation_kb.py         # Validación con parámetro kb
│   │   └── utils/                   # Utilidades auxiliares
│   │
│   ├── green_functions/              # Funciones de Green y validaciones
│   │   ├── greens_methods.py        # Métodos de función de Green
│   │   ├── greens_waveguides.py     # Green para guías de ondas
│   │   ├── lattice_sums.py          # Sumas de red
│   │   ├── trapped_modes_solver.py  # Solver de modos atrapados
│   │   ├── linton1992*.py           # Validación Linton (1992)
│   │   ├── linton2002*.py           # Validación Linton (2002)
│   │   ├── obstacles_polar.py       # Obstáculos en coordenadas polares
│   │   └── ellip.py                 # Funciones elípticas
│   │
│   ├── wave_equation/                # Ecuación de onda
│   │   ├── simple_wave.py           # Solver básico de onda
│   │   ├── simple_wave_eig.py       # Análisis de eigenvalores
│   │   ├── error_wave_analytic_numerical.py  # Análisis de error
│   │   └── wave_analytical/         # Soluciones analíticas
│   │
│   ├── fem_core/                     # Núcleo FEM
│   │   ├── numerical.py             # Métodos numéricos principales
│   │   ├── loc_matrices.py          # Ensamblaje de matrices locales
│   │   ├── simple_poisson.py        # Solver de Poisson
│   │   └── numerical_nh_*.py        # Condiciones de frontera no homogéneas
│   │
│   ├── fourier_analysis/             # Análisis de Fourier
│   │   ├── fourier.py               # Series de Fourier
│   │   └── fourier_nh_bc.py         # Fourier con BC no homogéneas
│   │
│   ├── grid_generation/              # Generación de mallas
│   │   ├── gen_p1grid.py            # Mallas triangulares P1
│   │   ├── mesh_create.py           # Creación de mallas
│   │   └── mesh_reading.py          # Lectura de mallas
│   │
│   ├── laplace/                      # Ecuación de Laplace (BEM)
│   │   ├── bem.py                   # Método de elementos de frontera
│   │   ├── bem_circle.py            # BEM para círculo
│   │   └── gauss_legendre.py        # Cuadratura de Gauss-Legendre
│   │
│   ├── theorem_1/                    # Implementación Teorema 1
│   │   ├── validation.py            # Validación principal
│   │   ├── validation_kb.py         # Validación con kb
│   │   ├── validation_a_vs_eps.py   # Análisis a vs ε
│   │   └── lattice_sums.py          # Sumas de red específicas
│   │
│   ├── theorem_2.3/                  # Implementación Teorema 2.3
│   │   ├── fem_theorem_2.3.py       # FEM para teorema 2.3
│   │   ├── fem_theorem_2.3_dtn.py   # FEM con DtN
│   │   ├── bem_theorem_2.3.py       # BEM para teorema 2.3
│   │   ├── dipole_theorem.py        # Teorema de dipolo
│   │   └── gauss_legendre.py        # Cuadratura numérica
│   │
│   └── matlab_implementations/       # Implementaciones MATLAB
│       ├── ellip.m                  # Solver elíptico
│       ├── gen_p1grid.m             # Generación P1
│       └── loc*.m                   # Matrices locales
│
├── meshes/                           # Archivos de malla (.msh, .npz)
├── solutions/                        # Archivos de solución (.h5, .xdmf)
├── docs/                            # Documentación LaTeX
│   ├── report.tex                   # Reporte principal
│   ├── DtN_matrix.tex              # Documentación DtN
│   ├── references.bib              # Referencias bibliográficas
│   └── *.tex                       # Otros documentos
│
├── figures/                         # Figuras generadas
├── semillero/                       # Material del semillero
│   ├── poster_ondas/               # Póster de presentación
│   └── programas/                  # Programas auxiliares
│
└── examples/                        # Scripts de ejemplo
    ├── ocr.py                      # OCR para documentos
    └── old_code.py                 # Código de referencia

```

## 🔧 Componentes Principales

### 1. Solvers de Helmholtz (`src/helmholtz/`)

Implementaciones de la ecuación de Helmholtz con diferentes enfoques:

- **`helmholtz_solution.py`**: Solver FEM principal usando DOLFINx
- **`helmholtz_solution_stability.py`**: Análisis de estabilidad numérica
- **`DtN_base.py`**: Condiciones de frontera Dirichlet-to-Neumann (DtN)
- **`lattice_sums.py`**: Funciones de Green mediante sumas de red
- **`validation_kb.py`**: Validación usando el parámetro kb

**Características:**
- Elementos finitos lineales (P1)
- Condiciones de frontera mixtas (Dirichlet/Neumann/DtN)
- Análisis de convergencia
- Validación contra soluciones analíticas

### 2. Funciones de Green (`src/green_functions/`)

Implementación de funciones de Green para guías de ondas:

- **`greens_methods.py`**: Métodos generales de función de Green
- **`greens_waveguides.py`**: Green específica para guías de ondas
- **`lattice_sums.py`**: Sumas de red usando funciones especiales
- **`trapped_modes_solver.py`**: Solver especializado para modos atrapados
- **`obstacles_polar.py`**: Obstáculos en coordenadas polares

**Validaciones implementadas:**
- Linton (1992): Círculos, elipses, cuadrados
- Linton (2002): Redes de obstáculos circulares
- Análisis de determinantes para eigenvalores

### 3. Ecuación de Onda (`src/wave_equation/`)

Solvers para la ecuación de onda:

- **`simple_wave.py`**: Solver básico con método de diferencias finitas
- **`simple_wave_eig.py`**: Análisis de eigenvalores y modos propios
- **`error_wave_analytic_numerical.py`**: Comparación numérico-analítica
- **`wave_analytical/`**: Soluciones analíticas de referencia

### 4. Núcleo FEM (`src/fem_core/`)

Implementación de elementos finitos desde cero:

- **`numerical.py`**: Métodos numéricos generales
- **`loc_matrices.py`**: Ensamblaje de matrices elementales
  - `locA`: Matriz de rigidez
  - `locB`: Matriz de convección
  - `locM`: Matriz de masa
- **`simple_poisson.py`**: Solver de ecuación de Poisson
- **`numerical_nh_*.py`**: Condiciones de frontera no homogéneas

### 5. Análisis de Fourier (`src/fourier_analysis/`)

Métodos de series de Fourier:

- **`fourier.py`**: Implementación de series de Fourier
- **`fourier_nh_bc.py`**: Fourier con condiciones de frontera no homogéneas

### 6. Generación de Mallas (`src/grid_generation/`)

Utilidades para crear y manipular mallas:

- **`gen_p1grid.py`**: Generación de mallas triangulares P1
- **`mesh_create.py`**: Creación de mallas con obstáculos
- **`mesh_reading.py`**: Lectura/escritura de archivos de malla

### 7. BEM para Laplace (`src/laplace/`)

Método de elementos de frontera:

- **`bem.py`**: Implementación BEM general
- **`bem_circle.py`**: BEM para dominios circulares
- **`gauss_legendre.py`**: Integración numérica

## 📐 Teoremas Implementados

### Teorema 1: Modos Atrapados en Guías de Ondas

**Ubicación:** `src/theorem_1/`

**Descripción:** Validación del Teorema 2.1 sobre existencia de modos atrapados en guías de ondas con obstáculos pequeños.

**Archivos principales:**
- `validation.py`: Validación usando transformación logarítmica s = -2log₁₀(σ)
- `validation_kb.py`: Validación con parámetro kb
- `validation_a_vs_eps.py`: Relación entre altura del obstáculo (a) y parámetro ε
- `validation_multiple.py`: Validación para múltiples valores de a

**Objetivo:** Determinar ε*_máx para el cual el teorema se cumple.

**Referencias:** [1] Artículo base del teorema

### Teorema 2.3: Soluciones Dipolo

**Ubicación:** `src/theorem_2.3/`

**Descripción:** Implementación del Teorema 2.3 sobre soluciones tipo dipolo.

**Archivos principales:**
- `fem_theorem_2.3.py`: Implementación FEM
- `fem_theorem_2.3_dtn.py`: FEM con condiciones DtN
- `bem_theorem_2.3.py`: Implementación BEM
- `dipole_theorem.py`: Análisis de dipolo
- Archivos de análisis: `*_analysis.png`, `*_scaling.png`

## 🚀 Instalación

### Requisitos Previos

- Python 3.8+
- pip o conda
- (Opcional) MATLAB R2018b+ para scripts MATLAB

### Dependencias Python

```bash
# Instalación con pip
pip install numpy scipy matplotlib
pip install meshio
pip install dolfinx  # FEniCSx
pip install mpi4py petsc4py
```

### Instalación con Conda (Recomendado para DOLFINx)

```bash
# Crear entorno
conda create -n ondas python=3.10
conda activate ondas

# Instalar dependencias
conda install -c conda-forge numpy scipy matplotlib
conda install -c conda-forge fenics-dolfinx
conda install -c conda-forge meshio

```

## 📊 Uso

### Ejemplo 1: Solver de Helmholtz Básico

```python
# Navegar a src/helmholtz/
cd src/helmholtz/

# Ejecutar solver principal
python helmholtz_solution.py
```

### Ejemplo 2: Validación Teorema 1

```python
# Navegar a theorem_1
cd src/theorem_1/

# Validación principal
python validation.py


# Validación con kb
python validation_kb.py

# Análisis a vs epsilon
python validation_a_vs_eps.py
```

### Ejemplo 3: Modos Atrapados

```python
# Navegar a green_functions
cd src/green_functions/

# Ejecutar solver de modos atrapados
python trapped_modes_solver.py

# Validación Linton 1992 para círculo
python linton1992_circle.py

# Validación Linton 2002

python linton2002_validation_theorem.py
```

### Ejemplo 4: Ecuación de Onda

```python
# Navegar a wave_equation
cd src/wave_equation/

# Solver básico
python simple_wave.py

# Análisis de error
python error_wave_analytic_numerical.py
```

### Ejemplo 5: Generación de Mallas

```python
# Navegar a grid_generation
cd src/grid_generation/

# Generar malla P1
python gen_p1grid.py

# Crear malla con obstáculo
python mesh_create.py
```

## 📦 Dependencias

### Python

| Paquete | Versión | Propósito |
|---------|---------|-----------|
| `numpy` | ≥1.20 | Computación numérica |
| `scipy` | ≥1.7 | Álgebra lineal, optimización |
| `matplotlib` | ≥3.3 | Visualización |
| `dolfinx` | ≥0.6 | Elementos finitos (FEniCSx) |
| `meshio` | ≥5.0 | I/O de mallas |
| `mpi4py` | ≥3.0 | Paralelización MPI |
| `petsc4py` | ≥3.17 | Álgebra lineal PETSc |
| `slepc4py` | ≥3.17 | Eigenvalores SLEPc |
| `sympy` | ≥1.9 | Matemática simbólica |

### MATLAB

- MATLAB base (versión R2018b o superior)
- (Opcional) Partial Differential Equation Toolbox

## 📚 Documentación

La documentación completa en LaTeX se encuentra en el directorio `docs/`:

- **`report.tex`**: Reporte principal del proyecto
- **`DtN_matrix.tex`**: Documentación de matrices DtN
- **`DtN_justification.tex`**: Justificación del método DtN
- **`references.bib`**: Referencias bibliográficas completas
- **`waveguide_*.tex`**: Diagramas de guías de ondas

### Compilar Documentación

```bash
cd docs/
pdflatex report.tex
bibtex report
pdflatex report.tex
pdflatex report.tex
```

## 🔬 Métodos Implementados

### Métodos de Elementos Finitos (FEM)
- Elementos P1 triangulares
- Ensamblaje de matrices locales (rigidez, masa, convección)
- Condiciones de frontera Dirichlet, Neumann y mixtas
- Eigenvalores y eigenvectores

### Métodos de Elementos de Frontera (BEM)
- Discretización de contornos
- Integración de funciones de Green
- Cuadratura de Gauss-Legendre
- Método directo e indirecto

### Sumas de Red (Lattice Sums)
- Funciones de Green periódicas
- Aceleración de convergencia
- Funciones especiales (Bessel, Hankel)
- Correcciones de singularidad

### Análisis Numérico
- Análisis de convergencia
- Estimación de errores
- Análisis de estabilidad
- Comparación con soluciones analíticas

## 📈 Resultados y Validaciones

El proyecto incluye validaciones contra:

1. **Linton (1992)**: Modos atrapados en guías de ondas
   - Obstáculos circulares, elípticos y cuadrados
   - Diferentes condiciones de frontera

2. **Linton (2002)**: Redes periódicas de obstáculos
   - Análisis de determinantes
   - Búsqueda de eigenvalores

3. **Soluciones Analíticas**: 
   - Ecuación de onda: d'Alembert
   - Helmholtz: Series de Fourier
   - Casos especiales con geometrías simples

## 🤝 Contribuciones

### Autores Principales
- Kevin Santiago Sepúlveda García
- agarz (colaborador principal)

### Cómo Contribuir

1. Fork el repositorio
2. Crea una rama para tu feature (`git checkout -b feature/NuevaCaracteristica`)
3. Commit tus cambios (`git commit -m 'Añadir nueva característica'`)
4. Push a la rama (`git push origin feature/NuevaCaracteristica`)
5. Abre un Pull Request

## 📄 Licencia

Este proyecto es material académico desarrollado para investigación en propagación de ondas y métodos numéricos.

## 📧 Contacto

Para preguntas o colaboraciones:
- Kevin Santiago Sepúlveda García
- Universidad: [Información de contacto]

## 🔗 Referencias Principales

[1] Linton, C. M. (1992). "Trapped modes in waveguides"

[2] Linton, C. M. (2002). "Lattice sums for the Helmholtz equation"

[3] FEniCSx Project: https://fenicsproject.org/

Ver `docs/references.bib` para referencias completas.

## 📝 Notas Adicionales

- Los archivos `.npz` y `.mat` contienen datos precalculados para validaciones
- Las soluciones se guardan en formato HDF5 (.h5) y XDMF para visualización en ParaView
- Los scripts incluyen parámetros configurables al inicio de cada archivo
- Algunas implementaciones tienen versiones alternativas en `old_files/`

## 🎯 Estado del Proyecto

- ✅ Solver FEM Helmholtz
- ✅ Solver BEM Laplace
- ✅ Funciones de Green con lattice sums
- ✅ Validación Teorema 1
- ✅ Validación Teorema 2.3
- ✅ Solver de modos atrapados
- ✅ Ecuación de onda
- 🔄 Análisis de estabilidad numérica (en progreso)
- 🔄 Optimización de rendimiento (en progreso)
- 📋 Interfaz gráfica de usuario (planeado)

---

**Última actualización:** 2025