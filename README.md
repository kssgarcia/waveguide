# Computational Wave Equations and Helmholtz Solvers

Este proyecto contiene implementaciones para resolver ecuaciones de ondas, ecuaciones de Helmholtz y métodos numéricos relacionados utilizando Python y MATLAB. El código se enfoca en métodos de elementos finitos (FEM), métodos de elementos de frontera (BEM), análisis de Fourier y validación contra soluciones analíticas.

## 📋 Tabla de Contenidos

- [Descripción General](#descripción-general)
- [Estructura del Proyecto](#estructura-del-proyecto)
- [Componentes Principales](#componentes-principales)
- [Instalación](#instalación)
- [Uso](#uso)
- [Dependencias](#dependencias)
- [Documentación](#documentación)

## 📖 Descripción General

Este proyecto implementa métodos numéricos para resolver problemas de propagación de ondas en guías de ondas bidimensionales con obstáculos. El enfoque principal está en:

1. **Ecuación de Helmholtz**: -∇²u - k²u = f(x,y)
2. **Ecuación de Onda**: ∂²u/∂t² = c²∇²u
3. **Ecuación de Poisson**: -∇²u = f(x,y)
4. **Modos Atrapados**: Soluciones localizadas en guías de ondas con obstáculos

El proyecto incluye validaciones contra resultados de referencia (Linton 1992, Linton 2002) y análisis de estabilidad numérica.

## 📁 Estructura del Proyecto

```
waveguide/
├── src/
│   ├── fem_core/                    # Núcleo de Elementos Finitos
│   │   ├── __init__.py
│   │   ├── numerical.py            # Métodos numéricos generales
│   │   ├── loc_matrices.py         # Matrices locales (rigidez, masa, convección)
│   │   ├── simple_poisson.py        # Solver de Poisson
│   │   ├── numerical_nh_bc.py       # Condiciones de frontera no homogéneas BC
│   │   ├── numerical_nh_dir.py      # Condiciones de frontera no homogéneas Dirichlet
│   │   └── numerical_nh_dir_2.py   # Condiciones de frontera no homogéneas Dirichlet v2
│   │
│   ├── fourier_analysis/            # Análisis de series de Fourier
│   │   ├── __init__.py
│   │   ├── fourier.py              # Series de Fourier
│   │   └── fourier_nh_bc.py        # Fourier con condiciones de frontera no homogéneas
│   │
│   ├── green_functions/             # Funciones de Green y validaciones
│   │   ├── __init__.py
│   │   ├── greens_methods.py        # Métodos generales de función de Green
│   │   ├── greens_waveguides.py     # Green para guías de ondas
│   │   ├── lattice_sums.py          # Sumas de red
│   │   ├── ellip.py                 # Funciones elípticas
│   │   ├── trapped_modes_solver.py  # Solver de modos atrapados
│   │   ├── obstacles_polar.py       # Obstáculos en coordenadas polares
│   │   ├── helmholtz_corrected.py   # Corrección de Helmholtz
│   │   ├── analytic_k2.py           # Solución analítica de k²
│   │   ├── linton1992.py            # Validación Linton 1992
│   │   ├── linton1992_circle.py     # Validación Linton 1992 - círculo
│   │   ├── linton1992_circle_2.py   # Validación Linton 1992 - círculo v2
│   │   ├── linton1992_circle_lattice.py   # Validación Linton 1992 - círculo lattice
│   │   ├── linton1992_square.py     # Validación Linton 1992 - cuadrado
│   │   ├── linton1992_square_lattice.py  # Validación Linton 1992 - cuadrado lattice
│   │   ├── linton1992_ellipse.py    # Validación Linton 1992 - elipse
│   │   ├── linton1992_neumann.py    # Validación Linton 1992 - Neumann
│   │   ├── linton2002_bisection.py  # Bisección Linton 2002
│   │   ├── linton2002_circle_lattice.py   # Linton 2002 - círculo lattice
│   │   ├── linton2002_circle_mesh_plot.py # Visualización malla circular
│   │   ├── linton2002_plot_determinant.py # Gráficos de determinante
│   │   ├── linton2002_validation_theorem.py    # Validación teorema
│   │   ├── linton2002_validation_theorem_2.py  # Validación teorema v2
│   │   ├── plot_kd_det_circle.py    # Gráfico kd determinante círculo
│   │   ├── plot_comparison.py       # Gráficos comparativos
│   │   ├── test_beta.py             # Pruebas de beta
│   │   ├── test_nested_functions.py # Pruebas de funciones anidadas
│   │   └── results.dat              # Datos de resultados
│   │
│   ├── grid_generation/             # Generación de mallas
│   │   ├── __init__.py
│   │   ├── gen_p1grid.py           # Generación de mallas P1
│   │   ├── drive_gen_p1grid.py     # Driver para generación P1
│   │   ├── mesh_reading.py         # Lectura de mallas
│   │   ├── mesh_create_center_xy.py # Crear malla con centro xy
│   │   ├── mesh_create_center_y.py  # Crear malla con centro y
│   │   ├── mesh_create_center_y_no_hole.py  # Malla sin agujero
│   │   ├── mesh_create_center_y_pertube.py  # Malla con perturbación
│   │   └── gen_mesh_center_y.py    # Generar malla centro y
│   │
│   ├── helmholtz/                   # Solvers de ecuación de Helmholtz
│   │   ├── __init__.py
│   │   ├── helmholtz_solution.py    # Solver principal FEM
│   │   ├── helmholtz_solution_stability.py  # Análisis de estabilidad
│   │   ├── DtN_base.py             # Condiciones Dirichlet-to-Neumann
│   │   ├── lattice_sums.py         # Sumas de red para función de Green
│   │   ├── validation_kb.py        # Validación con parámetro kb
│   │   ├── bem_potential_plot.py   # Visualización de potencial BEM
│   │   ├── utils/
│   │   │   ├── __init__.py
│   │   │   └── utils.py            # Utilidades auxiliares
│   │   ├── solutions/              # Archivos de solución
│   │   │   ├── Helmholtz_dolfinx_doc_example.h5
│   │   │   ├── Helmholtz_dolfinx_doc_example.xdmf
│   │   │   ├── hemlholtz_meshio_dolfinx.h5
│   │   │   └── hemlholtz_meshio_dolfinx.xdmf
│   │   └── old_files/              # Código legacy
│   │       ├── DtN.py
│   │       ├── DtN_newton.py
│   │       ├── DtN_newton_correction.py
│   │       ├── DtN_newton_sol_evaluation.py
│   │       ├── DtN_strip.py
│   │       ├── find_k.py
│   │       ├── find_k_dirichlet.py
│   │       ├── full_dirichlet.py
│   │       ├── generate_mesh_square_obstacle.py
│   │       ├── helmholtz_dolfinx_doc_example.py
│   │       ├── helmholtz_given_mesh.py
│   │       ├── helmholtz_meshio_dolfinx.py
│   │       ├── read_mesh.py
│   │       ├── simple_helmholtz_eig.py
│   │       ├── simple_helmholtz_eig_obstacle.py
│   │       ├── test_dolpinx.py
│   │       ├── test_github.py
│   │       └── check_integral_phi_varphi.py
│   │
│   ├── laplace/                     # Método de elementos de frontera (BEM)
│   │   ├── bem.py                  # Implementación BEM general
│   │   ├── bem_circle.py           # BEM para círculo
│   │   ├── gauss_legendre.py       # Cuadratura de Gauss-Legendre
│   │   └── plot_output.py          # Visualización de salida
│   │
│   ├── matlab_implementations/       # Implementaciones en MATLAB
│   │   ├── README.md
│   │   ├── gen_p1grid.m            # Generación de mallas P1
│   │   ├── drive_gen_p1grid.m      # Driver para generación
│   │   ├── ellip.m                 # Funciones elípticas
│   │   ├── locA.m                  # Matriz de rigidez
│   │   ├── locB.m                  # Matriz de convección
│   │   └── locM.m                  # Matriz de masa
│   │
│   ├── theorem_1/                   # Validación del Teorema 1
│   │   ├── README.md
│   │   ├── validation.py           # Validación principal
│   │   ├── validation_kb.py        # Validación con kb
│   │   ├── validation_a_vs_eps.py  # Análisis a vs ε
│   │   ├── validation_multiple.py   # Validación múltiple
│   │   ├── validation_single.py     # Validación caso único
│   │   ├── lattice_sums.py         # Sumas de red
│   │   ├── kevin.csv               # Datos de validación
│   │   └── kevin_kb.csv            # Datos de validación kb
│   │
│   ├── theorem_2.3/                 # Implementación del Teorema 2.3
│   │   ├── dipole_theorem.py       # Análisis de dipolo
│   │   ├── gauss_legendre.py       # Cuadratura numérica
│   │   ├── lattice_sums.py         # Sumas de red
│   │   ├── bem_theorem_2.3_3.py    # BEM teorema 2.3
│   │   ├── bem_theorem_2.3_3_mult.py  # BEM teorema 2.3 mult
│   │   ├── better_bem_therem_2.3_3_mult.py  # BEM mejorado mult
│   │   ├── old_scripts/            # Scripts legacy (no listar)
│   │   └── reviewer_validation_output/  # Resultados para revisores
│   │
│   └── wave_equation/               # Solvers de ecuación de onda
│       ├── __init__.py
│       ├── simple_wave.py          # Solver básico de onda
│       ├── simple_wave_eig.py      # Análisis de eigenvalores
│       ├── error_wave_analytic_numerical.py  # Análisis de error
│       ├── wave_3d_animation.gif   # Animación 3D
│       └── wave_analytical/        # Soluciones analíticas
│           ├── wave_analytical_1D_neumann.py
│           ├── wave_analytical_2D_dirichlet.py
│           ├── wave_analytical_2D_neumann.py
│           ├── wave_fin_diff_1D_neumann.py
│           ├── wave_fin_diff_2D_dirichlet.py
│           ├── wave_fin_diff_2D_neumann.py
│           └── f_tanh.py
│
├── meshes/                          # Archivos de malla (.msh)
│   ├── mesh_no_hole.msh
│   ├── mesh_with_hole.msh
│   ├── mesh_with_perturbed_hole.msh
│   └── mesh_with_perturbed_hole_refine.msh
│
├── solutions/                       # Archivos de solución (.h5, .xdmf)
│   ├── Helmholtz_dolfinx_doc_example.h5
│   ├── Helmholtz_dolfinx_doc_example.xdmf
│   ├── hemlholtz_meshio_dolfinx.h5
│   └── hemlholtz_meshio_dolfinx.xdmf
│
├── docs/                            # Documentación LaTeX
│   ├── report.tex
│   ├── reporte2.tex
│   ├── report2.tex
│   ├── DtN_matrix.tex
│   ├── DtN_matrix_2.tex
│   ├── DtN_justification.tex
│   ├── coefficients_mu_nu.tex
│   ├── laplace.tex
│   ├── solucion.tex
│   ├── square_polar.tex
│   ├── circle_polar.tex
│   ├── waveguide_tikz.tex
│   ├── waveguide_rectangle_obstacle.tex
│   ├── waveguide_cardiod.tex
│   ├── linton1992_figure.tex
│   ├── references.bib
│   ├── Guia_Plantilla_Latex_ColoquioUMNG_2025_Esp.tex
│   └── AG_Guia_Plantilla_Latex_ColoquioUMNG_2025_Esp.tex
│
├── semillero/                       # Material del semillero (MATLAB)
│   ├── movie_example.m
│   ├── movie_example_2.m
│   ├── onda_dAlembert.m
│   ├── test.m
│   ├── programas/
│   │   ├── semillero2.m
│   │   ├── semillero3.m
│   │   ├── semillero4.m
│   │   ├── semillero5.m
│   │   ├── semillero6.m
│   │   ├── SEMILELRO.m
│   │   ├── fourier.m
│   │   ├── fun.m
│   │   ├── f.m
│   │   └── plot_numerical.m
│   └── poster_ondas/
│       ├── poster_ondas.tex
│       ├── references.bib
│       └── beamerthemesharelatex.sty
│
├── examples/                        # Scripts de ejemplo
│   ├── old_code.py
│   └── ocr.py
│
└── document.tex                     # Documento principal
```

## 🔧 Componentes Principales

### 1. Núcleo FEM (`src/fem_core/`) — **IMPORTANTE**

Implementación de elementos finitos desde cero:

- **`numerical.py`**: Métodos numéricos generales
- **`loc_matrices.py`**: Ensamblaje de matrices elementales (rigidez, masa, convección)
- **`simple_poisson.py`**: Solver de ecuación de Poisson
- **`numerical_nh_bc.py`**, **`numerical_nh_dir.py`**, **`numerical_nh_dir_2.py`**: Condiciones de frontera no homogéneas

### 2. Funciones de Green (`src/green_functions/`) — **IMPORTANTE**

Implementación de funciones de Green para guías de ondas y validaciones:

- **`greens_methods.py`**: Métodos generales de función de Green
- **`greens_waveguides.py`**: Green específica para guías de ondas
- **`lattice_sums.py`**: Sumas de red usando funciones especiales
- **`ellip.py`**: Funciones elípticas
- **`trapped_modes_solver.py`**: Solver especializado para modos atrapados
- **`obstacles_polar.py`**: Obstáculos en coordenadas polares
- **`helmholtz_corrected.py`**: Corrección de Helmholtz

**Validaciones Linton (1992):**

- **`linton1992_circle.py`**, **`linton1992_circle_lattice.py`**, **`linton1992_circle_2.py`**: Círculos
- **`linton1992_ellipse.py`**: Elipses
- **`linton1992_square.py`**, **`linton1992_square_lattice.py`**: Cuadrados
- **`linton1992_neumann.py`**: Condiciones Neumann

**Validaciones Linton (2002):**

- **`linton2002_circle_lattice.py`**: Redes de círculos
- **`linton2002_bisection.py`**: Método de bisección
- **`linton2002_validation_theorem.py`**, **`linton2002_validation_theorem_2.py`**: Validación de teoremas
- **`linton2002_plot_determinant.py`**, **`linton2002_circle_mesh_plot.py`**: Visualización
- **`plot_kd_det_circle.py`**, **`plot_comparison.py`**: Gráficos comparativos

### 3. Helmholtz (`src/helmholtz/`) — **IMPORTANTE**

Solvers de la ecuación de Helmholtz con DOLFINx:

- **`helmholtz_solution.py`**: Solver FEM principal
- **`helmholtz_solution_stability.py`**: Análisis de estabilidad numérica
- **`DtN_base.py`**: Condiciones de frontera Dirichlet-to-Neumann (DtN)
- **`lattice_sums.py`**: Sumas de red para función de Green
- **`validation_kb.py`**: Validación con parámetro kb
- **`bem_potential_plot.py`**: Visualización de potencial BEM
- **`utils/utils.py`**: Utilidades auxiliares
- **`old_files/`**: Código legacy (DtN_newton.py, find_k.py, etc.)
- **`solutions/`**: Archivos de solución guardados

### 4. Teorema 1 (`src/theorem_1/`) — **IMPORTANTE**

Validación del Teorema 2.1 sobre modos atrapados en guías de ondas:

- **`validation.py`**: Validación principal
- **`validation_kb.py`**: Validación con parámetro kb
- **`validation_a_vs_eps.py`**: Análisis de relación a vs ε
- **`validation_multiple.py`**: Validación para múltiples valores
- **`validation_single.py`**: Validación para caso único
- **`lattice_sums.py`**: Sumas de red específicas
- **`kevin.csv`**, **`kevin_kb.csv`**: Datos de validación

### 5. Teorema 2.3 (`src/theorem_2.3/`) — **IMPORTANTE**

Implementación del Teorema 2.3 sobre soluciones tipo dipolo:

- **`fem_theorem_2.3*.py`**: Implementación FEM
- **`bem_theorem_2.3*.py`**: Implementación BEM
- **`dipole_theorem.py`**: Análisis de dipolo
- **`gauss_legendre.py`**: Cuadratura numérica
- **`lattice_sums.py`**: Sumas de red
- **`old_scripts/`**: Scripts legacy (mencionado, no listar archivos internos)
- **`reviewer_validation_output/`**: Resultados de validación para revisores

### 6. Otros Módulos (Resumidos)

- **`src/fourier_analysis/`**: Series de Fourier y análisis espectral (fourier.py, fourier_nh_bc.py)
- **`src/grid_generation/`**: Generación de mallas triangulares (gen_p1grid.py, mesh_create_*.py, mesh_reading.py)
- **`src/laplace/`**: Método BEM para ecuación de Laplace (bem.py, bem_circle.py, gauss_legendre.py)
- **`src/matlab_implementations/`**: Implementaciones MATLAB del núcleo FEM (gen_p1grid.m, locA.m, locM.m, ellip.m)
- **`src/wave_equation/`**: Solvers de ecuación de onda (simple_wave.py, simple_wave_eig.py, error_wave_analytic_numerical.py, wave_analytical/)

## 🚀 Instalación

### Requisitos Previos

- Python 3.8+
- pip o conda
- (Opcional) MATLAB R2018b+ para scripts MATLAB

### Dependencias Python

```bash
pip install numpy scipy matplotlib
pip install meshio
pip install dolfinx  # FEniCSx
pip install mpi4py petsc4py
```

### Instalación con Conda (Recomendado para DOLFINx)

```bash
conda create -n ondas python=3.10
conda activate ondas
conda install -c conda-forge numpy scipy matplotlib fenics-dolfinx meshio
```

## 📊 Uso

```bash
# Helmholtz
cd src/helmholtz/
python helmholtz_solution.py

# Modos atrapados
cd src/green_functions/
python trapped_modes_solver.py

# Validación Teorema 1
cd src/theorem_1/
python validation.py

# Ecuación de onda
cd src/wave_equation/
python simple_wave.py
```

## 📦 Dependencias

| Paquete | Propósito |
| --------- | ----------- |
| `numpy`, `scipy` | Computación numérica |
| `matplotlib` | Visualización |
| `dolfinx` | Elementos finitos (FEniCSx) |
| `meshio` | I/O de mallas |
| `mpi4py`, `petsc4py` | Paralelización |

## 📚 Documentación

La documentación LaTeX está en `docs/`:

- `report.tex` — Reporte principal
- `DtN_matrix.tex`, `DtN_justification.tex` — Documentación DtN
- `references.bib` — Referencias bibliográficas

## 🤝 Autores

- Kevin Santiago Sepúlveda García
- Colaboradores del semillero

## 📄 Licencia

Material académico para investigación en propagación de ondas y métodos numéricos.

## 🎯 Estado del Proyecto

- ✅ Solver FEM Helmholtz
- ✅ Condiciones DtN
- ✅ Funciones de Green con lattice sums
- ✅ Validación Linton (1992/2002)
- ✅ Solver de modos atrapados
- ✅ Análisis de estabilidad numérica
- ✅ Semillero MATLAB

---

**Última actualización:** 2025
