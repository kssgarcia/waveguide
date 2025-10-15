# Computational Wave Equations and Helmholtz Solvers

This project contains implementations for solving wave equations, Helmholtz equations, and related numerical methods using both Python and MATLAB. The codebase focuses on finite element methods (FEM), Fourier analysis, and validation against analytical solutions.

## Project Structure

```
├── src/                          # Main source code
│   ├── helmholtz/               # Helmholtz equation solvers
│   ├── wave_equation/           # Wave equation implementations
│   ├── fourier_analysis/        # Fourier series and analysis
│   ├── fem_core/               # Core FEM implementations
│   ├── grid_generation/        # Mesh and grid generation
│   ├── validation/             # Validation and benchmark tests
│   └── matlab_implementations/ # MATLAB versions of algorithms
├── meshes/                     # Mesh files (.msh format)
├── solutions/                  # Output solution files (.h5, .xdmf)
├── docs/                      # Documentation and LaTeX files
└── examples/                  # Example scripts and utilities
```

## Key Components

### Helmholtz Equation Solvers (`src/helmholtz/`)

- `helmholtz_dolfinx_doc_example.py` - DOLFINx-based Helmholtz solver with validation
- `simple_helmholtz.py` - Custom FEM implementation from scratch
- `simple_helmholtz_eig.py` - Eigenvalue problems for Helmholtz equation
- `helmholtz_meshio_dolfinx.py` - Mesh I/O integration with DOLFINx

### Wave Equation Analysis (`src/wave_equation/`)

- `simple_wave.py` - Basic wave equation solver
- `simple_wave_eig.py` - Eigenvalue analysis for wave equations
- `error_wave_analytic_numerical.py` - Error analysis between analytical and numerical solutions
- `wave_analytical/` - Analytical wave solutions
- `wave_3d_animation.gif` - 3D visualization of wave propagation

### Fourier Analysis (`src/fourier_analysis/`)

- `fourier.py` - Fourier series implementations
- `fourier_nh_bc.py` - Fourier methods with non-homogeneous boundary conditions

### FEM Core (`src/fem_core/`)

- `numerical.py` - Core numerical methods
- `loc_matrices.py` - Local matrix assembly routines
- `simple_poisson.py` - Poisson equation solver
- `numerical_nh_*.py` - Non-homogeneous boundary condition implementations

### Grid Generation (`src/grid_generation/`)

- `gen_p1grid.py` - P1 triangular grid generation
- `mesh_create.py` - Mesh creation utilities
- `mesh_reading.py` - Mesh I/O operations

### Validation (`src/validation/`)

- `linton1992*.py` - Validation against Linton (1992) benchmark results
- `ellip.py` - Elliptic function computations for validation

### MATLAB Implementations (`src/matlab_implementations/`)

- `ellip.m` - MATLAB version of elliptic solver
- `gen_p1grid.m` - MATLAB P1 grid generation
- `loc*.m` - Local matrix assembly functions
- `drive_gen_p1grid.m` - Driver script for grid generation

## Dependencies

### Python

- `numpy` - Numerical computations
- `scipy` - Scientific computing
- `matplotlib` - Plotting and visualization
- `dolfinx` - Modern FEniCS for finite element analysis
- `meshio` - Mesh I/O operations
- `mpi4py` - MPI parallelization
- `petsc4py` - PETSc linear algebra

### MATLAB

- Standard MATLAB toolboxes for numerical computation

## Usage

### Running Helmholtz Solvers

```python
# Example: DOLFINx Helmholtz solver
python src/helmholtz/helmholtz_dolfinx_doc_example.py

# Example: Custom FEM implementation
python src/helmholtz/simple_helmholtz.py
```

### Grid Generation

```python
# Generate P1 triangular grids
python src/grid_generation/gen_p1grid.py
```

### Validation Tests

```python
# Run Linton 1992 benchmark
python src/validation/linton1992.py
```

## Mathematical Background

This project implements numerical methods for:

1. **Helmholtz Equation**: -∇²u - k²u = f(x,y)
2. **Wave Equation**: ∂²u/∂t² = c²∇²u
3. **Poisson Equation**: -∇²u = f(x,y)

The implementations use:

- Linear finite elements (P1)
- Triangular mesh discretization
- Various boundary condition types (Dirichlet, Neumann, mixed)
- Error analysis and convergence studies

## Output Files

- Solution files are saved in `solutions/` directory
- Mesh files are stored in `meshes/` directory
- Visualization outputs include XDMF format for ParaView
