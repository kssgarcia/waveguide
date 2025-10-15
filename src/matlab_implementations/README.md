# MATLAB Implementations

This directory contains MATLAB versions of the numerical algorithms implemented in Python. These serve as cross-validation tools and provide compatibility with MATLAB-based workflows.

## Files

### Grid Generation

- `gen_p1grid.m` - Generate P1 triangular grids with uniform refinement
- `drive_gen_p1grid.m` - Driver script for grid generation

### Local Matrix Assembly

- `locA.m` - Element diffusion matrix for P1 elements
- `locB.m` - Element convection matrix
- `locM.m` - Element mass matrix

### Solvers

- `ellip.m` - Finite element solver for 2nd-order convection-diffusion problems

## Usage

Run the driver scripts to test the implementations:

```matlab
% Generate a refined triangular grid
[ne, np, p, conn, gbc] = gen_p1grid(3);

% Run the elliptic solver
ellip
```

## Cross-Validation

These MATLAB implementations should produce equivalent results to their Python counterparts in the corresponding directories. Use them to verify numerical accuracy and algorithm correctness.
