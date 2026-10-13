# Validaciones SVD previas

Esta carpeta conserva los scripts y resultados anteriores para comparación con la implementación principal de [Beyn](../beyn_theorems_validations/README.md). No incorpora automáticamente las correcciones ni los certificados de la auditoría de Beyn.

| Archivo | Función |
|---|---|
| `first_theorem_2.py` | Teorema 2.1 |
| `first_theorem_2_lattice_convergence.py` | Estudio de términos lattice |
| `first_theorem_2_internal_parameter_convergence.py` | Estudio de parámetros internos |
| `second_theorem_X.py` | Teorema 2.3(iii), simetría X |
| `second_theorem_Y_1.py` | Teorema 2.3(iv), primera variante |
| `second_theorem_Y_2.py` | Variante Y con ajuste de la posición del obstáculo |
| `lattice_sums.py`, `dipole_theorem.py`, `gauss_legendre.py` | Dependencias locales necesarias |

Los logs, figuras, tablas y carpetas `theorem_*` son evidencia histórica conservada. No se movieron ni se modificaron para no romper referencias o confundirlos con nuevas ejecuciones.

Estos scripts usan directorios de salida relativos al directorio de trabajo y no el harness versionado de Beyn. Para ejecutarlos, situarse aquí:

```sh
cd src/svd_theorems_validations
python first_theorem_2.py
```

Pueden escribir sobre sus directorios configurados; para la continuación auditada del proyecto se usa la implementación Beyn. No se ejecutó ningún experimento SVD durante la limpieza.
