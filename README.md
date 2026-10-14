# Validación numérica de modos atrapados en una guía de ondas

Este repositorio se dedica a la validación numérica del Teorema 2.1 y de los casos 2.3(iii), simetría X, y 2.3(iv), simetría Y, del documento [WG_03.03.25.pdf](docs/references/WG_03.03.25.pdf).

## Organización

| Carpeta | Contenido |
|---|---|
| [src/beyn_theorems_validations](src/beyn_theorems_validations/README.md) | Implementación principal con búsqueda de Beyn, certificados, harness, tests y resultados versionados |
| [src/svd_theorems_validations](src/svd_theorems_validations/README.md) | Implementaciones SVD y estudios previos de convergencia para comparación |
| [docs/references](docs/references/README.md) | PDFs de referencia |
| [docs/audits](docs/README.md) | Evidencia histórica de auditoría |

Los kernels `lattice_sums.py`, el harness y los auxiliares de cada implementación se conservan junto a sus scripts. Sus versiones no se unifican: hacerlo podría modificar la formulación o la identidad de resultados existentes. Las mallas y experimentos antiguos ajenos a estas validaciones se eliminaron; su historial permanece en Git.

## Entorno y pruebas

Python 3.10 o posterior. Desde la raíz:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s src/beyn_theorems_validations/tests
```

`requirements.txt` fija las versiones del entorno registrado en la última validación del Teorema 2.1. Las pruebas pequeñas no ejecutan el barrido científico completo.

## Scripts principales

| Caso | Beyn | SVD previo |
|---|---|---|
| Teorema 2.1 | `first_theorem_beyn.py` | `first_theorem_2.py` |
| Teorema 2.3(iii), X | `second_theorem_X.py` | `second_theorem_X.py` |
| Teorema 2.3(iv), Y | `second_theorem_Y.py` | `second_theorem_Y_1.py`, `second_theorem_Y_2.py` |

Para Beyn, desde la raíz:

```sh
python src/beyn_theorems_validations/first_theorem_beyn.py
python src/beyn_theorems_validations/second_theorem_X.py
python src/beyn_theorems_validations/second_theorem_Y.py
```

Cada comando ejecuta un experimento largo y crea una carpeta nueva en `results/<teorema>/runs/`. Consultar el README de Beyn para ejecuciones dirigidas, pruebas y estado de los manifiestos. Los scripts SVD previos usan rutas relativas: consultar su guía antes de ejecutarlos.

## Estado científico

El Teorema 2.1 tiene una auditoría cerrada con M=24: [informe final](src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b/results_audit.md). La evidencia está distribuida entre ejecuciones originales y validaciones dirigidas; no se presenta como un único barrido completo de la última versión.

El caso X dispone de certificados reforzados y validación del sector impar; ver [criterios y alcance](src/beyn_theorems_validations/docs/theorem_2_3_iii/CORRECTED_VALIDATION.md). La unicidad del problema completo sigue separada del conteo impar. Para Y, consultar su [auditoría y alcance](src/beyn_theorems_validations/docs/theorem_2_3_iv/CORRECTED_VALIDATION.md). Los resultados SVD previos se conservan como comparación, sin atribuirles las correcciones de Beyn. El apoyo numérico a resolución finita no constituye una prueba analítica de los teoremas.

Si se interrumpe una ejecución, usar el [runner de reanudación](src/beyn_theorems_validations/docs/RESUME_RUNS.md): recupera los casos compatibles y completa los pendientes en la misma carpeta.
