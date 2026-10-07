# Validaciones de Beyn

| Script científico | Teorema |
|---|---|
| `first_theorem_beyn.py` | 2.1 |
| `second_theorem_X.py` | 2.3 iii |
| `second_theorem_Y.py` | 2.3 iv |

`tests/` contiene regresiones pequeñas, `runners/` los ejecutores auxiliares y `docs/` la documentación por teorema. El Green compartido está en `lattice_sums.py`.

## Teorema 2.1: ejecuciones automáticas

Desde la raíz del repositorio:

```sh
python src/beyn_theorems_validations/first_theorem_beyn.py
```

Cada ejecución crea una carpeta única, sin editar el código:

```text
results/theorem_2_1/
└── runs/
    ├── latest_run.json
    ├── 2026-10-07_11-30-00_full_a1b2c3d4e5f6/
    │   ├── run_manifest.json
    │   ├── baseline_cases/
    │   ├── summary.csv
    │   └── ... figuras y demás resultados
    └── ... siguientes ejecuciones
```

La fecha usa America/Bogota. Un identificador aleatorio evita colisiones incluso si se inicia más de una ejecución en el mismo segundo. `runs/latest_run.json` identifica la ejecución más recientemente iniciada, no necesariamente la más recientemente terminada. Consultar su `run_manifest.json` para el estado: created, running, completed, interrupted o failed (los estados de etapa dan más detalle).

Alternativa con pruebas y representantes antes del barrido:

```sh
python src/beyn_theorems_validations/runners/validate_first_theorem_beyn.py --full-after-gate
```

Este comando crea una sola ejecución `gated`: las pruebas, los representantes y el barrido completo usan la misma carpeta. No hace falta ejecutar los dos comandos.

Los runners de convergencia y benchmark también crean carpetas propias. El de convergencia ejecuta primero su puerta representativa en esa misma carpeta.

**Una nueva invocación comienza una ejecución nueva:** no sobrescribe ni reutiliza salidas de otra ejecución. Los checkpoints con fingerprint/checksum se mantienen dentro de cada ejecución. El baseline conserva hasta cuatro procesos, M=24 y todos los criterios de certificación. Las cachés internas nunca permiten convertir un fallo o diagnóstico ambiguo en una certificación.

Sólo pruebas pequeñas, sin iniciar un experimento:

```sh
python -m unittest discover -s src/beyn_theorems_validations/tests -p test_first_theorem_beyn.py -v
```

## Limpieza realizada

A petición del usuario se eliminaron todos los resultados anteriores del Teorema 2.1, incluidos legacy, versiones corregidas y pilotos. No hay resultados científicos nuevos hasta ejecutar el script. Las referencias a ejecuciones antiguas en documentos históricos no representan archivos disponibles.

Los resultados de los Teoremas 2.3 iii y iv se conservan bajo `results/theorem_2_3_iii/` y `results/theorem_2_3_iv/`. Sus reglas de salida no se modificaron por este versionado.
