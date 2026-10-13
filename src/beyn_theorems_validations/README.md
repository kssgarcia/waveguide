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

A petición del usuario se eliminaron todos los resultados anteriores del Teorema 2.1, incluidos legacy, versiones corregidas y pilotos. Las ejecuciones posteriores a esa limpieza se conservan en `runs/`, con sus manifiestos e informes de auditoría. Las referencias a ejecuciones antiguas en documentos históricos no representan archivos disponibles.

Los resultados de los Teoremas 2.3 iii y iv se conservan bajo `results/theorem_2_3_iii/` y `results/theorem_2_3_iv/`. Cada nueva ejecución de X/Y se guarda automáticamente en su propia carpeta `runs/<fecha>_full_<id>/`; las salidas históricas se conservan. Ver [harness del Teorema 2.3](docs/THEOREM_2_3_HARNESS.md).

La [validación corregida de Y](docs/theorem_2_3_iv/CORRECTED_VALIDATION.md) distingue existencia física, screen muestreado y unicidad completa. Su runner permite `--stage smoke` (.02/.10/.20 con contraste modal independiente) y `--stage full`; el script principal conserva la ejecución completa y el versionado automático.

## Pruebas del harness de 2.3

```sh
python -m unittest discover -s src/beyn_theorems_validations/tests -p test_second_theorem_harness.py -v
```

X paraleliza los epsilons independientes; Y conserva la continuación secuencial. Ambos tienen manifiestos, checkpoints JSON con fingerprint/checksum, logs por caso y caché acotada de geometría. Sus cálculos y criterios científicos siguen siendo los propios de cada script.

Los tres scripts usan `execution_harness.py` para versionado y checkpoints. No se ejecuta directamente. El primero conserva sus manifiestos y criterios científicos específicos; los fingerprints incluyen el harness compartido.

## Ejecutar sólo lo pendiente del Teorema 2.1

```sh
python src/beyn_theorems_validations/first_theorem_beyn.py --only-pending
```

Ejecuta la rejilla de geometrías (incluidos los nuevos checks físicos y diagnósticos) y la altura crítica. Omite los diez baseline y el estudio interno de convergencia. Crea una carpeta nueva `<fecha>_pending_<id>/`, registra su alcance en `execution_scope.json` y termina con estado `pending-validation-completed`. No importa ni modifica los resultados anteriores: el informe de esta ejecución contiene sólo la evidencia recalculada, y la evidencia baseline anterior debe consultarse en su carpeta original. No requiere editar el código. El comando sin opciones sigue ejecutando la validación completa.

## Repetir sólo geometrías que terminaron con error

```sh
python src/beyn_theorems_validations/runners/retry_first_theorem_errors.py RUTA_A_LA_EJECUCION
```

Lee el CSV de esa ejecución, selecciona únicamente filas `error` y las calcula de nuevo con la implementación actual y la configuración científica guardada. Crea otra carpeta `retry-errors`, registra procedencia y checksum, conserva diagnósticos completos y no modifica la ejecución fuente. No repite baseline, altura crítica ni las otras geometrías. Un resultado ambiguo permanece ambiguo.

Desde Green/assembly v12 la reconstrucción a separación horizontal `|X| >= b` usa la expansión modal del mismo Green de Dirichlet (Linton 4.9), con control de convergencia. Así cubre los puntos lejanos que quedan fuera del disco de convergencia de la representación lattice; el ensamblaje BEM y sus parámetros permanecen iguales.

Para recalcular sólo los casos ambiguos y repetir la altura crítica:

```sh
python src/beyn_theorems_validations/runners/retry_first_theorem_errors.py RUTA_A_LA_EJECUCION --status ambiguous --critical-height
```

El runner admite ambos estados, pese a conservar su nombre original. Usa la configuración guardada, conserva los resultados fuente y genera una nueva carpeta `resolve-pending`. Que termine la ejecución no garantiza que todos los resultados se puedan certificar: revisar los estados y diagnósticos de la carpeta nueva.

## Auditoría terminada del Teorema 2.1

La ejecución `2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b` terminó y se verificó contra el fingerprint de la implementación actual v15: 20 casos antes ambiguos con apoyo numérico de no existencia, cuatro alturas críticas acotadas (ancho 0.0009375), tres representantes recertificados local/físicamente y 54 tests aprobados. M=24 y criterios de modos conservados.

Consultar el [informe final](results/theorem_2_1/runs/2026-10-11_07-34-49_resolve-pending_e921a3fdfe2b/results_audit.md). Es evidencia numérica a resolución finita; no se repitió el barrido completo ni se trasladaron certificados antiguos como si fueran v15. Las carpetas parciales marcadas `stopped-for-additional-audit` no son la ejecución final.

## Puerta científica del caso X (Teorema 2.3 iii)

```sh
python src/beyn_theorems_validations/runners/validate_second_theorem_X.py
```

Valida .02,.10,.20 y continúa el barrido sólo si pasan los certificados físicos y de unicidad del sector impar. Comparte los checkpoints firmados dentro de la nueva ejecución. Ver [criterios corregidos](docs/theorem_2_3_iii/CORRECTED_VALIDATION.md). El flag de unicidad del problema completo se mantiene separado: no se ha excluido el sector par ni los márgenes omitidos.

La ejecución X corregida v3 terminó con 10/10 modos físicos y conteos convergentes del sector impar. [Informe verificado](results/theorem_2_3_iii/runs/2026-10-13_14-59-29_full_0a01e51d6f48/results_audit.md). Las salidas anteriores permanecen intactas; la unicidad del problema completo continúa fuera del certificado impar.
