# Harness de ejecución del Teorema 2.3

Los scripts científicos y sus criterios se conservan. El harness compartido está en `execution_harness.py`; los tres scripts ahora lo importan. `lattice_sums.py` permanece sin cambios en esta integración.

## Ejecución

```sh
python src/beyn_theorems_validations/second_theorem_X.py
python src/beyn_theorems_validations/second_theorem_Y.py
```

Cada comando crea una carpeta única con fecha de America/Bogota e identificador aleatorio:

- X: `results/theorem_2_3_iii/runs/<fecha>_full_<id>/`
- Y: `results/theorem_2_3_iv/runs/<fecha>_full_<id>/`

`latest_run.json` dentro de cada `runs/` señala la más recientemente iniciada; `run_manifest.json` de esa ejecución registra created/running/completed/interrupted/failed, configuración, fingerprint y versión del harness. Los resultados anteriores quedan intactos en sus carpetas históricas.

## Eficiencia y checkpoints

X ejecuta los epsilons independientes con hasta cuatro procesos (también acotados por las CPUs disponibles) y genera los CSV finales en el orden original. `execution_workers=1` permite ejecución secuencial.

Y conserva el orden secuencial y la continuación mediante `previous_mode`: paralelizarla alteraría el algoritmo. No se usa `execution_workers` para paralelizar ese barrido.

Ambos reutilizan la geometría escalar del kernel mediante una caché acotada de 4096 entradas. La clave incluye posición, epsilon, configuración completa y, para Y, altura a. No se interpola Green ni se cambian ramas, cuadratura o tolerancias.

Los casos guardan logs en `epsilon_cases/<caso>/run.log`. El checkpoint `completed.json` se publica atómicamente sólo cuando la función de validación termina sin excepción. Su fingerprint incluye configuración, argumentos (incluida la semilla anterior en Y), fuentes del script, Green y harness. Esquema, fingerprint y checksum son obligatorios para reutilizarlo. El formato es JSON, sin pickle. Las banderas y diagnósticos se conservan; un caso ambiguo terminado no pasa a certificado por guardarse en caché.

Cada invocación normal crea una ejecución nueva y no reutiliza automáticamente checkpoints de otra carpeta. El trabajo MFS, los diagnósticos de convergencia y la generación de informes conservan su implementación anterior.

## Verificación realizada

- Cinco regresiones nuevas aprobadas: versionado sin sobreescritura; rechazo de caché alterada o incompatible; failure/interruption; semilla de Y en fingerprint; matrices con/sin caché idénticas.
- Matrices complejas M=24 de X y Y comparadas con copias anteriores a este cambio: igualdad exacta.
- Smoke de dos procesos en X: igualdad exacta con ejecución secuencial. Configuración pequeña exclusiva de la prueba, no usada como evidencia del teorema.
- Suite conjunta: 39 tests aprobados.
- No se ejecutaron barridos completos del Teorema 2.3.

Este cambio aporta robustez de ejecución y trazabilidad. No transfiere automáticamente los certificados científicos del Teorema 2.1 a 2.3, ni sustituye la auditoría matemática/numerical específica de X y Y.

## Integración del Teorema 2.1

El primer script delega también la asignación de carpetas, las escrituras atómicas y los checkpoints en el harness. Conserva sus manifiestos científicos, su validación del resultado, nombres baseline_cases y su política de procesos. El fingerprint incluye el fuente del harness; cambios en éste invalidan cachés incompatibles. La suite conjunta de 39 tests pasó después de la integración. Las salidas existentes no se migran ni se modifican.
