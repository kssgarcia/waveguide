# Validaciones de Beyn

Los scripts científicos están en la raíz. Las pruebas, los runners, la documentación y los resultados se guardan por separado.

```text
beyn_theorems_validations/
├── first_theorem_beyn.py       # Teorema 2.1
├── second_theorem_X.py         # Teorema 2.3 iii
├── second_theorem_Y.py         # Teorema 2.3 iv
├── lattice_sums.py             # Green compartido
├── tests/                     # Regresiones pequeñas
├── runners/                   # Ejecución por etapas y estudio dirigido
├── docs/
│   ├── theorem_2_1/            # Implementación actual e historial
│   ├── theorem_2_3_iii/
│   └── theorem_2_3_iv/
└── results/
    ├── theorem_2_1/
    │   ├── legacy/pre_fix/     # Salidas anteriores a las correcciones
    │   ├── corrected/v8/      # Barrido completo anterior
    │   ├── corrected/v9/      # Estudio dirigido con nueva cuadratura
    │   └── development/       # Intento interrumpido y piloto
    ├── theorem_2_3_iii/embedded_validation/
    ├── theorem_2_3_iv/
    │   ├── embedded_validation/
    │   └── embedded_validation_v2/
    └── relocation_map.json    # Equivalencia entre rutas antiguas y nuevas
```

## Qué ejecutar

Desde la raíz del repositorio, con el Python del entorno científico:

```sh
# Validación completa del Teorema 2.1:
python src/beyn_theorems_validations/first_theorem_beyn.py

# Alternativa: primero pruebas y representantes; barrido sólo si pasan:
python src/beyn_theorems_validations/runners/validate_first_theorem_beyn.py --full-after-gate

# Sólo las regresiones pequeñas:
python -m unittest discover -s src/beyn_theorems_validations/tests -p test_first_theorem_beyn.py -v
```

Los dos primeros comandos son alternativas: no hace falta ejecutar ambos.

El estudio dirigido se ejecuta con `runners/validate_beyn_convergence.py` después de una puerta representativa aprobada con el fingerprint vigente. Los scripts `second_theorem_X.py` y `second_theorem_Y.py` ejecutan sus teoremas por separado. Todos guardan sus salidas en rutas absolutas bajo `results/`, independientemente del directorio desde el cual se lancen.

## Estado de la evidencia del Teorema 2.1

| Carpeta | Alcance | Estado científico |
|---|---|---|
| `legacy/pre_fix` | Resultados anteriores | Conservados; no reutilizarlos como evidencia corregida |
| `corrected/v8` | 10 baseline, 99 geometrías y altura crítica | Existencia local apoyada; conteo global ambiguo |
| `corrected/v9` | 10 baseline y cuatro subcríticos | Unicidad numérica apoyada en los 10 baseline; ausencia subcrítica ambigua |
| `development` | Pilotos e intento interrumpido | Historial de desarrollo; no sustituye la evidencia final |

V9 aún no contiene un nuevo barrido completo de 99 geometrías ni un nuevo estudio completo de altura crítica. Consultar `corrected/v9/convergence_validation_report.md` y `docs/theorem_2_1/THEOREM_2_1_CORRECTED.md`.

Los resultados históricos se trasladaron sin modificar sus bytes. Por trazabilidad, sus manifiestos, logs e informes mantienen las rutas y fingerprints originales; `results/relocation_map.json` permite localizar esas rutas. La reorganización cambia las rutas del código y el fingerprint actual: las cachés anteriores no deben reutilizarse silenciosamente y el runner dirigido exigirá una nueva puerta vigente. No se recalcularon resultados científicos durante esta reorganización.
