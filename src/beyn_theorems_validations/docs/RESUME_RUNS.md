# Reanudar después de una interrupción

Desde `src/beyn_theorems_validations/`:

```sh
python runners/resume_validation.py --theorem 2.1
python runners/resume_validation.py --theorem X
python runners/resume_validation.py --theorem Y
```

Por defecto se selecciona la última carpeta de ese teorema. También se puede indicar una ejecución concreta:

```sh
python runners/resume_validation.py --theorem 2.1 --run results/theorem_2_1/runs/2026-10-13_17-52-48_full_5605ecd253a8
```

No cambies código/configuración durante una ejecución. La reanudación exige que el fingerprint coincida con la implementación y parámetros guardados; si cambian, rechaza la reutilización y hay que iniciar una ejecución nueva. Ejecutar nuevamente el script principal crea una carpeta nueva; el runner de reanudación continúa la carpeta existente.

2.1 recupera checkpoints del baseline y certificados atómicos del barrido de geometrías, registra su adopción con hashes y conserva estudios completos de sensibilidad/alturas críticas. Los casos pendientes se ejecutan y guardan individualmente. X/Y recuperan checkpoints firmados por epsilon; el caso incompleto vuelve a empezar. Los gráficos e informes se reconstruyen al completar los datos.

Los estados ambiguos siguen siendo ambiguos: checkpoint terminado significa que la evaluación terminó, no que el teorema esté certificado. Los errores no se reutilizan como casos válidos. No se reduce ninguna tolerancia ni se modifica la física.

El avance y los logs quedan en `resume_history/<fecha>/`. Ctrl+C deja el manifiesto como `interrupted`; puede volver a usarse el mismo comando. No se deben ejecutar dos procesos sobre la misma carpeta simultáneamente.
