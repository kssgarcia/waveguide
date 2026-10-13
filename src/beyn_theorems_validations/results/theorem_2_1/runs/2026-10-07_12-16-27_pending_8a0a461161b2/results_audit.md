# Auditoría de ejecución sólo pendiente

Proceso terminado con código 0 y estado pending-validation-completed. Baseline no repetido.

Conteos: {'ambiguous': 20, 'one': 76, 'error': 3}.

Las filas con unicidad apoyada conservan rango umbral=gap=1, convergencia S0/S1, residual normalizado bajo tolerancia y checks físicos ejecutados/aprobados. Los 20 subcríticos permanecen ambiguos. Las cuatro alturas críticas no se establecen.

Tres casos fallan en el check físico porque la evaluación de Green excede r<=.99d: no se certifican y quedan error. Son una limitación de evaluación fuera del disco de la serie, no evidencia de ausencia o multiplicidad. La siguiente corrección debe adaptar la evaluación del campo en esos puntos (por ejemplo mediante una representación modal validada), conservando los checks y sin ampliar arbitrariamente el dominio de la serie. No se cambió el código durante la ejecución.

- a=0.5, epsilon=0.449474: ValueError: greens_periodic requires r <= 0.99*d; got r=3.9632625115069926, d=4.0.
- a=0.5, epsilon=0.474737: ValueError: greens_periodic requires r <= 0.99*d; got r=3.960111678263245, d=4.0.
- a=0.55, epsilon=0.449474: ValueError: greens_periodic requires r <= 0.99*d; got r=3.9728120360011565, d=4.0.

El informe automático marca baseline not-run porque esta ejecución no importa resultados anteriores; la evidencia baseline sigue en su carpeta original.
