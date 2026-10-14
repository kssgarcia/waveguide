# Auditoría de las últimas corridas completas

Se identificaron las ejecuciones usando latest_run.json; se verificaron manifiestos, fingerprints exactos, las diez epsilon configuradas, tablas físicas y quince variantes de sensibilidad. No se repitieron los barridos. Los contrastes independientes X/Y se ejecutaron sobre estas últimas carpetas.

| Algoritmo | Corrida | Casos físicos aprobados | Residual normalizado máximo | BIE off-grid máximo |
|---|---|---:|---:|---:|
| theorem_2_1 | [2026-10-13_17-52-48_full_5605ecd253a8](/Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_1/runs/2026-10-13_17-52-48_full_5605ecd253a8) | 10/10 | 6.707e-10 | 7.112e-06 |
| theorem_2_3_iii | [2026-10-13_21-23-19_full_fe445a6af342](/Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_3_iii/runs/2026-10-13_21-23-19_full_fe445a6af342) | 10/10 | 9.294e-11 | 4.033e-03 |
| theorem_2_3_iv | [2026-10-13_20-16-45_full_63388d81f58d](/Users/kevinsepulveda/Documents/waveguide/src/beyn_theorems_validations/results/theorem_2_3_iv/runs/2026-10-13_20-16-45_full_63388d81f58d) | 10/10 | 1.293e-11 | 2.117e-03 |

## Conclusión y alcance

Las tres corridas están completas, compatibles con el código y sin fallo detectado en los controles físicos configurados. Los residuales BIE no son cero: X alcanza 0.004033 y Y 0.002117, ambos bajo el umbral 0.005; se declaran como error de reconstrucción medido a M=24. El contraste modal independiente confirma los operadores en .02/.10/.20 con diferencias inferiores a 5e-10.

2.1 incluye además 99 geometrías (79 modos físicamente certificados y 20 exclusiones numéricas) y cuatro brackets críticos. Su auditoría detallada verifica los certificados independientes para los candidatos SVD sospechosos.

X: existencia física 10/10 y unicidad en el sector impar dentro de los contornos estudiados. No está establecida la unicidad del problema completo. Y: existencia física 10/10; unicidad global no establecida. La verificación independiente de Y pasa también colocación desde a=0, H20/H40, cuadraturas alternativas, sobremuestreo y secciones desplazadas.

Y conserva como rechazadas tres variantes de baja precisión: paso 1e-5, paso 1e-4 y harmonic_order=5. Esto documenta sensibilidad, no un error de los casos principales ni inexistencia física.

## Figuras

Se inspeccionaron las nueve figuras principales de X y dieciocho de Y. Las 99 figuras principales de la misma corrida 2.1 fueron revisadas en su results_audit.md. Los gráficos son coherentes con las tablas y con las identidades de dispersión. Las curvas asintóticas concuerdan mejor al disminuir epsilon; su separación a epsilon mayor no invalida el término principal. Los screens de Y muestran el mínimo refinado aunque la malla de dibujo no lo resuelva. Los puntos únicos a M=24 no representan convergencia de malla.

La discrepancia interna del ejemplo 5.4 (-beta/12 impreso frente a -beta/2 al sustituir en 5.47) sigue documentada; los resultados de Y se comparan con la formulación (2.7). No se cambió la física para ajustar ese ejemplo.

Estas son validaciones numéricas a resolución fija y dominios finitos. No prueban uniformidad del Big-O, analiticidad ni exclusión arbitrariamente cerca de los umbrales.

Los hashes de las tablas y fingerprints quedan en verification.json.
