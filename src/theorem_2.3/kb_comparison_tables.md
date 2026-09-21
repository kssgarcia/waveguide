# Comparación KB: BEM vs Asintótico

## Teorema 2.1 — Modo Discreto

| ε     | kb_asym         | kb_BEM          | Δkb (abs)    | E_k (relativo) |
|-------|-----------------|-----------------|--------------|----------------|
| 0.05  | 1.570768582119 | 1.570769187085 | 6.050×10⁻⁷ | 3.85×10⁻⁷    |
| 0.07  | 1.570689740172 | 1.570694864633 | 5.124×10⁻⁶ | 3.26×10⁻⁶    |
| 0.09  | 1.570505049848 | 1.570529959870 | 2.491×10⁻⁵ | 1.59×10⁻⁵    |
| 0.11  | 1.570146262336 | 1.570232594781 | 8.633×10⁻⁵ | 5.50×10⁻⁵    |

---

## Teorema 2.3(iii) — BIC Simetría X

| ε     | kb_asym         | kb_BEM          | Δkb (abs)    | E_k (relativo) |
|-------|-----------------|-----------------|--------------|----------------|
| 0.05  | 3.140631984100 | 3.140675361277 | 4.338×10⁻⁵ | 1.38×10⁻⁵    |
| 0.07  | 3.137900540385 | 3.138273874883 | 3.733×10⁻⁴ | 1.19×10⁻⁴    |
| 0.09  | 3.131493237939 | 3.133273751158 | 1.781×10⁻³ | 5.68×10⁻⁴    |
| 0.11  | 3.119010674785 | 3.124990855494 | 5.980×10⁻³ | 1.92×10⁻³    |

---

## Teorema 2.3(iv) — BIC Simetría Y (a tuning)

| ε     | kb_asym         | kb_BEM          | Δkb (abs)    | E_k (relativo) |
|-------|-----------------|-----------------|--------------|----------------|
| 0.05  | 3.140622485938 | 3.140666069606 | 4.358×10⁻⁵ | 1.39×10⁻⁵    |
| 0.07  | 3.137864020328 | 3.138239627059 | 3.756×10⁻⁴ | 1.20×10⁻⁴    |
| 0.09  | 3.131393237609 | 3.133186483837 | 1.793×10⁻³ | 5.73×10⁻⁴    |
| 0.11  | 3.118786624544 | 3.124813903881 | 6.027×10⁻³ | 1.93×10⁻³    |

---

## Observaciones

- **E_k en modo discreto** es 2–3 órdenes de magnitud menor que en los BICs, porque el gap al cutoff es mucho más pequeño (ε⁴ vs ε²).
- **E_k crece con ε** en los tres casos, reflejando que la aproximación asintótica pierde precisión al alejarse del régimen ε ≪ 1.
- **BIC X vs BIC Y** tienen errores comparables, como ожидаемо dado que comparten la misma estructura asintótica.
- **E_k en BICs (~10⁻³ a 10⁻⁵)** vs **E_σ en σ (~10⁻² a 10⁻¹)**: el error relativo en kb parece más pequeño porque kb ≈ π ≈ 3.14, pero ambos miden la misma discrepancia.
