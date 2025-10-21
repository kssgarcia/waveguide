# Modos Atrapados para un Obstáculo Parametrizado en una Guía de Ondas Bidimensional

**Autor:** Kevin Santiago Sepúlveda García

**Fecha:** 20 de octubre de 2025

## Resumen

Este repositorio contiene la implementación numérica para la validación del Teorema 2.1 de [1], que predice la existencia de modos atrapados en una guía de ondas bidimensional con un obstáculo pequeño.

Utilizando el Método de Elementos de Frontera (BEM), calculamos los valores propios de la ecuación de Helmholtz para diferentes posiciones del obstáculo y parámetros de perturbación ε. El objetivo principal es determinar el valor máximo de ε, denotado ε_máx, para el cual el teorema se cumple. Además, analizamos cómo ε_máx depende de la altura _a_ del obstáculo.

## Estructura del Repositorio

### Scripts de Validación

- **`validation.py`** : Validación principal del teorema usando la transformación logarítmica _s = -2log₁₀(σ)_
- **`validation_kb.py`** : Validación del teorema utilizando el parámetro _kb_
- **`validation_a_vs_eps.py`** : Análisis de la relación entre la altura del obstáculo _a_ y el parámetro de perturbación _ε_

### Utilidades

- **`lattice_sums.py`** : Implementación de funciones auxiliares para el cálculo de la función de Green en la guía de ondas

### Datos

- **`kevin.csv`** : Datos de resultados numéricos
- **`kevin_kb.csv`** : Datos de resultados usando el parámetro _kb_

## Uso

Los scripts pueden ejecutarse directamente con Python:

bash

```bash
python validation.py
python validation_kb.py
python validation_a_vs_eps.py
```
