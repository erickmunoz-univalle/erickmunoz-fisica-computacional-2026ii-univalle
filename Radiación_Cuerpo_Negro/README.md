# Radiación de Cuerpo Negro

En esta carpeta se encuentran los códigos, archivos y notebooks trabajados en **Google Colab** para el Taller 8 del curso de Física Computacional (Univalle 2026-2).

### Contenido
- `constante_wien.dat` → Valores finales de las raíces x₅, x₃ y la constante de desplazamiento de Wien b.
- `diagrama de flujo.png ` → Diagrama lógico del algoritmo de búsqueda de ceros y cálculo del máximo espectral.
- `iteraciones_wien.dat` → Historial iteración a iteración de los métodos de bisección, Newton y secante para n = 5.
- `maximo_espectro.dat` → Errores absolutos del pico espectral solar en las cuatro mallas evaluadas.
- `modulo_ceros.f90` → Módulo en Fortran con la función de Wien, su derivada y los métodos de búsqueda de raíces.
- `modulo_romberg.f90` → Módulo en Fortran con integración de Romberg y los integrandos de Stefan-Boltzmann.
- `programa_sigma.f90` → Programa en Fortran para el cálculo numérico de la integral de Planck y la constante σ.
- `resultado_sigma.dat` → Resultados y conteo de evaluaciones de Romberg con y sin cambio de variable.
- `resultado_wien.dat` → Resumen de raíces e iteraciones obtenidas para n = 5, n = 3 y el caso inicial x₀ = 0.5.
- `taller8_colab.ipynb` → Notebook principal con el entorno de ejecución, pruebas y análisis en Colab.

Estos archivos muestran el proceso de programación, pruebas y resultados obtenidos durante el taller.
