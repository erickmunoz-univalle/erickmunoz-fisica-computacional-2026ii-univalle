"""Taller 6 -- Módulo de análisis de datos nucleares en Python."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

CARPETA = Path("taller_nucleo_atomico")


def _resolver_ruta(nombre_archivo):
    p = Path(nombre_archivo)
    if p.is_file():
        return p
    if (CARPETA / nombre_archivo).is_file():
        return CARPETA / nombre_archivo
    return p


def leer_valle(archivo="valle.dat"):
    """Lee el archivo valle.dat escrito por el programa en Fortran.

    Entrada: archivo (str o Path)
    Salida: tupla (A, Z, EB) con tres arreglos de NumPy.
    """
    ruta = _resolver_ruta(archivo)
    datos = np.loadtxt(ruta)
    return datos[:, 0], datos[:, 1], datos[:, 2]


def encontrar_maximo_enlace(a_arr, z_arr, eb_arr):
    """Encuentra el núcleo más ligado por nucleón (mínimo EB) mediante un ciclo.

    Cumple con la restricción del taller: no utiliza min, np.min ni np.argmin.
    """
    a_min = a_arr[0]
    z_min = z_arr[0]
    eb_min = eb_arr[0]

    for a_val, z_val, eb_val in zip(a_arr, z_arr, eb_arr):
        if eb_val < eb_min:
            eb_min = eb_val
            a_min = a_val
            z_min = z_val

    return a_min, z_min, eb_min


def guardar_maximo(a_max, z_max, eb_max, archivo="maximo_python.dat"):
    """Guarda el núcleo más ligado encontrado por Python en un archivo .dat."""
    ruta = CARPETA / archivo if CARPETA.is_dir() else Path(archivo)
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(f"{int(a_max):5d} {int(z_max):5d} {eb_max:16.8f}\n")


def graficar_curva_enlace(a_arr, z_arr, eb_arr, experimentales, archivo_salida="curva_enlace.png"):
    """Genera y guarda la gráfica de B/A vs A con los puntos experimentales."""
    AZUL, ROJO = "#162A54", "#B02A20"
    b_sobre_a_modelo = -eb_arr

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(a_arr, b_sobre_a_modelo, "-", color=AZUL, lw=2.0, label="Modelo (Curva de máximo enlace)")

    exp_a = [item[2] for item in experimentales]
    exp_ba = [item[3] for item in experimentales]
    ax.scatter(exp_a, exp_ba, color=ROJO, zorder=5, s=40, edgecolor="black", linewidth=0.5, label="Datos experimentales (AME2020)")

    ax.set_xlabel(r"Número de nucleones $A$", fontsize=11)
    ax.set_ylabel(r"Energía de enlace por nucleón $B/A$ [MeV]", fontsize=11)
    ax.set_title(r"Fórmula Semiempírica de Masa: Curva de Máximo Enlace vs. Datos", fontsize=12, pad=10)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(frameon=True, facecolor="white", edgecolor="none")

    plt.tight_layout()
    ruta_img = CARPETA / archivo_salida if CARPETA.is_dir() else Path(archivo_salida)
    plt.savefig(ruta_img, dpi=300)
    plt.show()
    plt.close()
