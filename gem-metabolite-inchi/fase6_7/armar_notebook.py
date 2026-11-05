"""Arma fase6_7_analisis.ipynb (para Colab) a partir de inchi_calidad.py + fase6_7_analisis.py.
Así hay una sola fuente de verdad: si cambia el .py, se vuelve a correr esto.
Uso:
  python armar_notebook.py                                   # solo la versión independiente
  python armar_notebook.py --original RUTA/DataMining_GEM.ipynb --salida RUTA/DataMining_GEM.ipynb
      # además, agrega las Fases 6 y 7 al final del notebook del tutor (sus celdas no se modifican)
"""
import argparse, json, os, re

AQUI = os.path.dirname(os.path.abspath(__file__))


def celda(tipo, texto):
    lineas = texto.strip("\n").split("\n")
    fuente = [l + "\n" for l in lineas[:-1]] + [lineas[-1]]
    c = {"cell_type": tipo, "metadata": {}, "source": fuente}
    if tipo == "code":
        c.update({"execution_count": None, "outputs": []})
    return c


def celdas_de_script(texto):
    out = []
    for bloque in re.split(r"^# %%", texto, flags=re.M)[1:]:
        cabecera, _, cuerpo = bloque.partition("\n")
        if cabecera.strip().startswith("[markdown]"):
            md = "\n".join(re.sub(r"^# ?", "", l) for l in cuerpo.strip("\n").split("\n"))
            out.append(celda("markdown", md))
        else:
            cuerpo = "\n".join(l for l in cuerpo.split("\n")
                               if not re.match(r"\s*from inchi_calidad import", l)
                               and "sys.path.insert" not in l)
            titulo = cabecera.strip()
            if titulo:
                out.append(celda("markdown", f"### {titulo}"))
            if cuerpo.strip():
                out.append(celda("code", cuerpo))
    return out


modulo = open(os.path.join(AQUI, "inchi_calidad.py"), encoding="utf-8").read()
script = open(os.path.join(AQUI, "fase6_7_analisis.py"), encoding="utf-8").read()

celdas = [
    celda("markdown", """# Fases 6 y 7: calidad de las anotaciones estructurales de metabolitos en GEMs

**Cómo correrlo en Colab**
1. Sube a una carpeta `resultados/` (junto al notebook) los 4 archivos de la Fase 3:
   `smallModel_tablas.xlsx`, `ecoli_core_model_tablas.xlsx`, `iDopaNeuroC_tablas.xlsx`, `Recon3D_301_tablas.xlsx`.
2. Ejecuta todo (Entorno de ejecución → Ejecutar todas). Tarda menos de un minuto.
3. Las salidas quedan en `resultados/fase6_7/`.

Las reglas y decisiones de método están explicadas en las celdas de texto y en `PLAN.md`."""),
    celda("code", "# En Colab, descomenta la siguiente línea:\n# !pip install -q rdkit scikit-learn openpyxl"),
    celda("markdown", "## Funciones (módulo `inchi_calidad`)\nCapas del InChI, neutralización, estereoquímica R/S, tipos de inconsistencia, variables y recomendación."),
    celda("code", modulo),
] + celdas_de_script(script)

nb = {"cells": celdas, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                     "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
ruta = os.path.join(AQUI, "fase6_7_analisis.ipynb")
json.dump(nb, open(ruta, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"{ruta}: {len(celdas)} celdas ({sum(c['cell_type'] == 'code' for c in celdas)} de código)")

# ------------------------------------------------------------------ notebook del tutor
# Toma DataMining_GEM.original.ipynb (sus Fases 1-5, sin tocar) y reemplaza su última celda ("Qué sigue")
# por las Fases 6 y 7.
_ap = argparse.ArgumentParser()
# Por omisión, la estructura del repo: gem-metabolite-inchi/{fase6_7, notebook}/
_ap.add_argument("--original", default=os.path.join(os.path.dirname(AQUI), "notebook", "DataMining_GEM.ipynb"))
_ap.add_argument("--salida", default=os.path.join(os.path.dirname(AQUI), "notebook", "DataMining_GEM.ipynb"))
_args = _ap.parse_args()
ORIGINAL = _args.original
if os.path.exists(ORIGINAL):
    tutor = json.load(open(ORIGINAL, encoding="utf-8"))
    textos = ["".join(c["source"]) for c in tutor["cells"]]
    # Si el notebook ya trae las Fases 6-7 (se está regenerando), se corta desde su encabezado;
    # si es el original del tutor, se reemplaza su última celda ("Qué sigue").
    corte = next((i for i, t in enumerate(textos) if t.startswith("## Fases 6 y 7")), None)
    if corte is None:
        assert textos[-1].startswith("## Qué sigue"), "La última celda del notebook del tutor cambió: revisar antes de reemplazar"
        corte = len(textos) - 1
    tutor["cells"] = tutor["cells"][:corte] + [None]   # el None se reemplaza abajo
    nuevas = [
        celda("markdown", """## Fases 6 y 7 — Agrupar metabolitos por tipo de problema y recomendar un InChI

Esta parte usa como insumo las tablas de la Fase 3 (`tablas1`, los InChIs por fuente) y revisa las variables de
la Fase 4 (`inchi_features`). Todo lo que sigue corre en menos de un minuto y guarda sus salidas en
`resultados/fase6_7/`: la tabla final de la Fase 7 en Excel y las figuras.

Primero van las funciones (capas del InChI, neutralización de la estructura, estereoquímica R/S, tipos de
inconsistencia, variables y regla de recomendación). Después, el análisis paso a paso. Requiere RDKit y
scikit-learn (Colab ya trae scikit-learn)."""),
        celda("code", modulo),
    ] + celdas_de_script(script)[1:]   # sin el título de la versión independiente
    tutor["cells"] = tutor["cells"][:-1] + nuevas
    assert len(tutor["cells"]) > len(nuevas)
    ruta_t = _args.salida
    json.dump(tutor, open(ruta_t, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"{ruta_t}: {len(tutor['cells'])} celdas (las primeras 31 son las del tutor, sin cambios)")
