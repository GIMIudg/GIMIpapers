# Fases 6 y 7 — calidad de las anotaciones estructurales de metabolitos en GEMs

Agrupa los metabolitos de un GEM según el tipo de problema en su anotación estructural (Fase 6) y
recomienda un InChI por metabolito con su nivel de confianza (Fase 7). Usa como entrada los InChIs que
las Fases 1–5 de `DataMining_GEM.ipynb` obtienen de cinco fuentes: el InChI del GEM, su SMILES,
PubChem, ChEBI y KEGG (y HMDB si hay `.mol` locales).

## Archivos
| Archivo | Qué es |
|---|---|
| `inchi_calidad.py` | Funciones: capas del InChI, neutralización, estereoquímica R/S, tipos de inconsistencia, variables y regla de recomendación. |
| `fase6_7_analisis.py` | El análisis completo, en celdas `# %%`. Corre solo (lee los Excel de la Fase 3) o dentro del notebook (usa `tablas1`). |
| `pruebas_inchi_calidad.py` | 72 casos con InChIs reales y conocidos. |
| `armar_notebook.py` | Genera las celdas de las Fases 6–7 dentro de `DataMining_GEM.ipynb` desde estos `.py`. **El notebook no se edita a mano.** |
| `PLAN.md` | Decisiones de método y por qué. |
| `RESULTADOS.md` | Resultados e interpretación química (corrida del 3–4 nov 2026). |

## Cómo se usa
```bash
pip install rdkit scikit-learn pandas openpyxl matplotlib
python pruebas_inchi_calidad.py                     # debe decir 72 bien, 0 mal
python armar_notebook.py        # actualiza ../notebook/DataMining_GEM.ipynb (las Fases 1–5 no se tocan)
```
Después se corre el notebook completo (Colab o local). Las salidas quedan en `resultados/fase6_7/`.

## Ideas centrales
- Las estructuras se comparan **neutralizadas** (RDKit), para que una diferencia de pH no cuente como otra molécula.
- La estereoquímica se compara por **configuración R/S y E/Z**, no como texto: separa *contradictoria* de *incompleta*.
- Los InChI no estándar con estereo relativa (`/s2`) se tratan como "absoluta no definida".
- MODELO y SMILES salen del mismo GEM y cuentan como un solo voto.

Autor de las Fases 6–7: Juan Carlos Kelly Naranjo (tesis, GIMI – UdeG).
