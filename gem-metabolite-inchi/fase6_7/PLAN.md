# Fases 6 y 7: plan (3-nov-2026)

Pedido del tutor (WhatsApp, 17-oct-2026): analizar los resultados que ya existen, sin meter cosas nuevas.
Fase 6 = análisis general + clustering interpretado químicamente. Fase 7 = tabla final por metabolito
con InChI recomendado y nivel de confianza; lo que no se pueda decidir va a "revisión manual".

## Qué se entrega ahora (sin el notebook del tutor)
- `inchi_calidad.py`: funciones (capas del InChI, tipos de inconsistencia, variables, recomendación).
- `fase6_7_analisis.py`: el análisis completo, en celdas `# %%` (se pasa tal cual a Colab).
- `pruebas_inchi_calidad.py`: batería de casos con InChIs conocidos.
- Salidas en `resultados/fase6_7/`: `tabla_final_fase7.xlsx` y figuras PNG.

## Qué queda para cuando llegue `DataMining_GEM.ipynb`
- Pasar el análisis a ese notebook (el tutor lo pidió ahí).
- Revisar si sus Fases 4–5 calculan variables que aquí no hay y, si es así, sumarlas.

## Datos de entrada
Tabla1_Extraccion de los 4 modelos (los 5 InChIs por metabolito: MODELO, SMILES, PUBCHEM, CHEBI, KEGG).
Se usan los InChIs, no los puntos de la Tabla 2/3: los puntos dicen "coincide con el modelo", no qué difiere.

## Decisiones de método (las tomé yo; marcadas para que el tutor las apruebe)
1. **Unidad de análisis = metabolito sin compartimento** (`atp[c]` y `atp[m]` son uno). Si no, Recon3D
   cuenta la misma estructura varias veces y los grupos se inflan. Se reporta cuántos se fusionaron.
2. **InChIs no estándar del modelo** (`InChI=1/`) se pasan a estándar con RDKit antes de comparar.
3. **Tipos de inconsistencia por capas anidadas**: fórmula → conectividad (`/c`) → hidrógenos (`/h`)
   → estereoquímica (`/b /t /m /s`, isótopos `/i`) → carga/protonación (`/q /p`). Si una capa más fina
   separa fuentes que la anterior juntaba, ese tipo está presente.
4. **Fuentes independientes**: MODELO y SMILES salen del mismo GEM; cuentan como una sola para la confianza.
5. **InChI recomendado** (propuesta):
   - Estructura de consenso ignorando protonación (estructura neutralizada, estereo compatible; ver 8–9).
   - Se entrega en la forma de protonación del GEM si alguna fuente del consenso la tiene.
   - Confianza **alta**: todas las fuentes coinciden y hay ≥2 bases externas.
   - **Media**: mayoría estricta con ≥1 base externa a favor.
   - **Revisión manual**: sin mayoría, solo fuentes del modelo, o una sola fuente.
6. **Prioridad de revisión alta**: el modelo difiere en estructura (no solo en protonación) del consenso
   de ≥2 bases externas. Es el caso típico de anotación equivocada en el GEM.
7. Clustering con variables escaladas; K-means y jerárquico (Ward), k de 2 a 8, se elige con silhouette
   y con que los grupos se puedan interpretar. **Resultado: K-means, k = 8** (ver RESULTADOS.md).

## Decisiones que salieron de revisar casos a mano (3-nov-2026)
8. **Estructura neutralizada** (Uncharger de RDKit) antes de comparar fórmula/conectividad/hidrógenos:
   el InChI estándar no siempre marca la protonación con /p (arsenito AsO3 q-3 vs ácido arsenioso AsH3O3).
9. **Estereoquímica por R/S y E/Z** (etiquetador CIP de RDKit), no leyendo /t /m: separa *contradictoria*
   (configuración opuesta) de *incompleta* (una fuente no define el centro). En el voto, una fuente que no
   define un centro es compatible con las que sí; no puede ganarle a una compatible más específica.
10. **InChI no estándar con /s2 (relativa) o /s3 (racémica)**: se quitan /t /m /s al convertir, porque RDKit
    inventa una configuración absoluta. Sin esto salían 498 contradicciones falsas en Recon3D.
11. **Anillo/cadena**: misma fórmula, distinta conectividad y distinto número de anillos = el mismo azúcar en
    forma abierta y cíclica; se reporta como tipo propio, no como "conectividad".
12. **Variables descartadas del clustering**: n_fuentes (ρ = 1.00 con n_externas), fracción con la fórmula
    del GEM (ρ = 0.91 con la de carga; se usa la de átomos pesados, que no depende del pH), nº de InChIs
    distintos (ρ = −0.94 con el consenso) y nº de capas (describe la molécula, no su anotación).
13. Fuera del clustering, por regla: sin InChI, una sola fuente y fórmula genérica en el GEM (con R).
