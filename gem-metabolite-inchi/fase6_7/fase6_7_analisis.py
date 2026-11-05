# %% [markdown]
# # Fases 6 y 7: calidad de las anotaciones estructurales de metabolitos en GEMs
# Entrada: Tabla1_Extraccion de la Fase 3 (los InChIs de MODELO, SMILES, PubChem, ChEBI y KEGG).
# Salida: `resultados/fase6_7/` con la tabla final (Fase 7), el resumen y las figuras.
# Reglas y decisiones de método: `fase6_7/PLAN.md`.

# %% Configuración
import os, sys
import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
RAIZ = os.path.dirname(AQUI) if os.path.basename(AQUI) == "fase6_7" else AQUI
sys.path.insert(0, os.path.join(RAIZ, "fase6_7"))
from inchi_calidad import FUENTES, TIPOS, base_id, variables, recomendar

# Dentro del notebook DataMining_GEM se usan las tablas de su Fase 3 (`tablas1`) y su carpeta de
# resultados (`OUT_DIR`). Fuera de él, se leen los Excel de la Fase 3 de `resultados/`.
EN_NOTEBOOK_TUTOR = "tablas1" in globals()
if EN_NOTEBOOK_TUTOR:
    SALIDA = os.path.join(OUT_DIR, "fase6_7")
    MODELOS_F67 = list(tablas1)
else:
    ENTRADA = os.path.join(RAIZ, "resultados")
    SALIDA = os.path.join(RAIZ, "resultados", "fase6_7")
    MODELOS_F67 = ["smallModel", "ecoli_core_model", "iDopaNeuroC", "Recon3D_301"]
os.makedirs(SALIDA, exist_ok=True)
SEMILLA = 42

# %% 1. Carga y fusión de compartimentos
filas, conteo = [], []
for m in MODELOS_F67:
    if EN_NOTEBOOK_TUTOR:
        t1 = tablas1[m].copy()
    else:
        t1 = pd.read_excel(os.path.join(ENTRADA, f"{m}_tablas.xlsx"), "Tabla1_Extraccion")
    for f in FUENTES:
        if f not in t1.columns:
            t1[f] = None
    t1["Modelo"] = m
    t1["Metabolito"] = t1["ID"].astype(str).map(base_id)
    comps = t1.groupby("Metabolito")["ID"].apply(lambda s: ", ".join(s.astype(str)))
    u = t1.drop_duplicates("Metabolito").copy()
    u["Compartimentos"] = u["Metabolito"].map(comps)
    conteo.append({"Modelo": m, "Filas_Fase3": len(t1), "Metabolitos_unicos": len(u)})
    filas.append(u)
datos = pd.concat(filas, ignore_index=True)
conteo = pd.DataFrame(conteo)
print(conteo.to_string(index=False))

# %% 2. Variables de consistencia y recomendación por metabolito
def por_fila(r):
    inchis = {f: r[f] for f in FUENTES}
    v = variables(inchis, r["Formula"], r["Charge"])
    rec = recomendar(inchis, r["Charge"])
    v.update({f"rec_{k}": x for k, x in rec.items()})
    return pd.Series(v)

var = datos.apply(por_fila, axis=1)
tabla = pd.concat([datos, var], axis=1)
print(tabla[["n_fuentes", "tipo_principal", "rec_confianza"]].describe(include="all").to_string())

# %% [markdown]
# ## Fase 6.0 — ¿Qué variables de la Fase 4 sirven para agrupar?
# Antes de agrupar se revisó qué miden las variables de `inchi_features`. Comparan cada capa del InChI como
# texto, tal cual llega de cada fuente. Al revisar casos a mano aparecieron cuatro efectos que no son
# diferencias químicas: (1) InChIs no estándar (`InChI=1/`) frente a estándar; (2) estereoquímica relativa
# (`/s2`), que no fija la configuración absoluta; (3) la protonación, que cambia `/p` y a veces la fórmula
# (p. ej. arsenito AsO3³⁻ frente a ácido arsenioso As(OH)3) sin cambiar la molécula; (4) una fuente sin capa
# de estereo cuenta como "sin dato", así que el consenso de estereo no detecta estereoquímica incompleta.
# La celda siguiente cuantifica cuánto pesa cada efecto. Por eso las Fases 6–7 calculan sus variables sobre
# InChIs estándar, estructuras neutralizadas (RDKit) y configuraciones R/S, y conservan las de la Fase 4
# como referencia.

# %% 3a. Cuánto pesan esos efectos en el "InChI completo distinto" de la Fase 4
if "inchi_features" in globals():
    from inchi_calidad import a_estandar, neutro
    _efecto = []
    for _m, _t in tablas1.items():
        for _, _r in _t.iterrows():
            _p = {f: _r[f] for f in FUENTES if f in _t.columns and isinstance(_r[f], str) and _r[f].startswith("InChI=")}
            if len(_p) < 2 or len(set(_p.values())) == 1:
                continue
            _e = {f: a_estandar(v) for f, v in _p.items()}
            _e = {f: v for f, v in _e.items() if v}
            _efecto.append({"Modelo": _m,
                            "solo notacion (estandar vs no estandar)": len(set(_e.values())) <= 1,
                            "solo protonacion": len(set(_e.values())) > 1 and len({neutro(v) for v in _e.values()}) == 1,
                            "estructura realmente distinta": len({neutro(v) for v in _e.values()}) > 1})
    _efecto = pd.DataFrame(_efecto)
    tabla_efectos = _efecto.groupby("Modelo").sum().astype(int)
    tabla_efectos.loc["Total"] = tabla_efectos.sum()
    tabla_efectos.insert(0, "InChI completo distinto (Fase 4)", tabla_efectos.sum(axis=1))
    print("Filas (metabolito × compartimento) con ≥2 fuentes cuyo InChI completo no coincide tal cual:")
    print(tabla_efectos.to_string())

# %% [markdown]
# ## Fase 6.1 — Análisis general
# Se comparan las fuentes capa por capa del InChI. Un tipo de inconsistencia está presente si su
# capa separa fuentes que la capa anterior todavía juntaba (fórmula → conectividad → hidrógenos →
# estereoquímica → carga/protonación). Un metabolito puede tener más de un tipo.

# %% 3. Consistentes vs inconsistentes y tipos más frecuentes
comparables = tabla[tabla["n_fuentes"] >= 2]
gen = []
for m in MODELOS_F67 + ["Total"]:
    t = tabla if m == "Total" else tabla[tabla["Modelo"] == m]
    c = t[t["n_fuentes"] >= 2]
    fila = {"Modelo": m, "Metabolitos": len(t), "Sin InChI": int((t["n_fuentes"] == 0).sum()),
            "Una fuente": int((t["n_fuentes"] == 1).sum()), "Comparables (≥2 fuentes)": len(c),
            "Consistentes (InChI idéntico)": int((c["tipo_principal"] == "consistente").sum()),
            "Con inconsistencias": int((c["tipo_principal"] != "consistente").sum())}
    for tp in TIPOS:
        fila[f"inc_{tp}"] = int(c[f"inc_{tp}"].sum())
    gen.append(fila)
general = pd.DataFrame(gen)
print(general.to_string(index=False))

# %% 4. ¿Qué fuente es la que discrepa? (frente a la estructura ganadora de cada metabolito)
from inchi_calidad import a_estandar, principal, contradice, especificidad, neutro, _ganador
from inchi_calidad import carga_neta
filas = []
for _, r in comparables[comparables["n_externas"] >= 2].iterrows():   # solo donde hay mayoría real
    p = {f: a_estandar(r[f]) for f in FUENTES}
    p = {f: v for f, v in p.items() if v}
    g = p[_ganador(p)[0]]
    for f, v in p.items():
        misma = principal(v) == principal(g)
        filas.append({"fuente": f.replace("INCHI_", ""),
                      "constitucion distinta": not misma,
                      "estereo contradictoria": misma and contradice(v, g),
                      "estereo menos definida": misma and not contradice(v, g) and especificidad(v) < especificidad(g),
                      "carga distinta a la del GEM": carga_neta(v) != r["Charge"]})
discrepa = (pd.DataFrame(filas).groupby("fuente").mean() * 100).round(1)
discrepa.columns = [f"% {c}" for c in discrepa.columns]
print(discrepa.to_string())

# %% 5. Distribución de los scores de consenso
import matplotlib
try:
    get_ipython  # noqa: F821  (existe dentro de Jupyter/Colab)
    EN_JUPYTER = True
except NameError:
    EN_JUPYTER = False
    matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cerrar(fig):
    """En el notebook muestra la figura; corriendo como script solo la cierra (ya quedó guardada)."""
    plt.show() if EN_JUPYTER else plt.close(fig)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for a, col, tit in [(ax[0], "consenso_inchi_completo", "InChI completo"),
                    (ax[1], "consenso_estructura", "Estructura (sin protonación)")]:
    for m in ["ecoli_core_model", "iDopaNeuroC", "Recon3D_301"]:
        a.hist(comparables.loc[comparables["Modelo"] == m, col], bins=np.linspace(0, 1, 11),
               alpha=0.5, label=m, density=True)
    a.set_title(f"Consenso: {tit}"); a.set_xlabel("fracción de fuentes que coinciden con la mayoría")
ax[0].set_ylabel("densidad"); ax[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(SALIDA, "fig1_consenso.png"), dpi=150); cerrar(fig)

fig, ax = plt.subplots(figsize=(7, 4))
(general.set_index("Modelo").loc[MODELOS_F67, [f"inc_{t}" for t in TIPOS]]
 .div(general.set_index("Modelo").loc[MODELOS_F67, "Comparables (≥2 fuentes)"], axis=0) * 100
 ).rename(columns=lambda c: c[4:]).plot.bar(ax=ax)
ax.set_ylabel("% de metabolitos comparables"); ax.set_title("Tipos de inconsistencia por modelo")
fig.tight_layout(); fig.savefig(os.path.join(SALIDA, "fig2_tipos_inconsistencia.png"), dpi=150); cerrar(fig)

# %% 6. Casos representativos
cols_caso = ["Modelo", "Metabolito", "Nombre", "Formula", "Charge", "n_fuentes", "tipos_inconsistencia",
             "consenso_estructura", "tanimoto_minimo", "rec_motivo"]
bien = comparables[(comparables["n_externas"] == 3) & (comparables["consenso_estructura"] == 1)
                   & (comparables["frac_carga_gem"] < 1)].head(5)
perfecto = comparables[(comparables["n_externas"] == 3) & (comparables["tipo_principal"] == "consistente")].head(5)
problema = tabla[tabla["rec_prioridad_revision"] == "alta"].sort_values("tanimoto_minimo").head(10)
print("Consistentes en todo:"); print(perfecto[cols_caso].to_string(index=False))
print("\nMisma estructura, solo cambia la protonación:"); print(bien[cols_caso].to_string(index=False))
print("\nProblemáticos (el GEM difiere del consenso externo):"); print(problema[cols_caso].to_string(index=False))

# %% [markdown]
# ## Fase 6.2 — Clustering no supervisado
# Unidad: metabolito con ≥2 fuentes (con una sola no hay consistencia que medir). Los que aparecen en
# varios modelos con los mismos datos se cuentan una vez. Variables elegidas tras explorar (ver PLAN.md):
# se quitaron las redundantes (|ρ| ≥ 0.9) y las que describen la molécula y no su anotación.

# %% 7. Variables y selección de k
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.decomposition import PCA

VARS = ["n_externas", "consenso_inchi_completo", "consenso_estructura", "consenso_capa_principal",
        "tanimoto_minimo", "frac_carga_gem", "frac_pesados_gem"]
llave = ["Metabolito", "Formula", "Charge"] + FUENTES
C = comparables[~comparables["formula_gem_generica"]].drop_duplicates(subset=llave).reset_index(drop=True)
print(f"Para clustering: {len(C)} metabolitos únicos (de {len(comparables)} filas comparables en {len(MODELOS_F67)} modelos; "
      f"{int(comparables['formula_gem_generica'].sum())} con fórmula genérica en el GEM van aparte)")
X = StandardScaler().fit_transform(C[VARS].astype(float))

sel = []
for k in range(2, 9):
    km = KMeans(n_clusters=k, n_init=20, random_state=SEMILLA).fit(X)
    wa = AgglomerativeClustering(n_clusters=k, linkage="ward").fit(X)
    for nombre, lab in [("K-means", km.labels_), ("Jerárquico (Ward)", wa.labels_)]:
        sel.append({"metodo": nombre, "k": k, "silhouette": silhouette_score(X, lab, random_state=SEMILLA),
                    "davies_bouldin": davies_bouldin_score(X, lab)})
sel = pd.DataFrame(sel)
print(sel.pivot(index="k", columns="metodo", values=["silhouette", "davies_bouldin"]).round(3).to_string())

fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
for a, met in zip(ax, ["silhouette", "davies_bouldin"]):
    for nombre, d in sel.groupby("metodo"):
        a.plot(d["k"], d[met], marker="o", label=nombre)
    a.axvline(8, color="grey", ls="--", lw=0.8); a.set_xlabel("k"); a.set_title(met)
ax[0].legend(); fig.tight_layout(); fig.savefig(os.path.join(SALIDA, "fig3_seleccion_k.png"), dpi=150); cerrar(fig)

# %% [markdown]
# **Elección: K-means con k = 8.** k = 2 tiene el mejor silhouette pero solo separa "hay bases externas" de
# "no hay". Entre k ≥ 3, k = 8 tiene el silhouette más alto y el Davies-Bouldin más bajo, y cada grupo tiene
# una lectura química clara. Con k = 7 se juntan "consistente con bases externas" y "consistente solo con el
# GEM", que no significan lo mismo. K-means supera a Ward en todos los k.

# %% 8. Grupos finales y su interpretación
K = 8
print(sel[(sel["metodo"] == "K-means") & (sel["k"] == K)].round(3).to_string(index=False))
modelo_km = KMeans(n_clusters=K, n_init=20, random_state=SEMILLA).fit(X)
C["grupo_num"] = modelo_km.labels_
perfil = C.groupby("grupo_num")[VARS + [f"inc_{t}" for t in TIPOS] + ["n_capas_media"]].mean()
perfil.insert(0, "n", C["grupo_num"].value_counts().sort_index())

def nombrar(p):
    """Reglas sobre el promedio del grupo (orden = de lo más específico a lo más general)."""
    if p["frac_pesados_gem"] < 0.5:
        return "La fórmula del GEM no corresponde a su estructura"
    if p["consenso_capa_principal"] < 0.9 or p["tanimoto_minimo"] < 0.6:
        return "Discrepancias de constitución (fórmula, conectividad o anillo/cadena)"
    if p["inc_estereo_contradictoria"] > 0.5:
        return "Estereoquímica contradictoria (posible enantiómero equivocado)"
    if p["frac_carga_gem"] < 0.2:
        return "Carga/fórmula del GEM ionizadas, estructura neutra"
    if p["inc_estereo_incompleta"] > 0.5:
        return "Estereoquímica incompleta, poca evidencia externa"
    if p["n_externas"] < 1:
        return "Información insuficiente (solo el GEM)"
    if p["frac_carga_gem"] < 0.7:
        return "Misma estructura, difiere solo la protonación"
    return "Alta consistencia"

perfil.insert(0, "grupo", perfil.apply(nombrar, axis=1))
if not perfil["grupo"].is_unique:   # con otra versión de RDKit podría pasar: se numeran en vez de detenerse
    print("Aviso: dos grupos recibieron el mismo nombre; se numeran. Revisar el perfil antes de interpretar.")
    perfil["grupo"] = [f"{g} ({i + 1})" if (perfil["grupo"] == g).sum() > 1 else g
                       for i, g in enumerate(perfil["grupo"])]
C["grupo"] = C["grupo_num"].map(perfil["grupo"])
print(perfil.round(2).to_string())

coords = PCA(n_components=2, random_state=SEMILLA).fit(X)
P = coords.transform(X)
fig, ax = plt.subplots(figsize=(8, 6))
rng = np.random.default_rng(SEMILLA)
for g, d in C.groupby("grupo"):
    i = d.index.values
    ax.scatter(P[i, 0] + rng.normal(0, 0.05, len(i)), P[i, 1] + rng.normal(0, 0.05, len(i)), s=8, alpha=0.6,
               label=f"{g} (n={len(i)})")
v = coords.explained_variance_ratio_ * 100
ax.set_xlabel(f"PC1 ({v[0]:.0f} %)"); ax.set_ylabel(f"PC2 ({v[1]:.0f} %)")
ax.set_title(f"Grupos de consistencia estructural (K-means, k = {K})")
ax.legend(fontsize=7, loc="best"); fig.tight_layout()
fig.savefig(os.path.join(SALIDA, "fig4_grupos_pca.png"), dpi=150); cerrar(fig)

hm = perfil.set_index("grupo")[VARS]
fig, ax = plt.subplots(figsize=(9, 4))
norm = (hm - hm.min()) / (hm.max() - hm.min()).replace(0, 1)
im = ax.imshow(norm, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(VARS))); ax.set_xticklabels(VARS, rotation=35, ha="right", fontsize=8)
ax.set_yticks(range(len(hm))); ax.set_yticklabels(hm.index, fontsize=8)
for (r, c_), val in np.ndenumerate(hm.values):
    ax.text(c_, r, f"{val:.2f}", ha="center", va="center", fontsize=7,
            color="black" if norm.values[r, c_] > 0.6 else "white")
ax.set_title("Perfil promedio de cada grupo (color escalado por columna)")
fig.tight_layout(); fig.savefig(os.path.join(SALIDA, "fig5_perfil_grupos.png"), dpi=150); cerrar(fig)

# %% [markdown]
# ### Interpretación química de los grupos
# Los nombres salen de reglas fijas sobre el promedio de cada grupo (función `nombrar`), así que la lectura
# es la misma en cada corrida. La tabla de perfiles y la figura 5 muestran los números en que se apoya.
#
# - **Alta consistencia.** Varias bases externas y todas coinciden, incluso en la carga. Suelen ser
#   moléculas pequeñas, neutras y sin grupos ionizables (agua, CO2, etanol). Anotación correcta.
# - **Misma estructura, difiere solo la protonación.** Varias bases externas; al neutralizar, la estructura
#   es idéntica, pero el GEM la guarda a pH fisiológico (fosfatos y carboxilatos desprotonados: NAD, CoA,
#   PEP) y las bases en forma neutra. No es error de anotación: es el estado químico que el modelo necesita
#   para balancear cargas en las reacciones. Explica la baja coincidencia de carga de la Fase 3.
# - **Estereoquímica contradictoria.** Alguna fuente da la configuración opuesta (R frente a S) en al menos
#   un centro. Es el problema más serio: el GEM o una base puede estar apuntando al enantiómero equivocado
#   (p. ej. `lac_D` en ecoli: su InChI es L-lactato).
# - **Discrepancias de constitución.** Fórmula o conectividad distintas: azúcares en forma abierta en una
#   fuente y cíclica en otra (mismo metabolito), poliglutamatos o polímeros representados con otro tamaño, e
#   IDs externos que apuntan a otro compuesto (p. ej. galactosamina frente a galactosamina-6-fosfato).
# - **Carga/fórmula del GEM ionizadas, estructura neutra.** Ninguna fuente, ni el InChI del propio GEM, tiene
#   la carga declarada en el GEM: los campos de fórmula y carga del modelo no corresponden a su estructura.
#   Es una inconsistencia interna del GEM, casi toda de Recon3D.
# - **La fórmula del GEM no corresponde a su estructura.** Los átomos pesados de la fórmula del GEM no son
#   los del InChI (p. ej. una acilcarnitina C20:4 con InChI de C24). Error de anotación en el GEM.
# - **Estereoquímica incompleta, poca evidencia externa.** El InChI del GEM trae estereoquímica relativa
#   (`/s2`) o no define centros que su SMILES sí define, y casi no hay bases externas para confirmar.
# - **Información insuficiente (solo el GEM).** El GEM es consistente consigo mismo, pero no hay bases
#   externas que lo validen; no se puede afirmar que la anotación sea correcta ni incorrecta.

# %% [markdown]
# ## Fase 7 — Tabla final por metabolito

# %% 9. Grupo, estado de la anotación y tabla final
from inchi_calidad import carga_neta, coincide_pesados_gem
tabla = tabla.merge(C[llave + ["grupo"]], on=llave, how="left")
tabla.loc[tabla["n_fuentes"] == 0, "grupo"] = "Sin InChI en ninguna fuente (regla)"
tabla.loc[tabla["n_fuentes"] == 1, "grupo"] = "Una sola fuente (regla)"
tabla.loc[(tabla["n_fuentes"] >= 2) & tabla["formula_gem_generica"], "grupo"] = "Fórmula genérica en el GEM (regla)"
assert tabla["grupo"].notna().all()

def estado(r):
    if r["rec_confianza"] == "sin_datos" or r["rec_motivo"] in ("solo fuentes del modelo", "una sola base externa"):
        return "Sin evidencia externa suficiente"
    if r["rec_confianza"] == "revision_manual":
        return "Inconsistente: fuentes en conflicto"
    if r["rec_gem_coincide"] is False:
        return "Inconsistente: el GEM difiere del consenso externo"
    try:
        carga_ok = carga_neta(r["rec_inchi"]) == int(r["Charge"])
    except (TypeError, ValueError):
        carga_ok = False
    pesados_ok = coincide_pesados_gem(r["rec_inchi"], r["Formula"]) is not False
    if not (carga_ok and pesados_ok):
        return "Inconsistente: carga/fórmula del GEM no coinciden con la estructura"
    return "Correctamente anotado"

tabla["Estado_anotacion"] = tabla.apply(estado, axis=1)
tabla["Fuentes"] = tabla.apply(lambda r: ", ".join(f[6:] for f in FUENTES if isinstance(r[f], str) and r[f].startswith("InChI")), axis=1)

final = tabla.rename(columns={
    "Formula": "Formula_GEM", "Charge": "Carga_GEM", "n_fuentes": "N_fuentes",
    "n_inchis_distintos": "N_InChIs_distintos", "consenso_inchi_completo": "Consenso_InChI_completo",
    "consenso_estructura": "Consenso_estructura", "tanimoto_minimo": "Tanimoto_minimo",
    "tipos_inconsistencia": "Tipo_inconsistencia", "grupo": "Grupo_clustering",
    "rec_inchi": "InChI_recomendado", "rec_confianza": "Confianza", "rec_motivo": "Motivo",
    "rec_fuentes_a_favor": "Fuentes_a_favor", "rec_prioridad_revision": "Prioridad_revision"})
COLS = ["Modelo", "Metabolito", "Nombre", "Compartimentos", "Formula_GEM", "Carga_GEM", "Fuentes", "N_fuentes",
        "N_InChIs_distintos"] + FUENTES + ["Consenso_InChI_completo", "Consenso_estructura", "Tanimoto_minimo",
        "Tipo_inconsistencia", "Grupo_clustering", "Estado_anotacion", "InChI_recomendado", "Confianza",
        "Fuentes_a_favor", "Motivo", "Prioridad_revision"]
if "inchi_features" in globals():   # variables de la Fase 4 del notebook, como referencia
    _f4 = inchi_features[["Modelo", "ID", "SourceSupport", "EstructuraPrincipalConsensus",
                          "InChICompletoConsensus", "StereoConsensus"]].rename(
        columns=lambda c: c if c in ("Modelo", "ID") else f"Fase4_{c}")
    final = final.merge(_f4, on=["Modelo", "ID"], how="left")
    COLS += [c for c in _f4.columns if c.startswith("Fase4_")]
final = final[COLS]
correctos = final[final["Estado_anotacion"] == "Correctamente anotado"]
inconsistentes = final[final["Estado_anotacion"].str.startswith("Inconsistente")]
print(f"Correctamente anotados: {len(correctos)} | Con inconsistencias: {len(inconsistentes)}")
print(inconsistentes["Estado_anotacion"].value_counts().to_string())
prioritarios = final[final["Prioridad_revision"] == "alta"].sort_values(["Tanimoto_minimo", "Modelo"])

resumen_estado = pd.crosstab(final["Modelo"], final["Estado_anotacion"], margins=True, margins_name="Total")
resumen_grupo = pd.crosstab(final["Grupo_clustering"], final["Modelo"], margins=True, margins_name="Total")
resumen_conf = pd.crosstab(final["Modelo"], final["Confianza"], margins=True, margins_name="Total")
print(resumen_estado.to_string()); print(); print(resumen_grupo.to_string()); print(); print(resumen_conf.to_string())

# %% 10. Excel de entrega
LEEME = pd.DataFrame([
    ("Resumen_estado", "Cuántos metabolitos por modelo están bien anotados, inconsistentes o sin evidencia."),
    ("Resumen_grupos", "Cuántos metabolitos de cada modelo caen en cada grupo."),
    ("Resumen_confianza", "Nivel de confianza de la recomendación por modelo."),
    ("Analisis_general", "Fase 6.1: consistentes vs inconsistentes y cuántos tienen cada tipo de inconsistencia."),
    ("Discrepancia_fuente", "Qué fuente se aparta de la estructura ganadora (metabolitos con ≥2 bases externas), en %."),
    ("Seleccion_k", "Silhouette y Davies-Bouldin para K-means y Ward, k = 2..8."),
    ("Perfil_grupos", "Promedio de cada variable por grupo (base de la interpretación química)."),
    ("Tabla_final", "Fase 7: una fila por metabolito (compartimentos fusionados)."),
    ("Correctamente_anotados", "Lista: el GEM coincide con el consenso externo y su carga/fórmula con la estructura."),
    ("Con_inconsistencias", "Lista: todos los que tienen algún Estado_anotacion 'Inconsistente: …'."),
    ("Casos_prioritarios", "Lista: la estructura del GEM difiere del consenso de ≥2 bases externas."),
    ("Efectos_en_Fase4", "Cuántas diferencias de 'InChI completo' de la Fase 4 son notación, protonación o estructura."),
    ("", ""),
    ("Consenso_InChI_completo", "Fracción de fuentes cuyo InChI completo (texto) es igual al más frecuente."),
    ("Consenso_estructura", "Fracción de fuentes compatibles con la estructura ganadora: misma constitución neutralizada y "
                            "sin configuraciones R/S o E/Z opuestas (una fuente que no define un centro sí es compatible)."),
    ("Tanimoto_minimo", "Similitud mínima entre pares de fuentes (Morgan r=2) sobre las estructuras neutralizadas con RDKit."),
    ("Tipo_inconsistencia", "formula, conectividad, anillo_cadena, hidrogenos, estereo_contradictoria, estereo_incompleta, carga. "
                            "Describe a todas las fuentes, aunque solo una minoritaria discrepe."),
    ("Estado_anotacion", "Juzga al GEM: 'Correctamente anotado' si coincide con el consenso externo y su carga/fórmula con la "
                         "estructura, aunque alguna fuente minoritaria discrepe (ver Tipo_inconsistencia)."),
    ("Confianza", "alta: todas compatibles y ≥2 bases externas; media: mayoría con ≥1 base externa; revision_manual: no se puede "
                  "decidir (Estado_anotacion dice si es por conflicto o por falta de bases externas); sin_datos: sin InChI."),
    ("Prioridad_revision", "alta: el GEM difiere del consenso de ≥2 bases externas; media: fuentes en conflicto; baja: el resto."),
    ("Fase4_*", "Variables de la Fase 4 (inchi_features) del primer compartimento, como referencia."),
], columns=["Hoja / columna", "Qué es"])

ruta = os.path.join(SALIDA, "tabla_final_fase7.xlsx")
with pd.ExcelWriter(ruta) as xw:
    LEEME.to_excel(xw, "Leeme", index=False)
    resumen_estado.to_excel(xw, "Resumen_estado")
    resumen_grupo.to_excel(xw, "Resumen_grupos")
    resumen_conf.to_excel(xw, "Resumen_confianza")
    conteo.merge(general, on="Modelo", how="right").to_excel(xw, "Analisis_general", index=False)
    discrepa.to_excel(xw, "Discrepancia_fuente")
    sel.round(3).to_excel(xw, "Seleccion_k", index=False)
    perfil.round(3).to_excel(xw, "Perfil_grupos", index=False)
    final.to_excel(xw, "Tabla_final", index=False)
    correctos.to_excel(xw, "Correctamente_anotados", index=False)
    inconsistentes.to_excel(xw, "Con_inconsistencias", index=False)
    prioritarios.to_excel(xw, "Casos_prioritarios", index=False)
    if "tabla_efectos" in globals():
        tabla_efectos.to_excel(xw, "Efectos_en_Fase4")
print("Guardado:", ruta)
