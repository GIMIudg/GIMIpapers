import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')           # ✅ CORRECCIÓN: necesario en entornos sin pantalla
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from collections import defaultdict
import warnings
import logging
import re
import sys

# Silenciar warnings de fuente (inofensivos)
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# ============================================================
# PLOS ONE — ESPECIFICACIONES GLOBALES
# ============================================================
# Ancho:       789–2250 px @ 300 dpi  →  2.63–7.5 pulgadas
# Alto máx.:   2625 px @ 300 dpi      →  8.75 pulgadas
# Resolución:  300 dpi
# Formato:     TIFF + compresión LZW
# Fuente:      Arial 8–12 pt  (fallback: DejaVu Sans si Arial no está instalada)
# Color:       RGB
# ============================================================

DPI        = 300
FONT_SIZE  = 10    # ticks, leyendas, anotaciones
FONT_TITLE = 12    # títulos  (máx. 12 pt según PLOS)
FONT_AXIS  = 11    # etiquetas de ejes
W_FULL     = 7.5   # página completa  (2250 px / 300 dpi)
W_HALF     = 5.2   # columna de texto (1560 px / 300 dpi)
H_MAX      = 8.75  # alto máximo

# ── Detección automática de Arial ──────────────────────────
_available  = {f.name for f in fm.fontManager.ttflist}
FONT_FAMILY = 'Arial' if 'Arial' in _available else 'DejaVu Sans'
print(f"  Fuente: {FONT_FAMILY}"
      + ("  ✔" if FONT_FAMILY == 'Arial'
         else "  ⚠ Arial no encontrada — usando DejaVu Sans (equivalente visual)"))

plt.rcParams.update({
    'font.family':     FONT_FAMILY,
    'font.size':       FONT_SIZE,
    'axes.titlesize':  FONT_TITLE,
    'axes.labelsize':  FONT_AXIS,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE,
    'legend.fontsize': FONT_SIZE,
    'figure.dpi':      DPI,
    'savefig.dpi':     DPI,
})

def save_plos(fig, filename):
    """Guarda como TIFF LZW, 300 dpi, fondo blanco — cumple PLOS One."""
    fig.savefig(
        filename,
        dpi=DPI,
        format='tiff',
        bbox_inches='tight',
        pil_kwargs={'compression': 'tiff_lzw'},
        facecolor='white',
        edgecolor='none',
    )
    size_mb = __import__('os').path.getsize(filename) / 1_048_576
    print(f"  ✅  {filename}  ({size_mb:.2f} MB)")
    plt.close(fig)

# ============================================================
# 1️⃣ CONFIGURACIÓN Y CARGA DE DATOS
# ============================================================

ID_COLUMN         = 'ModelName'
MIN_VALID_SAMPLES = 50

path1 = "/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis/src/Clustering and data analysis PYTHON/Clinical data analysis/ML models using clinical data/Results_clustering_UMAP_seleccion_reducida_conmerge-DATOSNUEVOS/pacientes_clusterizados_todos_sinfiltro.csv"
path2 = "/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis/src/Clustering and data analysis PYTHON/Metabolic data analysis/ML models using metabolic data/resultados_TumorPhenotype_PCA_metrics_actualizado_sinl2/PatientClusters_TumorPhenotype_PCA.csv"

def load_and_standardize_clusters(file_path, suffix, id_col='ModelName'):
    """Load, clean TCGA IDs and rename cluster columns."""
    try:
        df = pd.read_csv(file_path, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"Failed to read {file_path}. Error: {e}")

    df.columns = df.columns.str.strip()
    if id_col not in df.columns:
        raise ValueError(f"ID column '{id_col}' not found in: {file_path}")

    df[id_col] = (df[id_col].astype(str)
                  .str.split('_').str[0]
                  .str.split('.').str[0]
                  .str.slice(0, 16))

    cluster_cols    = [col for col in df.columns if col != id_col]
    rename_mapping  = {col: f"{col}_{suffix}" for col in cluster_cols}
    df = df.rename(columns=rename_mapping)

    return df[[id_col] + list(rename_mapping.values())].copy()

try:
    df1       = load_and_standardize_clusters(path1, suffix='C', id_col=ID_COLUMN)
    df2       = load_and_standardize_clusters(path2, suffix='M', id_col=ID_COLUMN)
    df_merged = df1.merge(df2, on=ID_COLUMN, how='inner')
    print(f"✅ Merge exitoso. Pacientes comunes: {len(df_merged)}")
except Exception as e:
    print(f"❌ ERROR: {e}")
    sys.exit()

cluster_cols_clinico     = [col for col in df_merged.columns if col.endswith('_C')]
cluster_cols_metabolitos = [col for col in df_merged.columns if col.endswith('_M')]

# ============================================================
# 2️⃣ CÁLCULO DE MÉTRICAS (ARI / AMI)
# ============================================================

def calculate_metrics_fixed(df, col1, col2):
    s1        = df[col1].replace(-1, np.nan)
    s2        = df[col2].replace(-1, np.nan)
    valid_idx = s1.dropna().index.intersection(s2.dropna().index)

    if len(valid_idx) < MIN_VALID_SAMPLES:
        return None

    labels1 = s1.loc[valid_idx].astype(int)
    labels2 = s2.loc[valid_idx].astype(int)

    if labels1.nunique() < 2 or labels2.nunique() < 2:
        return None

    return {
        'Clinical_Cluster':  col1,
        'Metabolic_Cluster': col2,
        'ARI': adjusted_rand_score(labels1, labels2),
        'AMI': adjusted_mutual_info_score(labels1, labels2),
        'N_Samples': len(valid_idx),
    }

results_list = []
for c_col in cluster_cols_clinico:
    for m_col in cluster_cols_metabolitos:
        res = calculate_metrics_fixed(df_merged, c_col, m_col)
        if res:
            results_list.append(res)

df_results = pd.DataFrame(results_list)

# ============================================================
# 3️⃣ TOP 10 — HEATMAP  (Fig 1)
# ============================================================

def extract_silhouette_clinico(col_name):
    """Extrae el Silhouette Score del nombre de columna clínica.
    Patrón: ..._K{k}_(S{score})_DB...
    Ejemplo: Cluster_KMeans_UMAP_C3_NN10_MD0.05_abc_S123_K2_S0.712_DB0.45_CH210_C
                                                                     ^^^^^^
    """
    m = re.search(r'_K\d+_(S[\d.]+)_DB', col_name)
    return float(m.group(1).replace('S', '')) if m else None

def extract_silhouette_metabolico(col_name):
    """Extrae el Silhouette Score del nombre de columna metabólica.
    Patrón: ..._K{k}_(S{score})_DB...
    Ejemplo: Cluster_KMeans_C3_W5_PCA_K3_S0.701_DB0.79_CH245_Seed42_M
                                          ^^^^^^
    """
    m = re.search(r'_K\d+_(S[\d.]+)_DB', col_name)
    return float(m.group(1).replace('S', '')) if m else None

def short_clinico(name):
    m = re.search(
        r'Cluster_(\w+)_UMAP_(C\d+)_NN(\d+)_MD([\d.]+)_\w+_(S\d+)_K(\d+)_(S[\d.]+)_DB([\d.]+)_CH(\d+)',
        name
    )
    if m:
        algo, cx, nn, md, seed, k, s, db, ch = m.groups()
        return f"{algo}  |  UMAP {cx} NN{nn} MD{md} {seed}  |  K{k} {s} DB{db} CH{ch}"
    return name.replace('_C', '').replace('_', ' ')[:60]

def short_metabolico(name):
    m = re.search(
        r'Cluster_(\w+)_(C\d+)_(W\d+)_(\w+)_K(\d+)_(S[\d.]+)_DB([\d.]+)_CH(\d+)_Seed(\d+)',
        name
    )
    if m:
        algo, comp, w, method, k, s, db, ch, seed = m.groups()
        return f"{algo} {comp} {w}  |  {method} K{k} {s} DB{db} CH{ch} Seed{seed}"
    return name.replace('_M', '').replace('_', ' ')[:60]

df_top10 = df_results.sort_values(by='ARI', ascending=False).head(10).copy()
df_top10['Clinical_Short']  = df_top10['Clinical_Cluster'].apply(short_clinico)
df_top10['Metabolic_Short'] = df_top10['Metabolic_Cluster'].apply(short_metabolico)

metabolic_label = df_top10['Metabolic_Short'].iloc[0]

# ✅ CORRECCIÓN: ordenar ascendente para que el mayor ARI quede arriba en el heatmap
df_top10 = df_top10.sort_values(by='ARI', ascending=True)

# ✅ CORRECCIÓN: tamaño en pulgadas a 300 dpi (≤ 7.5 × 8.75 in)
fig, ax = plt.subplots(figsize=(W_FULL, 6.5))

hm = sns.heatmap(
    df_top10[['ARI']].values,
    annot=df_top10['ARI'].values.reshape(-1, 1),
    fmt='.3f',
    cmap='magma',
    cbar_kws={'label': 'Adjusted Rand Index', 'shrink': 0.6},
    yticklabels=df_top10['Clinical_Short'],
    xticklabels=[metabolic_label],
    linewidths=0.6,
    linecolor='white',
    annot_kws={'size': FONT_SIZE, 'weight': 'bold'},   # ✅ ≤ 12 pt
    ax=ax,
)

# ✅ colorbar con fuente controlada
cbar = hm.collections[0].colorbar
cbar.ax.tick_params(labelsize=FONT_SIZE)
cbar.set_label('Adjusted Rand Index', size=FONT_SIZE)

ax.set_yticklabels(ax.get_yticklabels(), fontsize=FONT_SIZE, rotation=0)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=FONT_SIZE, rotation=15, ha='right')
ax.set_title('Top 10 Concordances between Clinical and Metabolic Clustering',
             fontsize=FONT_TITLE, fontweight='bold', pad=10)
ax.set_ylabel('Clinical Algorithm',  fontsize=FONT_AXIS, fontweight='bold')
ax.set_xlabel('Metabolic Algorithm', fontsize=FONT_AXIS, fontweight='bold')

plt.tight_layout()
save_plos(fig, 'Fig1_heatmap_top10.tif')   # ✅ TIFF LZW en lugar de PNG

# ============================================================
# 4️⃣ IDENTIFICACIÓN DE SUBGRUPOS — HEATMAP CONTINGENCIA  (Fig 2)
# ============================================================

# ✅ CORRECCIÓN: iloc[-1] es correcto porque df_top10 está en orden ascendente
best_pair  = df_top10.iloc[-1]
best_c_col = best_pair['Clinical_Cluster']
best_m_col = best_pair['Metabolic_Cluster']

df_best = df_merged[[ID_COLUMN, best_c_col, best_m_col]].copy()
df_best = df_best[(df_best[best_c_col] != -1) & (df_best[best_m_col] != -1)]

contingency_matrix = pd.crosstab(df_best[best_c_col], df_best[best_m_col])

fig, ax = plt.subplots(figsize=(W_HALF, 4.5))   # ✅ tamaño PLOS

hm2 = sns.heatmap(
    contingency_matrix,
    annot=True,
    fmt='d',
    cmap='YlGnBu',
    linewidths=0.5,
    linecolor='white',
    annot_kws={'size': FONT_SIZE, 'weight': 'bold'},
    cbar_kws={'label': 'Number of Patients', 'shrink': 0.8},
    ax=ax,
)
cbar2 = hm2.collections[0].colorbar
cbar2.ax.tick_params(labelsize=FONT_SIZE)
cbar2.set_label('Number of Patients', size=FONT_SIZE)

ax.set_title('Patient Distribution: Clinical vs Metabolic Clusters',
             fontsize=FONT_TITLE, fontweight='bold', pad=10)
ax.set_xlabel('Metabolic Clusters (M)', fontsize=FONT_AXIS, fontweight='bold')
ax.set_ylabel('Clinical Clusters (C)',  fontsize=FONT_AXIS, fontweight='bold')

plt.tight_layout()
save_plos(fig, 'Fig2_contingency.tif')

stacked    = contingency_matrix.stack()
c_best, m_best = stacked.idxmax()
n_patients = stacked.max()

print(f"\n💡 ANALYSIS RESULT:")
print(f"Highest agreement: Clinical Cluster '{c_best}' ↔ Metabolic Cluster '{m_best}'.")
print(f"Group size: {n_patients} / {len(df_merged)} patients.")

pacientes_core = df_best[
    (df_best[best_c_col] == c_best) & (df_best[best_m_col] == m_best)
]
pacientes_core[[ID_COLUMN]].to_csv('pacientes_core_correlacion_Paretoynormas.csv', index=False)
print("✅ 'pacientes_core_correlacion_Paretoynormas.csv' generado.")

# ============================================================
# 5️⃣ ANÁLISIS DE DIVERGENCIA — PIE CHART  (Fig 3)
# ============================================================

def categorize_patient(row):
    if   row[best_c_col] == c_best and row[best_m_col] == m_best:
        return 'Core (Concordant)'
    elif row[best_c_col] == c_best and row[best_m_col] != m_best:
        return 'Metabolic Divergence'
    elif row[best_c_col] != c_best and row[best_m_col] == m_best:
        return 'Clinical Divergence'
    else:
        return 'Total Discrepancy'

df_best['Analysis_Category'] = df_best.apply(categorize_patient, axis=1)
group_summary = df_best['Analysis_Category'].value_counts()

print("\n📊 COHORT DISTRIBUTION:")
print(group_summary)

labels = group_summary.index.tolist()
values = group_summary.values
colors = sns.color_palette('pastel', len(values))

# ✅ CORRECCIÓN: explode dinámico según nº de categorías reales
n_cat   = len(values)
explode = ([0.05, 0.15, 0.15, 0.10] + [0.05] * n_cat)[:n_cat]

fig, ax = plt.subplots(figsize=(W_HALF, 5.0))   # ✅ tamaño PLOS

wedges, texts, autotexts = ax.pie(
    values,
    autopct='%1.1f%%',
    startangle=140,
    colors=colors,
    explode=explode,
    pctdistance=0.72,
    wedgeprops={'linewidth': 0.8, 'edgecolor': 'white'},
)
for at in autotexts:
    at.set_fontsize(FONT_SIZE)
    at.set_fontweight('bold')

# ✅ CORRECCIÓN: anotaciones con fuente controlada (≤ 12 pt)
offsets = [(1.35, -0.85), (-1.50, 0.90), (-1.50, 0.65), (-1.50, 0.40)]
for i, w in enumerate(wedges):
    if i >= len(offsets):
        break
    ang    = (w.theta2 + w.theta1) / 2
    x, y   = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
    xt, yt = offsets[i]
    ax.annotate(
        labels[i],
        xy=(x * 0.92, y * 0.92),
        xytext=(xt, yt),
        arrowprops=dict(arrowstyle='-', lw=1.0, color='#444444'),
        ha='left' if xt > 0 else 'right',
        va='center',
        fontsize=FONT_SIZE,     # ✅ era 13 pt — ajustado a 10 pt (PLOS)
    )

ax.set_title('Cohort Composition: Concordance vs Divergence',
             fontsize=FONT_TITLE, fontweight='bold', pad=12)   # ✅ era 16 pt

plt.tight_layout()
save_plos(fig, 'Fig3_cohort_pie.tif')

divergent_patients = df_best[df_best['Analysis_Category'] != 'Core (Concordant)']
divergent_patients.to_csv('pacientes_divergentes_para_estudio_pyn.csv', index=False)
print(f"\n✅ {len(divergent_patients)} pacientes divergentes identificados.")

# ============================================================
# 6️⃣ 3D — RELACIÓN ENTRE ALGORITMOS GANADORES  (Fig 4)
# ============================================================

print("\n🚀 3D: relación entre algoritmos ganadores...")

from sklearn.metrics import silhouette_score

df_algo             = df_best[[best_c_col, best_m_col]].copy()
df_algo[best_c_col] = df_algo[best_c_col].astype(int)
df_algo[best_m_col] = df_algo[best_m_col].astype(int)

# ✅ Silhouette extraído del nombre de columna — es el valor original del entrenamiento
sil_c = extract_silhouette_clinico(best_c_col)
sil_m = extract_silhouette_metabolico(best_m_col)
sil_c_str = f"{sil_c:.3f}" if sil_c is not None else "N/A"
sil_m_str = f"{sil_m:.3f}" if sil_m is not None else "N/A"
print(f"  Silhouette Clinical  (del nombre de columna): {sil_c_str}")
print(f"  Silhouette Metabolic (del nombre de columna): {sil_m_str}")

counts = df_algo.groupby([best_c_col, best_m_col]).size().reset_index(name='Count')

fig = plt.figure(figsize=(4.0, 4.0))
ax  = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    counts[best_c_col],
    counts[best_m_col],
    counts['Count'],
    s=counts['Count'] * 2.5,
    c=counts['Count'],
    cmap='viridis',
    alpha=0.85,
    edgecolors='white',
    linewidths=0.4,
)
cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.10)
cb.ax.tick_params(labelsize=FONT_SIZE)
cb.set_label('Number of Patients', size=FONT_SIZE)

ax.set_xlabel('Clinical Cluster',  fontsize=FONT_SIZE, labelpad=8)
ax.set_ylabel('Metabolic Cluster', fontsize=FONT_SIZE, labelpad=8)
ax.set_zlabel('No. of Patients',   fontsize=FONT_SIZE, labelpad=8)
ax.tick_params(labelsize=FONT_SIZE - 1)
ax.set_title(
    f'3D Patient Distribution — Between Algorithms\n'
    f'(KMeans | UMAP C3 NN10 MD0.05 S123 | K2 | Sil C={sil_c_str} / M={sil_m_str})',
    fontsize=FONT_SIZE, fontweight='bold'
)
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
save_plos(fig, 'Fig4_3D_relationship.tif')

# ============================================================
# 7️⃣ 3D — ALGORITMO CLÍNICO  (Fig 5)
# ============================================================

print("\n🚀 3D: algoritmo clínico...")

df_clin = df_best[[ID_COLUMN, best_c_col]].copy()
df_clin = df_clin.sort_values(by=best_c_col).reset_index(drop=True)
df_clin['X'] = np.arange(len(df_clin))
df_clin['Y'] = df_clin[best_c_col].astype(int)
df_clin['Z'] = np.random.normal(0, 0.2, len(df_clin))

clusters_clin = np.sort(df_clin['Y'].unique())
palette_clin  = plt.cm.tab10(np.linspace(0, 0.6, len(clusters_clin)))

fig = plt.figure(figsize=(4.0, 4.0))
ax  = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_clin, palette_clin):
    subset = df_clin[df_clin['Y'] == c]
    ax.scatter(subset['X'], subset['Y'], subset['Z'],
               label=f'Cluster {c}', color=color,
               alpha=0.80, s=18, edgecolors='none')

ax.set_xlabel('Patients',            fontsize=FONT_SIZE, labelpad=8)
ax.set_ylabel('Cluster Label',       fontsize=FONT_SIZE, labelpad=8)
ax.set_zlabel('Jitter (arb. units)', fontsize=FONT_SIZE, labelpad=8)
ax.tick_params(labelsize=FONT_SIZE - 1)
ax.set_title(
    f'Clinical Algorithm — 3D Cluster Distribution\n'
    f'(KMeans | UMAP C3 NN10 MD0.05 S123 | K2 | Silhouette = {sil_c_str})',
    fontsize=FONT_SIZE, fontweight='bold'
)
ax.legend(title='Clusters', title_fontsize=FONT_SIZE,
          fontsize=FONT_SIZE, loc='upper left',
          framealpha=0.7, edgecolor='#cccccc')

plt.tight_layout()
save_plos(fig, 'Fig5_3D_clinical.tif')

# ============================================================
# 8️⃣ 3D — ALGORITMO METABÓLICO  (Fig 6)
# ============================================================

print("\n🚀 3D: algoritmo metabólico...")

df_met = df_best[[ID_COLUMN, best_m_col]].copy()
df_met = df_met.sort_values(by=best_m_col).reset_index(drop=True)
df_met['X'] = np.arange(len(df_met))
df_met['Y'] = df_met[best_m_col].astype(int)
df_met['Z'] = np.random.normal(0, 0.2, len(df_met))

clusters_met = np.sort(df_met['Y'].unique())
palette_met  = plt.cm.tab10(np.linspace(0.3, 0.9, len(clusters_met)))

fig = plt.figure(figsize=(4.0, 4.0))
ax  = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_met, palette_met):
    subset = df_met[df_met['Y'] == c]
    ax.scatter(subset['X'], subset['Y'], subset['Z'],
               label=f'Cluster {c}', color=color,
               alpha=0.80, s=18, edgecolors='none')

ax.set_xlabel('Patients',            fontsize=FONT_SIZE, labelpad=8)
ax.set_ylabel('Cluster Label',       fontsize=FONT_SIZE, labelpad=8)
ax.set_zlabel('Jitter (arb. units)', fontsize=FONT_SIZE, labelpad=8)
ax.tick_params(labelsize=FONT_SIZE - 1)
ax.set_title(
    f'Metabolic Algorithm — 3D Cluster Distribution\n'
    f'(KMeans | UMAP C3 NN10 MD0.05 S123 | K2 | Silhouette = {sil_m_str})',
    fontsize=FONT_SIZE, fontweight='bold'
)
ax.legend(title='Clusters', title_fontsize=FONT_SIZE,
          fontsize=FONT_SIZE, loc='upper left',
          framealpha=0.7, edgecolor='#cccccc')

plt.tight_layout()
save_plos(fig, 'Fig6_3D_metabolic.tif')

# ============================================================
# RESUMEN
# ============================================================
print("\n" + "=" * 55)
print("  TODAS LAS FIGURAS GENERADAS")
print("=" * 55)
print(f"  DPI      : {DPI}")
print(f"  Formato  : TIFF (LZW)")
print(f"  Fuente   : {FONT_FAMILY}, 8–12 pt")
print(f"  Color    : RGB")
print(f"  Archivos : Fig1–Fig6.tif")
print("=" * 55)