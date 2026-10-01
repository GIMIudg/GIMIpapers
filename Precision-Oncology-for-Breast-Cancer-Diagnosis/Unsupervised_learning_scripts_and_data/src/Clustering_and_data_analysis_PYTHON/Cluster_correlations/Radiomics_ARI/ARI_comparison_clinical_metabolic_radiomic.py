import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
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

warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# ============================================================
# PLOS ONE — ESPECIFICACIONES GLOBALES
# ============================================================
DPI        = 300
FONT_SIZE  = 10
FONT_TITLE = 12
FONT_AXIS  = 11
W_FULL     = 7.5
W_HALF     = 5.2
H_MAX      = 8.75

_available  = {f.name for f in fm.fontManager.ttflist}
FONT_FAMILY = 'Arial' if 'Arial' in _available else 'DejaVu Sans'
print(f"  Fuente: {FONT_FAMILY}"
      + ("  ✔" if FONT_FAMILY == 'Arial'
         else "  ⚠ Arial no encontrada — usando DejaVu Sans"))

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
# 1️⃣  CONFIGURACIÓN Y CARGA DE DATOS
# ============================================================

ID_COLUMN         = 'ModelName'
MIN_VALID_SAMPLES = 50

path1 = (
    "/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis"
    "/src/Clustering and data analysis PYTHON/Clinical data analysis"
    "/ML models using clinical data"
    "/Results_clustering_UMAP_seleccion_reducida_conmerge-DATOSNUEVOS"
    "/pacientes_clusterizados_todos_sinfiltro.csv"
)
path2 = (
    "/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis"
    "/src/Clustering and data analysis PYTHON/Metabolic data analysis"
    "/ML models using metabolic data"
    "/resultados_TumorPhenotype_PCA_metrics_actualizado_sinl2"
    "/PatientClusters_TumorPhenotype_PCA.csv"
)
path3 = (
    "/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis"
    "/src/Clustering and data analysis PYTHON/Radiomics data analysis"
    "/Results_clustering_PCA_seleccion_reducida_conmerge-DATOSNUEVOS"
    "/pacientes_clusterizados_todos_sinfiltro.csv"
)


def load_and_standardize_clusters(file_path, suffix, id_col='ModelName'):
    """
    Carga el CSV, limpia los IDs TCGA a 16 chars y añade sufijo a las
    columnas de cluster para distinguir el dominio (_C / _M / _R).
    """
    try:
        df = pd.read_csv(file_path, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"No se pudo leer {file_path}. Error: {e}")

    df.columns = df.columns.str.strip()
    if id_col not in df.columns:
        raise ValueError(f"Columna ID '{id_col}' no encontrada en: {file_path}")

    # Normalizar ID a 16 chars (TCGA-XX-XXXX-XXX)
    df[id_col] = (
        df[id_col].astype(str)
        .str.split('_').str[0]
        .str.split('.').str[0]
        .str.slice(0, 16)
    )

    cluster_cols   = [col for col in df.columns if col != id_col]
    rename_mapping = {col: f"{col}_{suffix}" for col in cluster_cols}
    df = df.rename(columns=rename_mapping)

    print(f"  [{suffix}] {len(df)} muestras | {len(cluster_cols)} columnas de cluster  ←  {file_path.split('/')[-1]}")
    return df[[id_col] + list(rename_mapping.values())].copy()


try:
    df_clin = load_and_standardize_clusters(path1, suffix='C', id_col=ID_COLUMN)
    df_meta = load_and_standardize_clusters(path2, suffix='M', id_col=ID_COLUMN)
    df_rad  = load_and_standardize_clusters(path3, suffix='R', id_col=ID_COLUMN)

    # Merge triple: inner en los tres dominios
    df_cm  = df_clin.merge(df_meta, on=ID_COLUMN, how='inner')
    df_merged = df_cm.merge(df_rad,  on=ID_COLUMN, how='inner')

    n_cm  = len(df_clin.merge(df_meta, on=ID_COLUMN, how='inner'))
    n_cr  = len(df_clin.merge(df_rad,  on=ID_COLUMN, how='inner'))
    n_mr  = len(df_meta.merge(df_rad,  on=ID_COLUMN, how='inner'))
    print(f"\n✅ Merge triple completado.")
    print(f"   C ∩ M : {n_cm} pacientes")
    print(f"   C ∩ R : {n_cr} pacientes")
    print(f"   M ∩ R : {n_mr} pacientes")
    print(f"   C ∩ M ∩ R : {len(df_merged)} pacientes  ← usado para todas las comparaciones")

except Exception as e:
    print(f"❌ ERROR: {e}")
    sys.exit()

# Separar columnas por dominio
cluster_cols_C = [col for col in df_merged.columns if col.endswith('_C')]
cluster_cols_M = [col for col in df_merged.columns if col.endswith('_M')]
cluster_cols_R = [col for col in df_merged.columns if col.endswith('_R')]

print(f"\n   Columnas _C : {len(cluster_cols_C)}")
print(f"   Columnas _M : {len(cluster_cols_M)}")
print(f"   Columnas _R : {len(cluster_cols_R)}")

# ============================================================
# 2️⃣  CÁLCULO DE MÉTRICAS (ARI / AMI) — TODAS LAS PAREJAS
# ============================================================

def calculate_metrics(df, col1, col2, pair_type):
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
        'Col_A':      col1,
        'Col_B':      col2,
        'Pair_Type':  pair_type,        # 'C_vs_M' | 'C_vs_R' | 'M_vs_R'
        'ARI':        adjusted_rand_score(labels1, labels2),
        'AMI':        adjusted_mutual_info_score(labels1, labels2),
        'N_Samples':  len(valid_idx),
    }


print("\n⏳ Calculando métricas ARI/AMI para los tres pares de dominios...")

results_list = []

# Clínico vs Metabólico
for c in cluster_cols_C:
    for m in cluster_cols_M:
        res = calculate_metrics(df_merged, c, m, 'C_vs_M')
        if res:
            results_list.append(res)

# Clínico vs Radiómico
for c in cluster_cols_C:
    for r in cluster_cols_R:
        res = calculate_metrics(df_merged, c, r, 'C_vs_R')
        if res:
            results_list.append(res)

# Metabólico vs Radiómico
for m in cluster_cols_M:
    for r in cluster_cols_R:
        res = calculate_metrics(df_merged, m, r, 'M_vs_R')
        if res:
            results_list.append(res)

df_results = pd.DataFrame(results_list)

print(f"✅ {len(df_results)} pares válidos calculados.")
for pt in ['C_vs_M', 'C_vs_R', 'M_vs_R']:
    sub = df_results[df_results['Pair_Type'] == pt]
    if not sub.empty:
        print(f"   {pt}: {len(sub)} pares | ARI máx = {sub['ARI'].max():.4f}")

# ============================================================
# HELPERS: etiquetas cortas por dominio
# ============================================================

def short_label_C(name):
    """Clínico: patrón UMAP — Cluster_ALGO_UMAP_Cx_NNx_MDx_Sx_Kx_Sx_DBx_CHx_C"""
    m = re.search(
        r'Cluster_(\w+?)_UMAP_(C\d+)_NN(\d+)_MD([\d.]+)_\w+_(S\d+)_K(\d+)_(S[\d.]+)_DB([\d.]+)_CH(\d+)',
        name
    )
    if m:
        algo, cx, nn, md, seed, k, s, db, ch = m.groups()
        return f"{algo} {cx} NN{nn} MD{md} {seed} K{k} {s}"
    return name.replace('_C', '').replace('_', ' ')[:55]


def short_label_M(name):
    """Metabólico: patrón PCA — Cluster_ALGO_Cx_Wx_METHOD_Kx_Sx_DBx_CHx_Seedx_M"""
    m = re.search(
        r'Cluster_(\w+?)_(C\d+)_(W\d+)_(\w+)_K(\d+)_(S[\d.]+)_DB([\d.]+)_CH(\d+)_Seed(\d+)',
        name
    )
    if m:
        algo, comp, w, method, k, s, db, ch, seed = m.groups()
        return f"{algo} {comp} {w} {method} K{k} {s}"
    return name.replace('_M', '').replace('_', ' ')[:55]


def short_label_R(name):
    """Radiómico: patrón PCA — Cluster_ALGO_PCA_Cx_Wx_Kx_Sx_DBx_CHx_Sx_R"""
    m = re.search(
        r'Cluster_(\w+?)_PCA_(C\d+)_(W\d+)_K(\d+)_(S[\d.]+)_DB([\d.]+)_CH(\d+)(?:_(S\d+))?',
        name
    )
    if m:
        groups = m.groups()
        algo, cx, w, k, s, db, ch = groups[:7]
        seed = groups[7] if groups[7] else ''
        return f"{algo} {cx} {w} K{k} {s} {seed}".strip()
    return name.replace('_R', '').replace('_', ' ')[:55]


def short_label(name):
    """Dispatcher automático según sufijo."""
    if name.endswith('_C'):
        return short_label_C(name)
    elif name.endswith('_M'):
        return short_label_M(name)
    elif name.endswith('_R'):
        return short_label_R(name)
    return name[:55]


def extract_silhouette(col_name):
    """Extrae Silhouette Score del nombre de cualquier columna de cluster."""
    m = re.search(r'_(S[\d.]+)_DB', col_name)
    return float(m.group(1).replace('S', '')) if m else None

# ============================================================
# 3️⃣  TOP 10 POR PAR — HEATMAP COMPARATIVO  (Fig 1)
#     Un panel por cada uno de los tres pares de dominios.
# ============================================================

PAIR_CONFIG = {
    'C_vs_M': {
        'title_y': 'Algoritmo Clínico (_C)',
        'title_x': 'Algoritmo Metabólico (_M)',
        'col_a':   'Col_A', 'col_b': 'Col_B',
    },
    'C_vs_R': {
        'title_y': 'Algoritmo Clínico (_C)',
        'title_x': 'Algoritmo Radiómico (_R)',
        'col_a':   'Col_A', 'col_b': 'Col_B',
    },
    'M_vs_R': {
        'title_y': 'Algoritmo Metabólico (_M)',
        'title_x': 'Algoritmo Radiómico (_R)',
        'col_a':   'Col_A', 'col_b': 'Col_B',
    },
}

fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 7.0))
fig.suptitle(
    'Top 10 Concordancias ARI entre Dominios de Clustering',
    fontsize=FONT_TITLE, fontweight='bold', y=1.01
)

best_per_pair = {}   # guardará el mejor par por tipo para usar en figuras siguientes

for ax, (pair_type, cfg) in zip(axes, PAIR_CONFIG.items()):
    sub = df_results[df_results['Pair_Type'] == pair_type].copy()

    if sub.empty:
        ax.set_visible(False)
        continue

    top10 = (sub.sort_values('ARI', ascending=False)
               .head(10)
               .sort_values('ARI', ascending=True)   # mayor ARI arriba en heatmap
               .copy())

    top10['Label_A'] = top10[cfg['col_a']].apply(short_label)
    top10['Label_B'] = top10[cfg['col_b']].apply(short_label)

    # Guardar el mejor par del tipo para figuras siguientes
    best_row = top10.iloc[-1]
    best_per_pair[pair_type] = {
        'col_a':   best_row[cfg['col_a']],
        'col_b':   best_row[cfg['col_b']],
        'label_a': best_row['Label_A'],
        'label_b': best_row['Label_B'],
        'ARI':     best_row['ARI'],
        'AMI':     best_row['AMI'],
    }

    x_label = top10['Label_B'].iloc[0]

    sns.heatmap(
        top10[['ARI']].values,
        annot=top10['ARI'].values.reshape(-1, 1),
        fmt='.3f',
        cmap='magma',
        cbar_kws={'label': 'ARI', 'shrink': 0.55},
        yticklabels=top10['Label_A'].tolist(),
        xticklabels=[x_label],
        linewidths=0.5,
        linecolor='white',
        annot_kws={'size': FONT_SIZE - 1, 'weight': 'bold'},
        ax=ax,
    )

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=FONT_SIZE - 1, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=FONT_SIZE - 1, rotation=20, ha='right')
    ax.set_title(pair_type.replace('_', ' '), fontsize=FONT_TITLE, fontweight='bold', pad=8)
    ax.set_ylabel(cfg['title_y'], fontsize=FONT_SIZE)
    ax.set_xlabel(cfg['title_x'], fontsize=FONT_SIZE)

plt.tight_layout()
save_plos(fig, 'Fig1_heatmap_top10_tres_dominios.tif')

# ============================================================
# 4️⃣  HEATMAP DE CONTINGENCIA — UN PANEL POR PAR  (Fig 2)
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 4.5))
fig.suptitle(
    'Distribución de Pacientes — Matrices de Contingencia',
    fontsize=FONT_TITLE, fontweight='bold', y=1.01
)

for ax, (pair_type, best) in zip(axes, best_per_pair.items()):
    col_a, col_b = best['col_a'], best['col_b']
    df_pair = df_merged[[ID_COLUMN, col_a, col_b]].copy()
    df_pair = df_pair[(df_pair[col_a] != -1) & (df_pair[col_b] != -1)]
    cont    = pd.crosstab(df_pair[col_a], df_pair[col_b])

    sns.heatmap(
        cont,
        annot=True, fmt='d',
        cmap='YlGnBu',
        linewidths=0.4, linecolor='white',
        annot_kws={'size': FONT_SIZE - 1, 'weight': 'bold'},
        cbar_kws={'label': 'N Pacientes', 'shrink': 0.75},
        ax=ax,
    )
    ax.set_title(pair_type.replace('_', ' '), fontsize=FONT_TITLE, fontweight='bold', pad=8)

    suffix_a = col_a[-1]   # C / M / R
    suffix_b = col_b[-1]
    ax.set_ylabel(f'Cluster ({suffix_a})', fontsize=FONT_SIZE)
    ax.set_xlabel(f'Cluster ({suffix_b})', fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE - 1)

plt.tight_layout()
save_plos(fig, 'Fig2_contingency_tres_dominios.tif')

# Identificar el par con mayor ARI global para el análisis de pacientes core
best_global_type = max(best_per_pair, key=lambda k: best_per_pair[k]['ARI'])
best_global      = best_per_pair[best_global_type]
best_col_a       = best_global['col_a']
best_col_b       = best_global['col_b']

print(f"\n💡 Par con mayor ARI global: {best_global_type}")
print(f"   Col A : {best_col_a}")
print(f"   Col B : {best_col_b}")
print(f"   ARI   : {best_global['ARI']:.4f} | AMI: {best_global['AMI']:.4f}")

df_best = df_merged[[ID_COLUMN, best_col_a, best_col_b]].copy()
df_best = df_best[(df_best[best_col_a] != -1) & (df_best[best_col_b] != -1)]

contingency_matrix = pd.crosstab(df_best[best_col_a], df_best[best_col_b])
stacked  = contingency_matrix.stack()
c_best, m_best = stacked.idxmax()
n_patients = stacked.max()

print(f"   Mayor concordancia intra-par: [{c_best}] ↔ [{m_best}]  ({n_patients} pacientes)")

pacientes_core = df_best[
    (df_best[best_col_a] == c_best) &
    (df_best[best_col_b] == m_best)
]
pacientes_core[[ID_COLUMN]].to_csv(
    'pacientes_core_correlacion_Paretoynormas.csv', index=False
)
print("✅  'pacientes_core_correlacion_Paretoynormas.csv' generado.")

# ============================================================
# 5️⃣  DISTRIBUCIÓN ARI POR DOMINIO — BOXPLOT  (Fig 3)
#     (reemplaza el pie chart original; muestra la distribución
#      completa de ARI para los tres pares, más informativo)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 4.5))

# ── Panel izquierdo: boxplot ARI por par ─────────────────────
pair_order  = ['C_vs_M', 'C_vs_R', 'M_vs_R']
pair_labels = ['Clínico\nvs\nMetabólico', 'Clínico\nvs\nRadiómico', 'Metabólico\nvs\nRadiómico']
palette_box = sns.color_palette('Set2', 3)

sns.boxplot(
    data=df_results,
    x='Pair_Type', y='ARI',
    order=pair_order,
    palette=palette_box,
    width=0.5,
    linewidth=1.2,
    ax=axes[0],
)
sns.stripplot(
    data=df_results,
    x='Pair_Type', y='ARI',
    order=pair_order,
    color='#333333', alpha=0.25, size=2.5, jitter=True,
    ax=axes[0],
)

axes[0].set_xticklabels(pair_labels, fontsize=FONT_SIZE)
axes[0].set_xlabel('Par de Dominios',      fontsize=FONT_AXIS, fontweight='bold')
axes[0].set_ylabel('ARI',                  fontsize=FONT_AXIS, fontweight='bold')
axes[0].set_title('Distribución de ARI\npor Par de Dominios',
                   fontsize=FONT_TITLE, fontweight='bold')
axes[0].axhline(0, color='red', lw=0.8, ls='--', alpha=0.6, label='ARI = 0 (azar)')
axes[0].legend(fontsize=FONT_SIZE - 1)

# ── Panel derecho: AMI por par ────────────────────────────────
sns.boxplot(
    data=df_results,
    x='Pair_Type', y='AMI',
    order=pair_order,
    palette=palette_box,
    width=0.5,
    linewidth=1.2,
    ax=axes[1],
)
sns.stripplot(
    data=df_results,
    x='Pair_Type', y='AMI',
    order=pair_order,
    color='#333333', alpha=0.25, size=2.5, jitter=True,
    ax=axes[1],
)

axes[1].set_xticklabels(pair_labels, fontsize=FONT_SIZE)
axes[1].set_xlabel('Par de Dominios',      fontsize=FONT_AXIS, fontweight='bold')
axes[1].set_ylabel('AMI',                  fontsize=FONT_AXIS, fontweight='bold')
axes[1].set_title('Distribución de AMI\npor Par de Dominios',
                   fontsize=FONT_TITLE, fontweight='bold')
axes[1].axhline(0, color='red', lw=0.8, ls='--', alpha=0.6, label='AMI = 0 (azar)')
axes[1].legend(fontsize=FONT_SIZE - 1)

plt.tight_layout()
save_plos(fig, 'Fig3_ARI_AMI_boxplot_dominios.tif')

# ============================================================
# 6️⃣  PIE CHART — CATEGORÍAS DE CONCORDANCIA  (Fig 4)
#     (basado en el mejor par global, igual que el original)
# ============================================================

def categorize_patient(row):
    if row[best_col_a] == c_best and row[best_col_b] == m_best:
        return 'Núcleo (Concordante)'
    elif row[best_col_a] == c_best and row[best_col_b] != m_best:
        return 'Divergencia B'
    elif row[best_col_a] != c_best and row[best_col_b] == m_best:
        return 'Divergencia A'
    else:
        return 'Discrepancia Total'

df_best['Analysis_Category'] = df_best.apply(categorize_patient, axis=1)
group_summary = df_best['Analysis_Category'].value_counts()

print("\n📊 DISTRIBUCIÓN DE LA COHORTE (par de mayor ARI):")
print(group_summary)

labels_pie = group_summary.index.tolist()
values_pie = group_summary.values
colors_pie = sns.color_palette('pastel', len(values_pie))
n_cat      = len(values_pie)
explode    = ([0.05, 0.15, 0.15, 0.10] + [0.05] * n_cat)[:n_cat]

fig, ax = plt.subplots(figsize=(W_HALF, 5.0))

wedges, texts, autotexts = ax.pie(
    values_pie,
    autopct='%1.1f%%',
    startangle=140,
    colors=colors_pie,
    explode=explode,
    pctdistance=0.72,
    wedgeprops={'linewidth': 0.8, 'edgecolor': 'white'},
)

for at in autotexts:
    at.set_fontsize(FONT_SIZE)
    at.set_fontweight('bold')

offsets = [(1.35, -0.85), (-1.50, 0.90), (-1.50, 0.65), (-1.50, 0.40)]
for i, w in enumerate(wedges):
    if i >= len(offsets):
        break
    ang     = (w.theta2 + w.theta1) / 2
    x, y    = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
    xt, yt  = offsets[i]
    ax.annotate(
        labels_pie[i],
        xy=(x * 0.92, y * 0.92),
        xytext=(xt, yt),
        arrowprops=dict(arrowstyle='-', lw=1.0, color='#444444'),
        ha='left' if xt > 0 else 'right',
        va='center',
        fontsize=FONT_SIZE,
    )

suf_a = best_col_a[-1]
suf_b = best_col_b[-1]
ax.set_title(
    f'Concordancia vs Divergencia — Par {best_global_type.replace("_", " ")}\n'
    f'(mejor ARI: {best_global["ARI"]:.3f})',
    fontsize=FONT_TITLE, fontweight='bold', pad=12
)

plt.tight_layout()
save_plos(fig, 'Fig4_cohort_pie.tif')

divergent_patients = df_best[df_best['Analysis_Category'] != 'Núcleo (Concordante)']
divergent_patients.to_csv('pacientes_divergentes_para_estudio_pyn.csv', index=False)
print(f"\n✅ {len(divergent_patients)} pacientes divergentes identificados.")

# ============================================================
# 7️⃣  3D — RELACIÓN ENTRE ALGORITMOS GANADORES (mejor par global)  (Fig 5)
# ============================================================

print("\n🚀 3D: relación entre algoritmos ganadores (mejor par global)...")

df_algo = df_best[[best_col_a, best_col_b]].copy()
df_algo[best_col_a] = df_algo[best_col_a].astype(int)
df_algo[best_col_b] = df_algo[best_col_b].astype(int)

sil_a     = extract_silhouette(best_col_a)
sil_b     = extract_silhouette(best_col_b)
sil_a_str = f"{sil_a:.3f}" if sil_a is not None else "N/A"
sil_b_str = f"{sil_b:.3f}" if sil_b is not None else "N/A"

counts = df_algo.groupby([best_col_a, best_col_b]).size().reset_index(name='Count')

fig = plt.figure(figsize=(4.5, 4.5))
ax  = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    counts[best_col_a], counts[best_col_b], counts['Count'],
    s=counts['Count'] * 2.5,
    c=counts['Count'], cmap='viridis',
    alpha=0.85, edgecolors='white', linewidths=0.4,
)

cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.10)
cb.ax.tick_params(labelsize=FONT_SIZE)
cb.set_label('Número de Pacientes', size=FONT_SIZE)

ax.set_xlabel(f'Cluster ({best_col_a[-1]})', fontsize=FONT_SIZE, labelpad=8)
ax.set_ylabel(f'Cluster ({best_col_b[-1]})', fontsize=FONT_SIZE, labelpad=8)
ax.set_zlabel('N Pacientes',                 fontsize=FONT_SIZE, labelpad=8)
ax.tick_params(labelsize=FONT_SIZE - 1)

ax.set_title(
    f'Distribución 3D — {best_global_type.replace("_", " ")}\n'
    f'Sil ({best_col_a[-1]})={sil_a_str}  |  Sil ({best_col_b[-1]})={sil_b_str}  |  ARI={best_global["ARI"]:.3f}',
    fontsize=FONT_SIZE, fontweight='bold'
)

ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
save_plos(fig, 'Fig5_3D_best_pair.tif')

# ============================================================
# 8️⃣  3D — DISTRIBUCIÓN POR DOMINIO: C / M / R  (Fig 6, 7, 8)
# ============================================================

def plot_3d_domain(df_src, id_col, cluster_col, domain_label,
                   sil_str, filename, palette_range=(0, 0.6)):
    df_d = df_src[[id_col, cluster_col]].copy()
    df_d = df_d[df_d[cluster_col] != -1].sort_values(by=cluster_col).reset_index(drop=True)
    df_d['X'] = np.arange(len(df_d))
    df_d['Y'] = df_d[cluster_col].astype(int)
    df_d['Z'] = np.random.normal(0, 0.2, len(df_d))

    clusters = np.sort(df_d['Y'].unique())
    palette  = plt.cm.tab10(np.linspace(*palette_range, len(clusters)))

    fig = plt.figure(figsize=(4.0, 4.0))
    ax  = fig.add_subplot(111, projection='3d')

    for c, color in zip(clusters, palette):
        subset = df_d[df_d['Y'] == c]
        ax.scatter(subset['X'], subset['Y'], subset['Z'],
                   label=f'Cluster {c}', color=color,
                   alpha=0.80, s=18, edgecolors='none')

    ax.set_xlabel('Pacientes',            fontsize=FONT_SIZE, labelpad=8)
    ax.set_ylabel('Etiqueta del Cluster', fontsize=FONT_SIZE, labelpad=8)
    ax.set_zlabel('Jitter (arb.)',        fontsize=FONT_SIZE, labelpad=8)
    ax.tick_params(labelsize=FONT_SIZE - 1)
    ax.set_title(
        f'Algoritmo {domain_label} — Distribución 3D\n'
        f'Silhouette = {sil_str}',
        fontsize=FONT_SIZE, fontweight='bold'
    )
    ax.legend(title='Clusters', title_fontsize=FONT_SIZE, fontsize=FONT_SIZE,
              loc='upper left', framealpha=0.7, edgecolor='#cccccc')

    plt.tight_layout()
    save_plos(fig, filename)


# Para cada par ganador, graficar los dos dominios implicados
for pair_type, best in best_per_pair.items():
    col_a      = best['col_a']
    col_b      = best['col_b']
    suf_a      = col_a[-1]        # C / M / R
    suf_b      = col_b[-1]
    domain_map = {'C': 'Clínico', 'M': 'Metabólico', 'R': 'Radiómico'}
    sil_a_s    = f"{extract_silhouette(col_a):.3f}" if extract_silhouette(col_a) else "N/A"
    sil_b_s    = f"{extract_silhouette(col_b):.3f}" if extract_silhouette(col_b) else "N/A"

    df_pair = df_merged[[ID_COLUMN, col_a, col_b]].copy()
    df_pair = df_pair[(df_pair[col_a] != -1) & (df_pair[col_b] != -1)]

    tag = pair_type.replace('_vs_', '_')
    print(f"\n🚀 3D: {pair_type}...")
    plot_3d_domain(df_pair, ID_COLUMN, col_a,
                   domain_map[suf_a], sil_a_s,
                   f'Fig_3D_{tag}_{suf_a}.tif',
                   palette_range=(0.0, 0.6))
    plot_3d_domain(df_pair, ID_COLUMN, col_b,
                   domain_map[suf_b], sil_b_s,
                   f'Fig_3D_{tag}_{suf_b}.tif',
                   palette_range=(0.3, 0.9))

# ============================================================
# RESUMEN FINAL
# ============================================================

print("\n" + "=" * 60)
print("  TODAS LAS FIGURAS GENERADAS")
print("=" * 60)
print(f"  DPI      : {DPI}")
print(f"  Formato  : TIFF (LZW)")
print(f"  Fuente   : {FONT_FAMILY}, 8–12 pt")
print(f"  Color    : RGB")
print(f"\n  Figuras fijas:")
print(f"    Fig1_heatmap_top10_tres_dominios.tif")
print(f"    Fig2_contingency_tres_dominios.tif")
print(f"    Fig3_ARI_AMI_boxplot_dominios.tif")
print(f"    Fig4_cohort_pie.tif")
print(f"    Fig5_3D_best_pair.tif")
print(f"\n  Figuras 3D por par ganador:")
for pair_type, best in best_per_pair.items():
    tag  = pair_type.replace('_vs_', '_')
    suf_a = best['col_a'][-1]
    suf_b = best['col_b'][-1]
    print(f"    Fig_3D_{tag}_{suf_a}.tif  |  Fig_3D_{tag}_{suf_b}.tif")

print(f"\n  Archivos CSV:")
print(f"    pacientes_core_correlacion_Paretoynormas.csv")
print(f"    pacientes_divergentes_para_estudio_pyn.csv")
print("=" * 60)