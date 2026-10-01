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
# 1️⃣  CONFIGURACIÓN Y CARGA DE DATOS (solo M y R)
# ============================================================

ID_COLUMN         = 'ModelName'
MIN_VALID_SAMPLES = 50

# --- Ajusta estas dos rutas a tus archivos ---
path_meta = (
    "/Users/eduardoruiz/Documents/GIMIpapers/Precision-Oncology-for-Breast-Cancer-Diagnosis/Unsupervised_learning_scripts_and_data/src/Clustering_and_data_analysis_PYTHON/Metabolic_data_analysis/ML_models_using_metabolic_data/results_TumorPhenotype_PCA_metrics_updated_nol2/PatientClusters_TumorPhenotype_PCA.csv"
)
path_rad = (
    "/Users/eduardoruiz/Documents/GIMIpapers/Precision-Oncology-for-Breast-Cancer-Diagnosis/Unsupervised_learning_scripts_and_data/src/Clustering_and_data_analysis_PYTHON/Radiomics_data_analysis/Results_clustering_PCA_radiomics/pacientes_clusterizados_todos_sinfiltro.csv"
)


def load_and_standardize_clusters(file_path, suffix, id_col='ModelName'):
    """
    Carga el CSV, limpia los IDs TCGA a 16 chars y añade sufijo a las
    columnas de cluster para distinguir el dominio (_M / _R).
    """
    try:
        df = pd.read_csv(file_path, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"No se pudo leer {file_path}. Error: {e}")

    df.columns = df.columns.str.strip()
    if id_col not in df.columns:
        raise ValueError(f"Columna ID '{id_col}' no encontrada en: {file_path}")

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
    df_meta = load_and_standardize_clusters(path_meta, suffix='M', id_col=ID_COLUMN)
    df_rad  = load_and_standardize_clusters(path_rad,  suffix='R', id_col=ID_COLUMN)

    # Merge simple: inner entre M y R
    df_merged = df_meta.merge(df_rad, on=ID_COLUMN, how='inner')

    print(f"\n✅ Merge M ∩ R completado.")
    print(f"   M ∩ R : {len(df_merged)} pacientes  ← usado para todas las comparaciones")

except Exception as e:
    print(f"❌ ERROR: {e}")
    sys.exit()

# Separar columnas por dominio
cluster_cols_M = [col for col in df_merged.columns if col.endswith('_M')]
cluster_cols_R = [col for col in df_merged.columns if col.endswith('_R')]

print(f"\n   Columnas _M : {len(cluster_cols_M)}")
print(f"   Columnas _R : {len(cluster_cols_R)}")

# ============================================================
# 2️⃣  CÁLCULO DE MÉTRICAS (ARI / AMI) — M vs R
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
        'Pair_Type':  pair_type,        # 'M_vs_R'
        'ARI':        adjusted_rand_score(labels1, labels2),
        'AMI':        adjusted_mutual_info_score(labels1, labels2),
        'N_Samples':  len(valid_idx),
    }


print("\n⏳ Calculando métricas ARI/AMI para Metabólico vs Radiómico...")

results_list = []
for m in cluster_cols_M:
    for r in cluster_cols_R:
        res = calculate_metrics(df_merged, m, r, 'M_vs_R')
        if res:
            results_list.append(res)

df_results = pd.DataFrame(results_list)

print(f"✅ {len(df_results)} pares válidos calculados.")
if not df_results.empty:
    print(f"   M_vs_R: {len(df_results)} pares | ARI máx = {df_results['ARI'].max():.4f}")
else:
    print("⚠️  No se encontraron pares válidos. Revisa MIN_VALID_SAMPLES o los datos de entrada.")
    sys.exit()

# ============================================================
# HELPERS: etiquetas cortas por dominio
# ============================================================

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
    if name.endswith('_M'):
        return short_label_M(name)
    elif name.endswith('_R'):
        return short_label_R(name)
    return name[:55]


def extract_silhouette(col_name):
    """Extrae Silhouette Score del nombre de cualquier columna de cluster."""
    m = re.search(r'_(S[\d.]+)_DB', col_name)
    return float(m.group(1).replace('S', '')) if m else None

# ============================================================
# 3️⃣  TOP 10 — HEATMAP COMPARATIVO M vs R  (Fig 1)
# ============================================================

top10 = (df_results.sort_values('ARI', ascending=False)
           .head(10)
           .sort_values('ARI', ascending=True)   # mayor ARI arriba en heatmap
           .copy())

top10['Label_A'] = top10['Col_A'].apply(short_label)
top10['Label_B'] = top10['Col_B'].apply(short_label)

best_row = top10.iloc[-1]
best = {
    'col_a':   best_row['Col_A'],
    'col_b':   best_row['Col_B'],
    'label_a': best_row['Label_A'],
    'label_b': best_row['Label_B'],
    'ARI':     best_row['ARI'],
    'AMI':     best_row['AMI'],
}

x_label = top10['Label_B'].iloc[0]

fig, ax = plt.subplots(figsize=(W_HALF, 6.0))
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
ax.set_title('Top 10 Concordancias ARI — Metabólico vs Radiómico',
             fontsize=FONT_TITLE, fontweight='bold', pad=10)
ax.set_ylabel('Algoritmo Metabólico (_M)', fontsize=FONT_SIZE)
ax.set_xlabel('Algoritmo Radiómico (_R)',  fontsize=FONT_SIZE)

plt.tight_layout()
save_plos(fig, 'Fig1_heatmap_top10_M_vs_R.tif')

# ============================================================
# 4️⃣  HEATMAP DE CONTINGENCIA — mejor par M vs R  (Fig 2)
# ============================================================

col_a, col_b = best['col_a'], best['col_b']
df_pair = df_merged[[ID_COLUMN, col_a, col_b]].copy()
df_pair = df_pair[(df_pair[col_a] != -1) & (df_pair[col_b] != -1)]
contingency_matrix = pd.crosstab(df_pair[col_a], df_pair[col_b])

fig, ax = plt.subplots(figsize=(W_HALF, 4.5))
sns.heatmap(
    contingency_matrix,
    annot=True, fmt='d',
    cmap='YlGnBu',
    linewidths=0.4, linecolor='white',
    annot_kws={'size': FONT_SIZE - 1, 'weight': 'bold'},
    cbar_kws={'label': 'N Pacientes', 'shrink': 0.75},
    ax=ax,
)
ax.set_title('Matriz de Contingencia — Metabólico vs Radiómico',
             fontsize=FONT_TITLE, fontweight='bold', pad=8)
ax.set_ylabel(f'Cluster (M)', fontsize=FONT_SIZE)
ax.set_xlabel(f'Cluster (R)', fontsize=FONT_SIZE)
ax.tick_params(labelsize=FONT_SIZE - 1)

plt.tight_layout()
save_plos(fig, 'Fig2_contingency_M_vs_R.tif')

print(f"\n💡 Mejor par M vs R:")
print(f"   Col A (M) : {col_a}")
print(f"   Col B (R) : {col_b}")
print(f"   ARI       : {best['ARI']:.4f} | AMI: {best['AMI']:.4f}")

stacked  = contingency_matrix.stack()
m_best, r_best = stacked.idxmax()
n_patients = stacked.max()
print(f"   Mayor concordancia intra-par: [{m_best}] ↔ [{r_best}]  ({n_patients} pacientes)")

pacientes_core = df_pair[
    (df_pair[col_a] == m_best) &
    (df_pair[col_b] == r_best)
]
pacientes_core[[ID_COLUMN]].to_csv(
    'pacientes_core_metabolico_radiomico.csv', index=False
)
print("✅  'pacientes_core_metabolico_radiomico.csv' generado.")

# ============================================================
# 5️⃣  DISTRIBUCIÓN ARI / AMI — BOXPLOT  (Fig 3)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 4.5))

sns.boxplot(data=df_results, y='ARI', ax=axes[0], width=0.35, color=sns.color_palette('Set2', 1)[0])
sns.stripplot(data=df_results, y='ARI', ax=axes[0], color='#333333', alpha=0.25, size=2.5, jitter=True)
axes[0].set_xlabel('Metabólico vs Radiómico', fontsize=FONT_AXIS, fontweight='bold')
axes[0].set_ylabel('ARI', fontsize=FONT_AXIS, fontweight='bold')
axes[0].set_title('Distribución de ARI', fontsize=FONT_TITLE, fontweight='bold')
axes[0].axhline(0, color='red', lw=0.8, ls='--', alpha=0.6, label='ARI = 0 (azar)')
axes[0].legend(fontsize=FONT_SIZE - 1)

sns.boxplot(data=df_results, y='AMI', ax=axes[1], width=0.35, color=sns.color_palette('Set2', 2)[1])
sns.stripplot(data=df_results, y='AMI', ax=axes[1], color='#333333', alpha=0.25, size=2.5, jitter=True)
axes[1].set_xlabel('Metabólico vs Radiómico', fontsize=FONT_AXIS, fontweight='bold')
axes[1].set_ylabel('AMI', fontsize=FONT_AXIS, fontweight='bold')
axes[1].set_title('Distribución de AMI', fontsize=FONT_TITLE, fontweight='bold')
axes[1].axhline(0, color='red', lw=0.8, ls='--', alpha=0.6, label='AMI = 0 (azar)')
axes[1].legend(fontsize=FONT_SIZE - 1)

plt.tight_layout()
save_plos(fig, 'Fig3_ARI_AMI_boxplot_M_vs_R.tif')

# ============================================================
# 6️⃣  PIE CHART — CATEGORÍAS DE CONCORDANCIA  (Fig 4)
# ============================================================

def categorize_patient(row):
    if row[col_a] == m_best and row[col_b] == r_best:
        return 'Núcleo (Concordante)'
    elif row[col_a] == m_best and row[col_b] != r_best:
        return 'Divergencia R'
    elif row[col_a] != m_best and row[col_b] == r_best:
        return 'Divergencia M'
    else:
        return 'Discrepancia Total'

df_pair['Analysis_Category'] = df_pair.apply(categorize_patient, axis=1)
group_summary = df_pair['Analysis_Category'].value_counts()

print("\n📊 DISTRIBUCIÓN DE LA COHORTE (M vs R):")
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

ax.set_title(
    f'Concordancia vs Divergencia — Metabólico vs Radiómico\n'
    f'(mejor ARI: {best["ARI"]:.3f})',
    fontsize=FONT_TITLE, fontweight='bold', pad=12
)

plt.tight_layout()
save_plos(fig, 'Fig4_cohort_pie_M_vs_R.tif')

divergent_patients = df_pair[df_pair['Analysis_Category'] != 'Núcleo (Concordante)']
divergent_patients.to_csv('pacientes_divergentes_M_vs_R.csv', index=False)
print(f"\n✅ {len(divergent_patients)} pacientes divergentes identificados.")

# ============================================================
# 7️⃣  3D — RELACIÓN ENTRE ALGORITMOS GANADORES  (Fig 5)
# ============================================================

print("\n🚀 3D: relación entre algoritmos ganadores M vs R...")

df_algo = df_pair[[col_a, col_b]].copy()
df_algo[col_a] = df_algo[col_a].astype(int)
df_algo[col_b] = df_algo[col_b].astype(int)

sil_a     = extract_silhouette(col_a)
sil_b     = extract_silhouette(col_b)
sil_a_str = f"{sil_a:.3f}" if sil_a is not None else "N/A"
sil_b_str = f"{sil_b:.3f}" if sil_b is not None else "N/A"

counts = df_algo.groupby([col_a, col_b]).size().reset_index(name='Count')

fig = plt.figure(figsize=(4.5, 4.5))
ax  = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    counts[col_a], counts[col_b], counts['Count'],
    s=counts['Count'] * 2.5,
    c=counts['Count'], cmap='viridis',
    alpha=0.85, edgecolors='white', linewidths=0.4,
)

cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.10)
cb.ax.tick_params(labelsize=FONT_SIZE)
cb.set_label('Número de Pacientes', size=FONT_SIZE)

ax.set_xlabel('Cluster (M)', fontsize=FONT_SIZE, labelpad=8)
ax.set_ylabel('Cluster (R)', fontsize=FONT_SIZE, labelpad=8)
ax.set_zlabel('N Pacientes', fontsize=FONT_SIZE, labelpad=8)
ax.tick_params(labelsize=FONT_SIZE - 1)

ax.set_title(
    f'Distribución 3D — Metabólico vs Radiómico\n'
    f'Sil (M)={sil_a_str}  |  Sil (R)={sil_b_str}  |  ARI={best["ARI"]:.3f}',
    fontsize=FONT_SIZE, fontweight='bold'
)

ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
save_plos(fig, 'Fig5_3D_best_pair_M_vs_R.tif')

# ============================================================
# 8️⃣  3D — DISTRIBUCIÓN POR DOMINIO INDIVIDUAL (M y R)  (Fig 6, 7)
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


print("\n🚀 3D: dominio Metabólico...")
plot_3d_domain(df_pair, ID_COLUMN, col_a, 'Metabólico', sil_a_str,
               'Fig6_3D_metabolico.tif', palette_range=(0.0, 0.6))

print("\n🚀 3D: dominio Radiómico...")
plot_3d_domain(df_pair, ID_COLUMN, col_b, 'Radiómico', sil_b_str,
               'Fig7_3D_radiomico.tif', palette_range=(0.3, 0.9))

# ============================================================
# RESUMEN FINAL
# ============================================================

print("\n" + "=" * 60)
print("  TODAS LAS FIGURAS GENERADAS (Metabólico vs Radiómico)")
print("=" * 60)
print(f"  DPI      : {DPI}")
print(f"  Formato  : TIFF (LZW)")
print(f"  Fuente   : {FONT_FAMILY}, 8–12 pt")
print(f"  Color    : RGB")
print(f"\n  Figuras:")
print(f"    Fig1_heatmap_top10_M_vs_R.tif")
print(f"    Fig2_contingency_M_vs_R.tif")
print(f"    Fig3_ARI_AMI_boxplot_M_vs_R.tif")
print(f"    Fig4_cohort_pie_M_vs_R.tif")
print(f"    Fig5_3D_best_pair_M_vs_R.tif")
print(f"    Fig6_3D_metabolico.tif")
print(f"    Fig7_3D_radiomico.tif")
print(f"\n  Archivos CSV:")
print(f"    pacientes_core_metabolico_radiomico.csv")
print(f"    pacientes_divergentes_M_vs_R.csv")
print("=" * 60)