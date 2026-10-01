
import re, os, warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from collections import defaultdict
import sys
from pathlib import Path
warnings.filterwarnings('ignore')

# ============================================================
# 0️⃣  PLOS ONE STYLE CONFIGURATION 🎨
# ============================================================

PALETTE = {
    'core':      '#2E86AB',
    'divergent': '#E84855',
    'normal':    '#3BB273',
    'tumor':     '#F18F01',
    'neutral':   '#A9A9A9',
    'accent1':   '#7B2D8B',
    'accent2':   '#F7B731',
}
CLUSTER_COLORS = [PALETTE['core'], PALETTE['divergent'], PALETTE['accent1'],
                  PALETTE['accent2'], PALETTE['neutral']]

# Separate palettes so clinical and metabolic plots are visually distinct
CLINICAL_COLORS  = ['#2E86AB', '#1B4F72', '#5DADE2', '#85C1E9', '#154360']   # cool / blue tones
METABOLIC_COLORS = ['#E84855', '#F18F01', '#F7B731', '#7B2D8B', '#C0392B']   # warm / red-orange tones

def _resolve_plos_font(preferred='Arial'):
    available = {f.name for f in fm.fontManager.ttflist}
    candidates = [preferred, 'Helvetica', 'Liberation Sans', 'Nimbus Sans',
                  'Arimo', 'DejaVu Sans']
    for name in candidates:
        if name in available:
            if name != preferred:
                print(f'Warning: "{preferred}" not installed. Using "{name}" instead. '
                      f'Install Arial (e.g. "ttf-mscorefonts-installer" on Linux) '
                      f'or embed the font in the final EPS/TIFF for strict PLOS compliance.')
            return name
    return 'DejaVu Sans'

FONT_FAMILY = _resolve_plos_font('Arial')

DPI = 300
FONT_SIZE  = 9
FONT_TITLE = 10
FONT_AXIS  = 9

W_HALF = 3.27
W_FULL = 6.83
H_MAX  = 8.75

plt.rcParams.update({
    'font.family':      FONT_FAMILY,
    'font.size':        FONT_SIZE,
    'axes.titlesize':   FONT_TITLE,
    'axes.labelsize':   FONT_AXIS,
    'xtick.labelsize':  FONT_SIZE,
    'ytick.labelsize':  FONT_SIZE,
    'legend.fontsize':  FONT_SIZE,
    'figure.dpi':       DPI,
    'savefig.dpi':      DPI,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.grid':          True,
    'grid.alpha':         0.3,
    'axes.linewidth':     0.8,
    'lines.linewidth':    1.2,
    'patch.linewidth':    0.8,
    'savefig.facecolor':  'white',
    'figure.facecolor':   'white',
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

RESULTS_DIR = str(Path(__file__).resolve().parent / 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

def savefig(name, fmt='tiff', also_pdf=False):
    """Save current figure as TIFF (PLOS ONE compliant) with optional PDF copy."""
    fig = plt.gcf()
    w_in, h_in = fig.get_size_inches()
    w_px, h_px = w_in * DPI, h_in * DPI
    if not (789 <= w_px <= 2250):
        print(f'Warning: width of "{name}" = {w_px:.0f}px at {DPI}dpi '
              f'(PLOS range: 789-2250 px). Adjust figsize.')
    if h_px > 2625:
        print(f'Warning: height of "{name}" = {h_px:.0f}px exceeds PLOS max (2625 px). '
              f'Adjust figsize.')

    base = os.path.join(RESULTS_DIR, name)

    if fmt in ('tiff', 'tif'):
        path = base + '.tif'
        plt.savefig(path, format='tiff', dpi=DPI, bbox_inches='tight',
                    facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
    else:
        path = base + f'.{fmt}'
        plt.savefig(path, format=fmt, dpi=DPI, bbox_inches='tight',
                    facecolor='white')

    if also_pdf:
        pdf_path = base + '.pdf'
        plt.savefig(pdf_path, format='pdf', dpi=DPI, bbox_inches='tight',
                    facecolor='white')
        print(f'  Saved: {pdf_path}')

    plt.show()
    size_mb = os.path.getsize(path) / (1024 * 1024)
    print(f'  Saved: {path}  ({size_mb:.2f} MB, {w_px:.0f}x{h_px:.0f}px @{DPI}dpi)')
    if size_mb > 10:
        print(f'  Warning: file > 10 MB, PLOS ONE may reject it. Reduce dpi or simplify figure.')

print('PLOS ONE style configured')
print(f'Font: {FONT_FAMILY} | DPI: {DPI} | Width: {W_HALF}in (1 col) / {W_FULL}in (2 col)')
print(f'Cluster colors: Core={PALETTE["core"]} | Divergent={PALETTE["divergent"]}')

# ============================================================
# 1️⃣ CONFIGURATION AND DATA LOADING 🔑
# ============================================================

ID_COLUMN = 'ModelName'
MIN_VALID_SAMPLES = 50

# Root data directory (3 levels up from this script → Unsupervised_learning_scripts_and_data/Clinical_data_and_models_ids/)
_SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = _SCRIPT_DIR.parents[2] / "Clinical_data_and_models_ids"

# Input file paths
# These point to clustering output CSVs produced by the clinical and metabolic pipeline scripts.
# path1 → output of ML_clinical_data_UMAP_reduced_selection.py (Results_clustering_UMAP_reduced_selection/)
# path2 → output of UMAP_FBA_pFBA_WL1_only.py or PCA_FBA_pFBA_WL1_only.py
_SRC_DIR = _SCRIPT_DIR.parents[0]   # Clustering_and_data_analysis_PYTHON/
path1 = _SRC_DIR / "Clinical_data_analysis" / "ML_models_using_clinical_data" / "Results_clustering_UMAP_reduced_selection" / "all_patients_clustered_no_filter.csv"
path2 = _SRC_DIR / "Metabolic_data_analysis" / "ML_models_using_metabolic_data" / "results_TumorPhenotype_PCA_metrics_updated_nol2" / "PatientClusters_TumorPhenotype_PCA.csv"


def load_and_standardize_clusters(file_path, suffix, id_col='ModelName'):
    """Loads, cleans TCGA IDs and renames cluster columns."""
    try:
        df = pd.read_csv(file_path, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"Failed to read {file_path}. Error: {e}")

    df.columns = df.columns.str.strip()
    if id_col not in df.columns:
        raise ValueError(f"ID column '{id_col}' not found in: {file_path}")

    # Clean TCGA IDs (maintain consistent format)
    df[id_col] = df[id_col].astype(str).str.split('_').str[0].str.split('.').str[0].str.slice(0, 16)

    cluster_cols = [col for col in df.columns if col != id_col]
    rename_mapping = {col: f"{col}_{suffix}" for col in cluster_cols}
    df = df.rename(columns=rename_mapping)

    return df[[id_col] + list(rename_mapping.values())].copy()


# Load and merge
try:
    df1 = load_and_standardize_clusters(path1, suffix='C', id_col=ID_COLUMN)
    df2 = load_and_standardize_clusters(path2, suffix='M', id_col=ID_COLUMN)
    df_merged = df1.merge(df2, on=ID_COLUMN, how='inner')
    print(f"✅ Successful merge. Common patients: {len(df_merged)}")
except Exception as e:
    print(f"❌ ERROR: {e}")
    sys.exit()

# Identify columns
cluster_cols_clinical = [col for col in df_merged.columns if col.endswith('_C')]
cluster_cols_metabolic = [col for col in df_merged.columns if col.endswith('_M')]

# =================================================================
# 2️⃣ CROSS-PIPELINE CONCORDANCE ANALYSIS (ARI / AMI)
# =================================================================
#
# PURPOSE
# -------
# Quantify how well unsupervised clustering partitions obtained from
# *clinical* features agree with those obtained from *metabolic* (flux-
# based) features.  Agreement is measured with two complementary,
# chance-corrected indices:
#
#   • Adjusted Rand Index (ARI)
#     Counts all pairs of samples and checks whether the two
#     partitions agree on placing them in the same or different
#     clusters.  The raw Rand Index is then adjusted for chance
#     using a hypergeometric null model:
#
#         ARI = (RI - E[RI]) / (max(RI) - E[RI])
#
#     Range: −1 to 1.  ARI = 1 means perfect agreement;
#     ARI ≈ 0 means agreement no better than random; ARI < 0
#     means agreement worse than random.
#
#   • Adjusted Mutual Information (AMI)
#     Measures the information shared between the two partitions,
#     normalised and corrected for chance:
#
#         AMI = (MI(U,V) - E[MI]) / (max(H(U),H(V)) - E[MI])
#
#     where MI is mutual information and H is Shannon entropy.
#     Same range and interpretation as ARI.  AMI is more sensitive
#     to differences in cluster *size distributions* than ARI.
#
# METHODOLOGY
# -----------
#   1. Every clinical clustering column (suffix '_C') is compared
#      against every metabolic clustering column (suffix '_M') in
#      an exhaustive pairwise fashion.
#   2. Noise labels (cluster == −1, produced by DBSCAN / HDBSCAN)
#      are excluded from each comparison by treating them as NaN.
#   3. A pair is skipped if fewer than MIN_VALID_SAMPLES (50)
#      patients remain after noise removal, or if either partition
#      contains only one unique label (metrics are undefined).
#   4. Results are collected into a DataFrame; the top 10 pairs by
#      ARI are retained for downstream analysis and visualisation.
#
# REFERENCES
#   - Hubert & Arabie (1985). Comparing partitions. J. Classif.
#   - Vinh, Epps & Bailey (2010). Information Theoretic Measures
#     for Clusterings Comparison. JMLR 11, 2837–2854.
# =================================================================

def calculate_metrics_fixed(df, col1, col2):
    """Compute ARI and AMI between two cluster-label columns.

    Parameters
    ----------
    df : pd.DataFrame
        Merged dataframe containing both clinical and metabolic
        cluster assignments for each patient.
    col1, col2 : str
        Column names holding integer cluster labels.  The value −1
        is treated as noise / unassigned and excluded before any
        metric computation.

    Returns
    -------
    dict or None
        Dictionary with keys 'Clinical_Cluster', 'Metabolic_Cluster',
        'ARI', 'AMI', and 'N_Samples'.  Returns None if the pair
        does not meet the minimum-sample or minimum-cluster criteria.
    """
    # Step 1 – Mask noise labels (−1) produced by density-based algorithms
    s1 = df[col1].replace(-1, np.nan)
    s2 = df[col2].replace(-1, np.nan)

    # Step 2 – Keep only patients with valid labels in BOTH partitions
    valid_idx = s1.dropna().index.intersection(s2.dropna().index)

    # Step 3 – Skip if too few valid samples for a meaningful comparison
    if len(valid_idx) < MIN_VALID_SAMPLES:
        return None

    labels1 = s1.loc[valid_idx].astype(int)
    labels2 = s2.loc[valid_idx].astype(int)

    # Step 4 – Skip degenerate cases (single-cluster partitions)
    if labels1.nunique() < 2 or labels2.nunique() < 2:
        return None

    # Step 5 – Compute chance-corrected concordance metrics
    return {
        'Clinical_Cluster': col1,
        'Metabolic_Cluster': col2,
        'ARI': adjusted_rand_score(labels1, labels2),
        'AMI': adjusted_mutual_info_score(labels1, labels2),
        'N_Samples': len(valid_idx)
    }

# --- Exhaustive pairwise comparison ---
# Every clinical model × every metabolic model → one ARI + AMI value.
print(f"\n🔄 Computing ARI/AMI for {len(cluster_cols_clinical)} clinical "
      f"× {len(cluster_cols_metabolic)} metabolic model pairs...")

results_list = []
for c_col in cluster_cols_clinical:
    for m_col in cluster_cols_metabolic:
        res = calculate_metrics_fixed(df_merged, c_col, m_col)
        if res:
            results_list.append(res)

df_results = pd.DataFrame(results_list)

print(f"   Valid pairs evaluated: {len(df_results)}")
print(f"   ARI range: [{df_results['ARI'].min():.4f}, {df_results['ARI'].max():.4f}]")
print(f"   AMI range: [{df_results['AMI'].min():.4f}, {df_results['AMI'].max():.4f}]")

# Rank by ARI and retain the top 10 most concordant pairs
df_top10 = df_results.sort_values(by='ARI', ascending=False).head(10)

# =================================================================
# 3️⃣ TOP 10 CONCORDANCE HEATMAP
# =================================================================
# Visualise the 10 clinical–metabolic algorithm pairs with the
# highest ARI as a heatmap.  Each cell shows the ARI value for
# one (clinical model, metabolic model) pair.  Higher values
# (brighter in the 'magma' colormap) indicate stronger agreement
# between the two independent clustering pipelines.

pivot_top10 = df_top10.pivot_table(index='Clinical_Cluster', columns='Metabolic_Cluster', values='ARI')

fig, ax = plt.subplots(figsize=(W_FULL, 5.0))
sns.heatmap(pivot_top10, annot=True, fmt=".4f", cmap='magma',
            linewidths=0.5, linecolor='white', ax=ax,
            annot_kws={'fontsize': FONT_SIZE})
ax.set_title('Top 10 Correlations: Clinical vs Metabolic (ARI)')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
savefig("Fig01_heatmap_top10_ARI")

# =================================================================
# 4️⃣ PATIENT SUBGROUP IDENTIFICATION
# =================================================================
# From the top-ranked pair (highest ARI), build a contingency
# matrix to see exactly how patients are distributed across the
# clinical × metabolic cluster combinations.  The cell with the
# largest count identifies the "Core" subgroup — patients for
# whom both pipelines converge on the same biological grouping.

# 1. Select the best-concordance model pair
best_pair = df_top10.iloc[0]
best_c_col = best_pair['Clinical_Cluster']
best_m_col = best_pair['Metabolic_Cluster']

# 2. Exclude noise-labelled patients (cluster == −1) from both models
df_best = df_merged[[ID_COLUMN, best_c_col, best_m_col]].copy()
df_best = df_best[(df_best[best_c_col] != -1) & (df_best[best_m_col] != -1)]

# 3. Contingency matrix: rows = clinical clusters, cols = metabolic clusters.
#    Each cell counts the number of patients assigned to that combination.
contingency_matrix = pd.crosstab(df_best[best_c_col], df_best[best_m_col])

fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))
sns.heatmap(contingency_matrix, annot=True, fmt='d', cmap='YlGnBu',
            linewidths=0.5, linecolor='white', ax=ax,
            annot_kws={'fontsize': FONT_SIZE})
ax.set_title(f'Patient Distribution:\n{best_c_col} vs {best_m_col}')
ax.set_xlabel('Metabolic Clusters (M)')
ax.set_ylabel('Clinical Clusters (C)')
plt.tight_layout()
savefig("Fig02_contingency_matrix")

# 4. Extract the "Core" (the cell with the most patients)
stacked = contingency_matrix.stack()
c_best, m_best = stacked.idxmax()
n_patients = stacked.max()

print(f"\n💡 ANALYSIS RESULT:")
print(f"The highest agreement is between Clinical Cluster '{c_best}' and Metabolic '{m_best}'.")
print(f"This group contains {n_patients} patients out of {len(df_merged)} total.")

# 5. Export IDs
patients_core = df_best[(df_best[best_c_col] == c_best) & (df_best[best_m_col] == m_best)]
patients_core[[ID_COLUMN]].to_csv(
    os.path.join(RESULTS_DIR, "core_patients_correlation_ParetoAndNorms.csv"), index=False)

print(f"✅ File 'core_patients_correlation_ParetoAndNorms.csv' generated.")



# =================================================================
# 6️⃣ DIVERGENCE ANALYSIS (THE "REBEL" PATIENTS) 🚀
# =================================================================

# 1. Create a Category column for all patients
def categorize_patient(row):
    if row[best_c_col] == c_best and row[best_m_col] == m_best:
        return 'Core (Concordant)'
    elif row[best_c_col] == c_best and row[best_m_col] != m_best:
        return 'Metabolic Divergence'  # Same clinical, different metabolism
    elif row[best_c_col] != c_best and row[best_m_col] == m_best:
        return 'Clinical Divergence'   # Same metabolism, different clinical
    else:
        return 'Total Discrepancy'     # Different in both

df_best['Analysis_Category'] = df_best.apply(categorize_patient, axis=1)

# 2. Statistical summary of groups
group_summary = df_best['Analysis_Category'].value_counts()
print("\n📊 COHORT DISTRIBUTION:")
print(group_summary)

# 3. Visualization of the cohort composition
pie_colors = [PALETTE['core'], PALETTE['divergent'], PALETTE['accent1'], PALETTE['neutral']]
# Ensure we have enough colors for all categories
while len(pie_colors) < len(group_summary):
    pie_colors.append(PALETTE['accent2'])
pie_colors = pie_colors[:len(group_summary)]

fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))
wedges, texts, autotexts = ax.pie(
    group_summary.values, labels=group_summary.index,
    autopct='%1.1f%%', startangle=140, colors=pie_colors,
    explode=[0.03] * len(group_summary),
    textprops={'fontsize': FONT_SIZE},
    pctdistance=0.75)
for at in autotexts:
    at.set_fontsize(FONT_SIZE - 1)
ax.set_title('Cohort Composition:\nConcordance vs Divergence')
plt.tight_layout()
savefig("Fig03_cohort_composition_pie")

# 4. Save the "Special" patients (Divergent)
# These are the ones who might have a different treatment response than expected
patients_divergent = df_best[df_best['Analysis_Category'] != 'Core (Concordant)']
patients_divergent.to_csv(
    os.path.join(RESULTS_DIR, "divergent_patients_study_pyn.csv"), index=False)

print(f"\n✅ {len(patients_divergent)} divergent patients identified.")
print("File 'divergent_patients_study_pyn.csv' ready for in-depth clinical analysis.")

# =================================================================
# 7️⃣ 3D VISUALIZATION – CONVERGENCE BETWEEN THE TWO BEST ALGORITHMS
# =================================================================

print("\n🚀 3D visualization of the relationship between winning algorithms...")

df_algo = df_best[[best_c_col, best_m_col]].copy()

# Convert to integers
df_algo[best_c_col] = df_algo[best_c_col].astype(int)
df_algo[best_m_col] = df_algo[best_m_col].astype(int)

# Count patients per combination
counts = df_algo.groupby([best_c_col, best_m_col]).size().reset_index(name='Count')

fig = plt.figure(figsize=(W_FULL, 5.5))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    counts[best_c_col],
    counts[best_m_col],
    counts['Count'],
    s=counts['Count'] * 2,
    c=counts['Count'],
    cmap='viridis',
    edgecolors='k', linewidths=0.3,
    alpha=0.85
)

ax.set_xlabel(f'Clinical Cluster\n({best_c_col})')
ax.set_ylabel(f'Metabolic Cluster\n({best_m_col})')
ax.set_zlabel('Number of Patients')
ax.set_title('3D Patient Distribution Between Algorithms')
fig.colorbar(sc, ax=ax, shrink=0.5, label='Patient Count')

plt.tight_layout()
savefig("Fig04_3D_algorithm_convergence", also_pdf=True)

print("✅ Figure saved: Fig04_3D_algorithm_convergence")

# =================================================================
# 7️⃣b 3D VISUALIZATION – CLINICAL ALGORITHM (COLORS BY CLUSTER)
# =================================================================

print("\n🚀 Visualizing clinical algorithm...")

df_clin = df_best[[ID_COLUMN, best_c_col]].copy()
df_clin = df_clin.sort_values(by=best_c_col).reset_index(drop=True)

df_clin['X'] = np.arange(len(df_clin))
df_clin['Y'] = df_clin[best_c_col].astype(int)
df_clin['Z'] = np.random.normal(0, 0.2, len(df_clin))

clusters_clin = np.sort(df_clin['Y'].unique())
n_clin = len(clusters_clin)
clin_colors = (CLINICAL_COLORS * ((n_clin // len(CLINICAL_COLORS)) + 1))[:n_clin]

fig = plt.figure(figsize=(W_FULL, 5.5))
ax = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_clin, clin_colors):
    subset = df_clin[df_clin['Y'] == c]
    ax.scatter(
        subset['X'], subset['Y'], subset['Z'],
        label=f'Cluster {c}',
        color=color,
        alpha=0.85,
        s=20, edgecolors='k', linewidths=0.2
    )

ax.set_title(f"Clinical Algorithm\n3D Cluster Distribution ({best_c_col})")
ax.set_xlabel('Patients (sorted)')
ax.set_ylabel('Clinical Cluster')
ax.set_zlabel('Jitter (visualization)')
ax.legend(title="Clusters", loc='best', fontsize=FONT_SIZE - 1,
          title_fontsize=FONT_SIZE)

plt.tight_layout()
savefig("Fig05_3D_clinical_clusters", also_pdf=True)

# =================================================================
# 8️⃣ 3D VISUALIZATION – METABOLIC ALGORITHM (COLORS BY CLUSTER)
# =================================================================

print("\n🚀 Visualizing metabolic algorithm...")

df_met = df_best[[ID_COLUMN, best_m_col]].copy()
df_met = df_met.sort_values(by=best_m_col).reset_index(drop=True)

df_met['X'] = np.arange(len(df_met))
df_met['Y'] = df_met[best_m_col].astype(int)
df_met['Z'] = np.random.normal(0, 0.2, len(df_met))

clusters_met = np.sort(df_met['Y'].unique())
n_met = len(clusters_met)
met_colors = (METABOLIC_COLORS * ((n_met // len(METABOLIC_COLORS)) + 1))[:n_met]

fig = plt.figure(figsize=(W_FULL, 5.5))
ax = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_met, met_colors):
    subset = df_met[df_met['Y'] == c]
    ax.scatter(
        subset['X'], subset['Y'], subset['Z'],
        label=f'Cluster {c}',
        color=color,
        alpha=0.85,
        s=20, edgecolors='k', linewidths=0.2
    )

ax.set_title(f"Metabolic Algorithm\n3D Cluster Distribution ({best_m_col})")
ax.set_xlabel('Patients (sorted)')
ax.set_ylabel('Metabolic Cluster')
ax.set_zlabel('Jitter (visualization)')
ax.legend(title="Clusters", loc='best', fontsize=FONT_SIZE - 1,
          title_fontsize=FONT_SIZE)

plt.tight_layout()
savefig("Fig06_3D_metabolic_clusters", also_pdf=True)

print("✅ All figures generated with PLOS ONE formatting.")
