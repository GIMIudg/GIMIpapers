
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from collections import defaultdict
import sys
from pathlib import Path

# ============================================================
# 1️⃣ CONFIGURATION AND DATA LOADING 🔑
# ============================================================

ID_COLUMN = 'ModelName'
MIN_VALID_SAMPLES = 50

# Root data directory (3 levels up from this script → Unsupervised_learning_scripts_and_data/Clinical_data_and_models_ids/)
_SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = _SCRIPT_DIR.parents[2] / "Clinical_data_and_models_ids"

# Input file paths
# Update these paths to point to the clustering output CSVs from the clinical and metabolic scripts.
path1 = _SCRIPT_DIR.parents[2] / "src" / "Clustering_and_data_analysis_PYTHON" / "Clinical_data_analysis" / "ML_models_using_clinical_data" / "Results_clustering_UMAP_seleccion_variables" / "pacientes_clusterizados_todos_sinfiltro.csv"
path2 = _SCRIPT_DIR.parents[2] / "src" / "Clustering_and_data_analysis_PYTHON" / "Metabolic_data_analysis" / "ML_models_using_metabolic_data" / "resultados_TumorPhenotype_PCA_metrics" / "PatientClusters_TumorPhenotype_PCA.csv"


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
# 2️⃣ METRICS CALCULATION (ARI/AMI) 📊
# =================================================================

def calculate_metrics_fixed(df, col1, col2):
    s1 = df[col1].replace(-1, np.nan)
    s2 = df[col2].replace(-1, np.nan)
    valid_idx = s1.dropna().index.intersection(s2.dropna().index)

    if len(valid_idx) < MIN_VALID_SAMPLES: return None

    labels1 = s1.loc[valid_idx].astype(int)
    labels2 = s2.loc[valid_idx].astype(int)

    if labels1.nunique() < 2 or labels2.nunique() < 2: return None

    return {
        'Clinical_Cluster': col1,
        'Metabolic_Cluster': col2,
        'ARI': adjusted_rand_score(labels1, labels2),
        'AMI': adjusted_mutual_info_score(labels1, labels2),
        'N_Samples': len(valid_idx)
    }

results_list = []
for c_col in cluster_cols_clinical:
    for m_col in cluster_cols_metabolic:
        res = calculate_metrics_fixed(df_merged, c_col, m_col)
        if res: results_list.append(res)

df_results = pd.DataFrame(results_list)
df_top10 = df_results.sort_values(by='ARI', ascending=False).head(10)

# =================================================================
# 3️⃣ TOP 10 VISUALIZATION 📈
# =================================================================

pivot_top10 = df_top10.pivot_table(index='Clinical_Cluster', columns='Metabolic_Cluster', values='ARI')

plt.figure(figsize=(12, 8))
sns.heatmap(pivot_top10, annot=True, fmt=".4f", cmap='magma')
plt.title('Top 10 Correlations: Clinical vs Metabolic (ARI)', fontsize=15)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("heatmap_top10_oncologia_paretoynormas.png")
plt.show()

# =================================================================
# 4️⃣ PATIENT SUBGROUP IDENTIFICATION 🕵️‍♂️
# =================================================================

# 1. Get the best model pair
best_pair = df_top10.iloc[0]
best_c_col = best_pair['Clinical_Cluster']
best_m_col = best_pair['Metabolic_Cluster']

# 2. Filter patients without noise in these two models
df_best = df_merged[[ID_COLUMN, best_c_col, best_m_col]].copy()
df_best = df_best[(df_best[best_c_col] != -1) & (df_best[best_m_col] != -1)]

# 3. Contingency Matrix (Confusion Matrix)
contingency_matrix = pd.crosstab(df_best[best_c_col], df_best[best_m_col])

#[Image of a confusion matrix heatmap]

plt.figure(figsize=(8, 6))
sns.heatmap(contingency_matrix, annot=True, fmt='d', cmap='YlGnBu')
plt.title(f'Patient Distribution:\n{best_c_col} vs {best_m_col}')
plt.xlabel('Metabolic Clusters (M)')
plt.ylabel('Clinical Clusters (C)')
plt.show()

# 4. Extract the "Core" (the cell with the most patients)
stacked = contingency_matrix.stack()
c_best, m_best = stacked.idxmax()
n_patients = stacked.max()

print(f"\n💡 ANALYSIS RESULT:")
print(f"The highest agreement is between Clinical Cluster '{c_best}' and Metabolic '{m_best}'.")
print(f"This group contains {n_patients} patients out of {len(df_merged)} total.")

# 5. Export IDs
patients_core = df_best[(df_best[best_c_col] == c_best) & (df_best[best_m_col] == m_best)]
patients_core[[ID_COLUMN]].to_csv("core_patients_correlation_ParetoAndNorms.csv", index=False)

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
plt.figure(figsize=(10, 6))
colors = sns.color_palette("pastel")[0:4]
group_summary.plot(kind='pie', autopct='%1.1f%%', startangle=140, colors=colors, explode=(0.1, 0.1, 0.1, 0.1))
plt.title('Cohort Composition: Concordance vs Divergence')
plt.ylabel('')
plt.savefig("composicion_grupos_tesis.png")
plt.show()

# 4. Save the "Special" patients (Divergent)
# These are the ones who might have a different treatment response than expected
patients_divergent = df_best[df_best['Analysis_Category'] != 'Core (Concordant)']
patients_divergent.to_csv("divergent_patients_study_pyn.csv", index=False)

print(f"\n✅ {len(patients_divergent)} divergent patients identified.")
print("File 'divergent_patients_study_pyn.csv' ready for in-depth clinical analysis.")

# =================================================================
# 7️⃣ 3D VISUALIZATION – CONVERGENCE BETWEEN THE TWO BEST ALGORITHMS
# =================================================================

from mpl_toolkits.mplot3d import Axes3D

print("\n🚀 3D visualization of the relationship between winning algorithms...")

df_algo = df_best[[best_c_col, best_m_col]].copy()

# Convert to integers
df_algo[best_c_col] = df_algo[best_c_col].astype(int)
df_algo[best_m_col] = df_algo[best_m_col].astype(int)

# Count patients per combination
counts = df_algo.groupby([best_c_col, best_m_col]).size().reset_index(name='Count')

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(
    counts[best_c_col],
    counts[best_m_col],
    counts['Count'],
    s=counts['Count'] * 2,
    alpha=0.8
)

ax.set_xlabel(f'Clinical Cluster ({best_c_col})')
ax.set_ylabel(f'Metabolic Cluster ({best_m_col})')
ax.set_zlabel('Number of Patients')

ax.set_title('3D Patient Distribution Between Algorithms')

plt.tight_layout()
plt.savefig("relacion_3D_algoritmos_ganadores.png")
plt.show()

print("✅ Figure saved: relacion_3D_algoritmos_ganadores.png")

# =================================================================
# 7️⃣ 3D VISUALIZATION – CLINICAL ALGORITHM (COLORS BY CLUSTER)
# =================================================================

print("\n🚀 Visualizing clinical algorithm...")

df_clin = df_best[[ID_COLUMN, best_c_col]].copy()
df_clin = df_clin.sort_values(by=best_c_col).reset_index(drop=True)

df_clin['X'] = np.arange(len(df_clin))
df_clin['Y'] = df_clin[best_c_col].astype(int)
df_clin['Z'] = np.random.normal(0, 0.2, len(df_clin))

clusters_clin = np.sort(df_clin['Y'].unique())
colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_clin)))

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_clin, colors):
    subset = df_clin[df_clin['Y'] == c]
    ax.scatter(
        subset['X'], subset['Y'], subset['Z'],
        label=f'Cluster {c}',
        color=color,
        alpha=0.85,
        s=45
    )

ax.set_title(
    f"Clinical Algorithm\n3D Cluster Distribution\n({best_c_col})",
    fontsize=14,
    fontweight='bold'
)

ax.set_xlabel('Patients (sorted)')
ax.set_ylabel('Clinical Cluster')
ax.set_zlabel('Jitter (visualization)')
ax.legend(title="Clusters", loc='best')

plt.tight_layout()
plt.savefig("clusters_3D_clinico.png", dpi=300)
plt.show()

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
colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_met)))

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

for c, color in zip(clusters_met, colors):
    subset = df_met[df_met['Y'] == c]
    ax.scatter(
        subset['X'], subset['Y'], subset['Z'],
        label=f'Cluster {c}',
        color=color,
        alpha=0.85,
        s=45
    )

ax.set_title(
    f"Metabolic Algorithm\n3D Cluster Distribution\n({best_m_col})",
    fontsize=14,
    fontweight='bold'
)

ax.set_xlabel('Patients (sorted)')
ax.set_ylabel('Metabolic Cluster')
ax.set_zlabel('Jitter (visualization)')
ax.legend(title="Clusters", loc='best')

plt.tight_layout()
plt.savefig("clusters_3D_metabolico.png", dpi=300)
plt.show()

print("✅ Figures generated with complete titles and legends.")
