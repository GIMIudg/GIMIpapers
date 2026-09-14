# ============================================================
# 🚀 Pipeline: Clinical & Demographic Clustering (PCA Manifold Analysis)
# Source Data: TCGA-BRCA Clinical, Survival, and Molecular Metadata
# Method: Deterministic Multi-Configuration PCA & Benchmark Clustering Suite
# Quality Metrics: Silhouette, Davies-Bouldin, and Calinski-Harabasz Scores
# Integration: Master Multi-Omics Harmonization & Multi-Tier ID Merging
# ============================================================

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

import sklearn
from sklearn.model_selection import ParameterGrid
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, OneHotEncoder, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.cluster import (
    KMeans, AgglomerativeClustering, Birch,
    DBSCAN, MeanShift, AffinityPropagation
)
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.metrics import (
    silhouette_score, davies_bouldin_score, calinski_harabasz_score
)
from sklearn.impute import SimpleImputer

warnings.filterwarnings("ignore")

# Conditional import for optional HDBSCAN package
try:
    import hdbscan
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False
    print("⚠️  Warning: Optional module 'hdbscan' not found. Omitting HDBSCAN.")

# ============================================================
# ⚙️ CONFIGURATION & PATH SETUP
# ============================================================

# Obtiene la raíz del proyecto (directorio donde reside este script o su carpeta superior)
BASE_DIR = Path(__file__).resolve().parent

# Construcción de rutas relativas e independientes del sistema operativo
DATA_DIR = BASE_DIR / "Clinical_data_and_models_ids"
PATH_FEATURES = DATA_DIR / "Clinical_data_and_models_ids" / "FeatureMatrix_TumorPhenotype_norm2agregado.csv"
PATH_CLINICAL = DATA_DIR / "TCGA-BRCA.clinical.tsv"
PATH_SURVIVAL = DATA_DIR / "TCGA-BRCA.survival.tsv.gz"

OUT_DIR = BASE_DIR / "results_TumorPhenotype_PCA"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PATH_BASE = "/Users/eduardoruiz/Documents/MCBCI/MCBCI2/Sistemas metabólicos/Proyecto_Tesis/Datos_actual/"
COL_SAMPLE_ID = 'sample'
COL_SAMPLE_TYPE = 'sample_type.samples'
OUTPUT_DIR = "Results_clustering_PCA_seleccion_reducida_conmerge-NUEVOSDATOS"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Random seeds tested across all stochastic clustering algorithms
SEEDS_TO_TEST = [42, 123, 100]
RANDOM_SEED = SEEDS_TO_TEST[0]
np.random.seed(RANDOM_SEED)

# ============================================================
# 2. DATA INGESTION & SAMPLE STANDARDIZATION
# ============================================================

def load_data(filename: str) -> pd.DataFrame:
    """Loads CSV, TSV, GZ, or Excel files from the base directory path."""
    full_path = os.path.join(PATH_BASE, filename)
    try:
        if filename.endswith(('.xlsx', '.xls')):
            return pd.read_excel(full_path)
        elif filename.endswith('.gz'):
            return pd.read_csv(full_path, sep="\t", compression='gzip')
        else:
            return pd.read_csv(full_path, sep="\t")
    except FileNotFoundError:
        print(f"❌ ERROR: File '{filename}' not found at {PATH_BASE}")
        return pd.DataFrame()
    except Exception as e:
        print(f"❌ ERROR while reading '{filename}': {e}")
        return pd.DataFrame()

df_clinical     = load_data("TCGA-BRCA.clinical.tsv")
df_survival     = load_data("TCGA-BRCA.survival.tsv.gz")
df_metadata_raw = load_data("MetaData.xlsx")
df_model_names  = load_data("Model's_ids.txt")

for name, df_check in [("Clinical", df_clinical), ("Survival", df_survival), ("Model names", df_model_names)]:
    if df_check.empty:
        raise ValueError(f"❌ Essential dataset '{name}' could not be loaded. Verify file path.")

print(f"Datasets loaded: Clinical ({len(df_clinical)}), Survival ({len(df_survival)}), Metadata ({len(df_metadata_raw)})")


def extract_sample_id(filename: str) -> str:
    """Standardizes TCGA barcodes to 16-character sample identifiers (TCGA-XX-XXXX-XX)."""
    match = re.search(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}[A-Z0-9]?)', filename)
    if match:
        return match.group(0)[:16]
    return filename.split('_')[0].strip()[:16]


# Build index cohort based on model identifiers
lista_modelos_unicos = df_model_names.iloc[:, 0].dropna().astype(str).tolist()
model_sample_ids = [extract_sample_id(m) for m in lista_modelos_unicos]
df_modelos_base = pd.DataFrame({
    COL_SAMPLE_ID: model_sample_ids,
    'modelo_path_completo': lista_modelos_unicos,
})

# ============================================================
# 3. CLINICAL, SURVIVAL & METADATA HARMONIZATION
# ============================================================

metadata_cols_to_keep = [
    'Menopausal Status', 'Cancer Type', 'ER', 'PR', 'HER2', 
    'Subtype', 'Genetic Ancestry', 'Survival Status', 
    'Survival Time (years)', 'Sex'
]

if not df_metadata_raw.empty and 'hidden' in df_metadata_raw.columns:
    df_metadata_clean = df_metadata_raw.copy()
    df_metadata_clean['temp_id'] = (
        df_metadata_clean['hidden'].astype(str)
        .str.replace(r'\.', '-', regex=True)
        .str.slice(0, 16)
    )
    df_metadata_clean.rename(columns={'temp_id': COL_SAMPLE_ID}, inplace=True)
    df_metadata_clean = df_metadata_clean.drop_duplicates(subset=[COL_SAMPLE_ID], keep='first')
    available_meta_cols = [c for c in metadata_cols_to_keep if c in df_metadata_clean.columns]
    df_metadata_clean = df_metadata_clean[[COL_SAMPLE_ID] + available_meta_cols].copy()

    # Create 12-character patient identifier for fallback imputation
    df_metadata_clean['patient_id'] = df_metadata_clean[COL_SAMPLE_ID].str.slice(0, 12)

    # Reference table indexed strictly by Patient ID
    df_metadata_by_patient = (
        df_metadata_clean
        .drop_duplicates(subset=['patient_id'], keep='first')
        [['patient_id'] + available_meta_cols]
        .copy()
    )

    print(f"✅ Metadata ingested: {len(df_metadata_clean)} samples | "
          f"{df_metadata_clean['patient_id'].nunique()} unique patients")
else:
    print("❌ Metadata not available. Proceeding without molecular subtyping metadata.")
    df_metadata_clean      = pd.DataFrame()
    df_metadata_by_patient = pd.DataFrame()
    available_meta_cols    = []

if 'sample' in df_clinical.columns:
    df_clinical['sample'] = df_clinical['sample'].apply(lambda x: extract_sample_id(str(x)))
if 'sample' in df_survival.columns:
    df_survival['sample'] = df_survival['sample'].apply(lambda x: extract_sample_id(str(x)))

df_survival_clean = (
    df_survival[['sample', 'OS.time', 'OS']]
    .drop_duplicates(subset=['sample'], keep='first')
    .rename(columns={'sample': COL_SAMPLE_ID})
)

df_merged_clinical = pd.merge(
    df_modelos_base,
    df_clinical.drop(columns=['id', 'case_id'], errors='ignore'),
    on=COL_SAMPLE_ID, how='left'
)
df_final = pd.merge(df_merged_clinical, df_survival_clean, on=COL_SAMPLE_ID, how='left')

# MERGE LEVEL 1: Exact 16-character sample ID match
if not df_metadata_clean.empty:
    cols_overlap = [col for col in available_meta_cols if col in df_final.columns]
    df_final = pd.merge(
        df_final.drop(columns=cols_overlap, errors='ignore'),
        df_metadata_clean.drop(columns=['patient_id']),
        on=COL_SAMPLE_ID, how='left'
    )

# MERGE LEVEL 2: Patient ID (12-character) fallback for missing clinical variables
if not df_metadata_by_patient.empty:
    df_final['patient_id'] = df_final[COL_SAMPLE_ID].str.slice(0, 12)
    missing_mask = df_final[available_meta_cols].isna().any(axis=1)
    n_missing_before = missing_mask.sum()

    df_fill = pd.merge(
        df_final.loc[missing_mask, ['patient_id']],
        df_metadata_by_patient,
        on='patient_id', how='left'
    )
    df_fill.index = df_final.index[missing_mask]

    for col in available_meta_cols:
        if col in df_fill.columns:
            df_final.loc[missing_mask, col] = df_final.loc[missing_mask, col].fillna(df_fill[col])

    df_final = df_final.drop(columns=['patient_id'])
    n_missing_after = df_final[available_meta_cols].isna().any(axis=1).sum()

    print("\n🔁 Patient-Level Metadata Fallback Resolution:")
    print(f"    Rows with incomplete metadata BEFORE : {n_missing_before}")
    print(f"    Rows with incomplete metadata AFTER  : {n_missing_after}")
    print(f"    Successfully recovered patient rows  : {n_missing_before - n_missing_after}")

    print("\n📊 Missing count by metadata column:")
    for col in available_meta_cols:
        if col in df_final.columns:
            n_na = df_final[col].isna().sum()
            print(f"    {col:<25} → {n_na} missing ({n_na/len(df_final)*100:.1f}%)")

df_filtered = df_final.copy()

# ============================================================
# 4. CLINICAL DESCRIPTOR SELECTION & FEATURE ENGINEERING
# ============================================================
descriptores_iniciales = [
    'ajcc_pathologic_stage.diagnoses',
    'ajcc_pathologic_t.diagnoses',
    'ajcc_pathologic_n.diagnoses',
    'ajcc_pathologic_m.diagnoses',
    'morphology.diagnoses',
    'primary_diagnosis.diagnoses',
    'treatment_type.treatments.diagnoses',
    'treatment_or_therapy.treatments.diagnoses',
    'prior_treatment.diagnoses',
    'sample_type.samples',
    'tissue_type.samples',
    'age_at_diagnosis.diagnoses',
]

descriptores_alta_res = [
    'Menopausal Status', 'Cancer Type', 'ER', 'PR', 'HER2', 'Subtype'
]

descriptores_finales = descriptores_iniciales + descriptores_alta_res
final_cols = [c for c in descriptores_finales if c in df_filtered.columns]

if not final_cols:
    cols_to_exclude = ['submitter_id', 'sample', 'modelo_path_completo']
    final_cols = [col for col in df_filtered.columns if col not in cols_to_exclude]

if not final_cols:
    raise ValueError("❌ No valid clinical descriptor columns found in DataFrame.")

id_cols_to_keep = [c for c in ['submitter_id', 'sample', 'modelo_path_completo'] if c in df_filtered.columns]
df_aug = df_filtered[final_cols + id_cols_to_keep].copy()

# Binary flag: Prior therapeutic intervention
def prior_treatment_flag(r: pd.Series) -> int:
    cols = ['prior_malignancy.diagnoses', 'prior_treatment.diagnoses', 'progression_or_recurrence.diagnoses']
    for c in cols:
        if c in r.index and pd.notna(r[c]):
            val = str(r[c]).lower().strip()
            if val in ['yes', 'true', 'had prior treatment', 'recurrence', 'progression']:
                return 1
    return 0

df_aug['Prior_Treatment_Flag'] = df_aug.apply(prior_treatment_flag, axis=1)

# Binary flag: Metastatic dissemination or advanced nodal burden
if all(c in df_aug.columns for c in ['ajcc_pathologic_m.diagnoses', 'ajcc_pathologic_n.diagnoses']):
    df_aug['Metastasis_Flag'] = df_aug.apply(
        lambda r: 1 if ('m1' in str(r['ajcc_pathologic_m.diagnoses']).lower() or
                        'n2' in str(r['ajcc_pathologic_n.diagnoses']).lower() or
                        'n3' in str(r['ajcc_pathologic_n.diagnoses']).lower()) else 0,
        axis=1
    )
else:
    df_aug['Metastasis_Flag'] = 0

# Molecular Breast Cancer Classification based on receptor status
receptor_cols = ['ER', 'PR', 'HER2']
if all(c in df_aug.columns for c in receptor_cols):
    def classify_molecular_subtype(row: pd.Series) -> str:
        er   = str(row['ER']).lower().strip()   if pd.notna(row['ER'])   else 'na'
        pr   = str(row['PR']).lower().strip()   if pd.notna(row['PR'])   else 'na'
        her2 = str(row['HER2']).lower().strip() if pd.notna(row['HER2']) else 'na'
        
        is_er_pos   = er   in ['positive', '+']
        is_pr_pos   = pr   in ['positive', '+']
        is_her2_pos = her2 in ['positive', '+', 'amplified', 'equivocal']
        
        if not is_er_pos and not is_pr_pos and not is_her2_pos:
            return 'Triple_Negative'
        elif is_her2_pos and (is_er_pos or is_pr_pos):
            return 'Luminal_HER2+'
        elif is_her2_pos and not (is_er_pos or is_pr_pos):
            return 'HER2_Enriched'
        elif is_er_pos or is_pr_pos:
            return 'Luminal_HR+'
        else:
            return 'Unknown'
            
    df_aug['ER_PR_HER2_Combo'] = df_aug.apply(classify_molecular_subtype, axis=1)
else:
    df_aug['ER_PR_HER2_Combo'] = 'Unknown'

# ============================================================
# 5. PREPROCESSING & ENCODING FOR MACHINE LEARNING
# ============================================================

label_cols_candidates = ['Molecular_Subtype', 'ER_PR_HER2_Combo', 'Subtype']
label_encoders = {}

for col in label_cols_candidates:
    if col in df_aug.columns and df_aug[col].dtype in ['object', 'category']:
        le = LabelEncoder()
        df_aug[f"{col}_encoded"] = le.fit_transform(df_aug[col].fillna('Unknown').astype(str))
        label_encoders[col] = le

id_and_meta_cols = set(id_cols_to_keep + ['submitter_id', 'sample', 'modelo_path_completo'])
numeric_cols = [
    c for c in df_aug.select_dtypes(include=['int64', 'float64', 'float32', 'int32', 'uint8']).columns
    if c not in id_and_meta_cols
]
label_orig_cols = [
    col for col in label_cols_candidates
    if col in df_aug.columns and f"{col}_encoded" in df_aug.columns
]
categorical_cols = [
    c for c in df_aug.select_dtypes(include=['object', 'category']).columns
    if c not in id_and_meta_cols and c not in label_orig_cols
]

# Impute numeric descriptors with zero
cols_to_impute_zero = [c for c in numeric_cols if df_aug[c].isna().any()]
if cols_to_impute_zero:
    imputer_num = SimpleImputer(strategy='constant', fill_value=0)
    df_aug.loc[:, cols_to_impute_zero] = imputer_num.fit_transform(df_aug.loc[:, cols_to_impute_zero])

# Fill missing categorical descriptors with explicit string
for col in categorical_cols:
    df_aug[col] = df_aug[col].fillna('Missing').astype(str)

sklearn_version = tuple(int(x) for x in sklearn.__version__.split(".")[:2])
ohe_kwargs = {'handle_unknown': 'ignore', 'sparse_output': False} if sklearn_version >= (1, 2) \
             else {'handle_unknown': 'ignore', 'sparse': False}

transformers = []
if numeric_cols:
    transformers.append(('num', StandardScaler(), numeric_cols))
if categorical_cols:
    transformers.append(('cat', OneHotEncoder(**ohe_kwargs), categorical_cols))

if not transformers:
    raise ValueError("❌ No valid numeric or categorical features available for preprocessing.")

preprocessor = ColumnTransformer(transformers, remainder='drop')
X_scaled = preprocessor.fit_transform(df_aug)

final_columns = []
if 'num' in preprocessor.named_transformers_:
    final_columns.extend(numeric_cols)
if 'cat' in preprocessor.named_transformers_ and preprocessor.named_transformers_['cat'] is not None:
    final_columns.extend(preprocessor.named_transformers_['cat'].get_feature_names_out(categorical_cols))

final_columns = np.array(final_columns)
print(f"✅ Preprocessing pipeline complete: {len(final_columns)} features encoded for ML.")

# ============================================================
# 6. CLUSTERING BENCHMARKING & EVALUATION ENGINE
# ============================================================

K_VALORES = range(2, 10)

param_grids = {
    "KMeans": {'n_clusters': K_VALORES},
    "Agglomerative": {'n_clusters': K_VALORES, 'linkage': ['ward', 'average', 'complete']},
    "Birch": {'n_clusters': K_VALORES},
    "GMM": {'n_components': K_VALORES},
    "BayesianGaussianMixture": {'n_components': K_VALORES},
    "DBSCAN": {'eps': [0.5, 1.0, 1.5], 'min_samples': [5, 10, 20]},
    "MeanShift": {'bandwidth': [None]},
    "AffinityPropagation": {'damping': [0.5, 0.9]},
}

if HDBSCAN_AVAILABLE:
    param_grids["HDBSCAN"] = {'min_cluster_size': [5, 10, 20], 'min_samples': [None, 5, 10]}

alg_classes = {
    'KMeans': KMeans,
    'Agglomerative': AgglomerativeClustering,
    'Birch': Birch,
    'GMM': GaussianMixture,
    'DBSCAN': DBSCAN,
    'MeanShift': MeanShift,
    'AffinityPropagation': AffinityPropagation,
    'BayesianGaussianMixture': BayesianGaussianMixture,
}
if HDBSCAN_AVAILABLE:
    alg_classes['HDBSCAN'] = hdbscan.HDBSCAN

alg_classes = {k: v for k, v in alg_classes.items() if v is not None}
param_grids = {k: v for k, v in param_grids.items() if v}


def clustering_evaluation_table(X: np.ndarray, algorithms: dict, param_grids: dict, seeds_list: list):
    """
    Evaluates clustering algorithms across hyperparameter grids and seeds.
    Ranks solutions by Silhouette Score while imposing an upper bound on outlier/noise ratio.
    """
    results = {}
    summary_rows = []
    MAX_NOISE_PCT = 10.0

    for alg_name, alg_class in algorithms.items():
        best_score = -1
        best_params = best_labels = best_seed = None
        best_n_clusters = best_noise_pct = 0

        current_grid_dict = param_grids.get(alg_name, [{}])
        param_grid = (
            ParameterGrid(current_grid_dict)
            if isinstance(current_grid_dict, dict)
            else ParameterGrid([current_grid_dict])
        )

        for current_seed in seeds_list:
            if alg_name not in ['GMM', 'KMeans', 'BayesianGaussianMixture', 'Birch'] and current_seed != seeds_list[0]:
                continue
            for param_val in param_grid:
                try:
                    if alg_name in ['GMM', 'KMeans', 'BayesianGaussianMixture', 'Birch']:
                        model = alg_class(**param_val, random_state=current_seed)
                    else:
                        model = alg_class(**param_val)

                    if alg_name == 'MeanShift':
                        model.fit(X)
                        labels = model.labels_
                    else:
                        labels = model.fit_predict(X)

                    valid_mask = labels != -1
                    n_clusters = len(set(labels[valid_mask]))
                    noise_pct = (labels == -1).sum() / len(labels) * 100 if -1 in labels else 0.0

                    if noise_pct > MAX_NOISE_PCT:
                        continue

                    score = (
                        silhouette_score(X[valid_mask], labels[valid_mask])
                        if n_clusters > 1 and valid_mask.sum() > 1
                        else -1
                    )

                    if score > best_score:
                        best_score, best_params, best_labels = score, param_val, labels
                        best_n_clusters, best_noise_pct, best_seed = n_clusters, noise_pct, current_seed

                except Exception:
                    continue

        if best_labels is not None:
            try:
                mask = best_labels != -1
                db_score = davies_bouldin_score(X[mask], best_labels[mask]) if best_n_clusters > 1 else np.nan
                ch_score = calinski_harabasz_score(X[mask], best_labels[mask]) if best_n_clusters > 1 else np.nan
            except Exception:
                db_score = ch_score = np.nan

            results[alg_name] = {
                "best_score": best_score, "best_params": best_params, "best_labels": best_labels,
                "n_clusters": best_n_clusters, "noise_pct": best_noise_pct,
                "best_seed": best_seed, "db_score": db_score, "ch_score": ch_score
            }
            summary_rows.append({
                "Algorithm": alg_name, "Silhouette Score": best_score,
                "Davies-Bouldin Score": db_score, "Calinski-Harabasz Score": ch_score,
                "Best Params": best_params, "Clusters": best_n_clusters,
                "Noise %": best_noise_pct, "Best Seed": best_seed,
                "Unique Labels": np.unique(best_labels)
            })

    summary_df = pd.DataFrame(summary_rows).sort_values(by="Silhouette Score", ascending=False)
    print("\n=== Clustering Performance Summary ===")
    print(summary_df.to_string(index=False))
    return results, summary_df


# ============================================================
# 7. GENERATE MULTI-CONFIGURATION PCA EMBEDDINGS
# ============================================================
X_input = X_scaled
n_samples, n_features = X_input.shape
max_components = min(n_samples, n_features)

PCA_N_COMPONENTS_LIST = [2, 3, 5, 10, 15, 20, 50]
PCA_N_COMPONENTS_LIST = [c for c in PCA_N_COMPONENTS_LIST if c <= max_components]
PCA_WHITEN_LIST       = [False, True]
total_emb = len(PCA_N_COMPONENTS_LIST) * len(PCA_WHITEN_LIST)

print("\n🔸 Computing PCA latent space projections:")
print(f"    Input features  : {n_features} clinical/demographic metrics")
print(f"    Samples         : {n_samples}")
print(f"    Max components  : {max_components}")
print(f"    Total setups    : {len(PCA_N_COMPONENTS_LIST)} components × {len(PCA_WHITEN_LIST)} whiten = {total_emb} embeddings")

embedding_matrices = {}

for n_comp in PCA_N_COMPONENTS_LIST:
    for whiten in PCA_WHITEN_LIST:
        emb_name = f"PCA_C{n_comp}_W{int(whiten)}"
        try:
            pca_model = PCA(n_components=n_comp, whiten=whiten, random_state=RANDOM_SEED)
            embedding_matrices[emb_name] = pca_model.fit_transform(X_input)
            var_exp = pca_model.explained_variance_ratio_.cumsum()[-1]
            print(f"    ✅ {emb_name} → Cumulative Explained Variance: {var_exp:.3f}")
        except Exception as e:
            print(f"    ❌ Error computing {emb_name}: {e}")
            continue

print(f"\n✅ {len(embedding_matrices)} PCA embeddings successfully generated.")

# ============================================================
# 8. EXECUTION: CLUSTERING ACROSS PCA EMBEDDINGS
# ============================================================

df_ids = (
    df_aug[['sample']].copy().rename(columns={'sample': 'ModelName'})
    if 'sample' in df_aug.columns
    else df_aug[[id_cols_to_keep[0]]].copy().rename(columns={id_cols_to_keep[0]: 'ModelName'})
)

df_clusters_master = df_ids.copy()
all_clustering_results = []

for emb_name_base, X_emb in embedding_matrices.items():
    print(f"\n🔹 Evaluating Clustering on Embedding: {emb_name_base}")
    results, summary_df = clustering_evaluation_table(
        X_emb, alg_classes, param_grids, seeds_list=SEEDS_TO_TEST
    )

    for index, row in summary_df.iterrows():
        if row['Algorithm'] not in results:
            continue

        seed_str      = f"_S{row['Best Seed']}" if row['Best Seed'] is not None else ""
        noise_pct_val = row['Noise %'] if pd.notna(row['Noise %']) else 0.0
        db_val        = row['Davies-Bouldin Score']    if pd.notna(row['Davies-Bouldin Score'])    else 0.0
        ch_val        = row['Calinski-Harabasz Score'] if pd.notna(row['Calinski-Harabasz Score']) else 0.0
        score_str     = f"S{row['Silhouette Score']:.2f}_DB{db_val:.2f}_CH{ch_val:.0f}"
        full_name     = f"{row['Algorithm']}_{emb_name_base}_K{row['Clusters']}_{score_str}{seed_str}"

        all_clustering_results.append({
            'Configuration': full_name, 'Algorithm': row['Algorithm'],
            'Embedding': emb_name_base, 'Silhouette Score': row['Silhouette Score'],
            'Davies-Bouldin Score': db_val, 'Calinski-Harabasz Score': ch_val,
            'Clusters': row['Clusters'], 'Noise_Pct': noise_pct_val,
            'Labels': results[row['Algorithm']]['best_labels'],
            'Embedding_Key': emb_name_base, 'Best Seed': row['Best Seed']
        })

        if row['Silhouette Score'] > 0:
            df_clusters_master[f"Cluster_{full_name}"] = results[row['Algorithm']]['best_labels']

output_filename_master = os.path.join(OUTPUT_DIR, "pacientes_clusterizados_todos_sinfiltro.csv")
df_clusters_master.to_csv(output_filename_master, index=False)
print(f"\n💾 Master cluster labels exported to: {output_filename_master}")

# ============================================================
# 9. FILTERING & EXPORT OF HIGH-CONFIDENCE MODELS
# ============================================================

df_all_results = pd.DataFrame(all_clustering_results)
SILHOUETTE_THRESHOLD = 0.1

df_selected_models = (
    df_all_results[df_all_results['Silhouette Score'] >= SILHOUETTE_THRESHOLD]
    .sort_values(
        by=['Silhouette Score', 'Calinski-Harabasz Score', 'Davies-Bouldin Score'],
        ascending=[False, False, True]
    )
    .copy()
    .reset_index(drop=True)
)

print(f"\n✨ Identified {len(df_selected_models)} clustering partitions with Silhouette >= {SILHOUETTE_THRESHOLD:.1f}.")

if not df_selected_models.empty:
    df_selected_labels = df_ids.copy()
    for index, row in df_selected_models.iterrows():
        config_name, labels = row['Configuration'], row['Labels']
        if len(labels) == df_ids.shape[0]:
            df_selected_labels[config_name] = labels
        else:
            print(f"⚠️  Label count for {config_name} ({len(labels)}) does not match sample count ({df_ids.shape[0]}). Skipping.")

    output_filename_selected = os.path.join(OUTPUT_DIR, "pacientes_clusterizados_seleccion_sinfiltro2.csv")
    df_selected_labels.to_csv(output_filename_selected, index=False)
    print(f"\n💾 Filtered cluster assignments ({len(df_selected_models)} models) saved to: {output_filename_selected}")
else:
    print(f"\n⚠️  No clustering solutions satisfied the Silhouette threshold >= {SILHOUETTE_THRESHOLD:.1f}.")

# ============================================================
# 10. VISUALIZATION (TOP 5 CLUSTERING CONFIGURATIONS)
# ============================================================

df_top_to_plot = df_selected_models.head(5).copy()

def generate_pca_plots(df_top_models: pd.DataFrame, embedding_matrices: dict, output_dir: str):
    """
    Renders projection plots for the top clustering configurations:
      - 2 components: 2D Scatter Plot
      - 3 components: 3D Projection
      - >3 components: Pairwise scatter grid of the first 3 Principal Components
    """
    if df_top_models.empty:
        return

    print(f"\n📈 Rendering PCA cluster visualizations for top {len(df_top_models)} models...")

    for i, row in df_top_models.reset_index(drop=True).iterrows():
        emb_key     = row['Embedding']
        labels      = row['Labels']
        config_name = row['Configuration']

        if emb_key not in embedding_matrices:
            print(f"    ❌ Embedding matrix '{emb_key}' missing from dictionary. Skipping.")
            continue

        X_pca  = embedding_matrices[emb_key]
        n_dims = X_pca.shape[1]

        title = (
            f"TOP {i+1}: {config_name}\n"
            f"Silhouette: {row['Silhouette Score']:.3f} | K={int(row['Clusters'])} | "
            f"Noise: {row['Noise_Pct']:.1f}%"
        )

        safe_name = re.sub(r'[\\/*?:"<>|]', '_', config_name)
        filename  = os.path.join(output_dir, f"TOP_{i+1}_{safe_name}.png")

        df_plot = pd.DataFrame(
            X_pca[:, :min(n_dims, 3)],
            columns=[f'PC{j+1}' for j in range(min(n_dims, 3))]
        )
        df_plot['Cluster'] = labels.astype(str)

        cluster_list = sorted(df_plot['Cluster'].unique())
        non_noise    = [c for c in cluster_list if c != '-1']
        palette      = sns.color_palette("tab10", n_colors=max(len(non_noise), 1))
        color_map    = {c: palette[j] for j, c in enumerate(non_noise)}
        if '-1' in cluster_list:
            color_map['-1'] = 'gray'

        # 2D Projection
        if n_dims == 2:
            plt.figure(figsize=(10, 8))
            sns.scatterplot(
                x='PC1', y='PC2', hue='Cluster', data=df_plot,
                palette=color_map, legend="full", alpha=0.7, s=50
            )
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(title, fontsize=12)
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"    > 2D visualization saved: {filename}")

        # 3D Projection
        elif n_dims == 3:
            fig = plt.figure(figsize=(12, 10))
            ax  = fig.add_subplot(111, projection='3d')
            for cl in cluster_list:
                sub = df_plot[df_plot['Cluster'] == cl]
                ax.scatter(
                    sub['PC1'], sub['PC2'], sub['PC3'],
                    label=f'Cluster {cl}',
                    color=color_map.get(cl, 'black'), alpha=0.7, s=50
                )
            ax.set_title(title, fontsize=12)
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC2')
            ax.set_zlabel('PC3')
            ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"    > 3D visualization saved: {filename}")

        # Pairwise PC1-PC3 grid for higher-dimensional spaces
        else:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            pairs = [('PC1', 'PC2'), ('PC1', 'PC3'), ('PC2', 'PC3')]
            for ax, (px, py) in zip(axes, pairs):
                for cl in cluster_list:
                    sub = df_plot[df_plot['Cluster'] == cl]
                    ax.scatter(
                        sub[px], sub[py],
                        label=f'Cluster {cl}',
                        color=color_map.get(cl, 'black'), alpha=0.7, s=30
                    )
                ax.set_xlabel(px)
                ax.set_ylabel(py)
                ax.set_title(f'{px} vs {py}')
            axes[-1].legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
            fig.suptitle(title, fontsize=11)
            plt.tight_layout()
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"    > Pairwise PC1-3 projection saved: {filename}")


generate_pca_plots(df_top_to_plot, embedding_matrices, OUTPUT_DIR)

# ============================================================
# 11. UNIFIED CLINICAL & CLUSTERING MASTER DATASET EXPORT
# ============================================================

print("\n🔗 Compiling harmonized master dataset...")

sample_col = 'sample' if 'sample' in df_aug.columns else id_cols_to_keep[0]

# Essential clinical & metadata columns to preserve in the master record
PRIORITY_COLS = [
    # Molecular subtyping metadata
    'Menopausal Status', 'Cancer Type', 'ER', 'PR', 'HER2', 'Subtype', 'Genetic Ancestry',
    # Survival metrics
    'OS.time', 'OS',
    # Demographics
    'gender.demographic', 'sex.demographic', 'race.demographic', 'ethnicity.demographic',
    'age_at_diagnosis.diagnoses', 'days_to_birth.demographic',
    # Staging
    'ajcc_pathologic_stage.diagnoses', 'ajcc_pathologic_t.diagnoses',
    'ajcc_pathologic_n.diagnoses',     'ajcc_pathologic_m.diagnoses',
    # Receptors and molecular profiles
    'ER_PR_HER2_Combo', 'ER_PR_HER2_Combo_encoded', 'Subtype_encoded',
    # Engineered binary indicators
    'Prior_Treatment_Flag', 'Metastasis_Flag',
]

df_clinica_final = df_aug.drop_duplicates(subset=[sample_col]).copy()

# Re-attach survival metrics if excluded during feature preprocessing
survival_cols = [
    c for c in ['OS.time', 'OS']
    if c in df_final.columns and c not in df_clinica_final.columns
]
if survival_cols:
    df_clinica_final = pd.merge(
        df_clinica_final,
        df_final[[sample_col] + survival_cols].drop_duplicates(subset=[sample_col]),
        on=sample_col, how='left'
    )

# Re-attach molecular metadata columns if dropped or missing in df_aug
meta_cols_missing_in_aug = [
    c for c in available_meta_cols
    if c not in df_clinica_final.columns or df_clinica_final[c].isna().all()
]
if meta_cols_missing_in_aug and not df_final.empty:
    df_clinica_final = df_clinica_final.drop(columns=meta_cols_missing_in_aug, errors='ignore')
    df_clinica_final = pd.merge(
        df_clinica_final,
        df_final[[sample_col] + meta_cols_missing_in_aug].drop_duplicates(subset=[sample_col]),
        on=sample_col, how='left'
    )

# Master join: cluster labels with harmonized clinical metadata
df_final_unificado = pd.merge(
    df_clusters_master,
    df_clinica_final,
    left_on='ModelName', right_on=sample_col, how='inner'
)

if sample_col in df_final_unificado.columns and sample_col != 'ModelName':
    df_final_unificado = df_final_unificado.drop(columns=[sample_col])

# Audit priority column retention and completeness
print("\n📋 Priority column completeness audit:")
found_cols   = []
missing_cols = []

for col in PRIORITY_COLS:
    candidates = [c for c in df_final_unificado.columns if col.lower() in c.lower()]
    if col in df_final_unificado.columns:
        n_na  = df_final_unificado[col].isna().sum()
        n_tot = len(df_final_unificado)
        pct   = n_na / n_tot * 100
        status = f"✅ {col:<35} → {n_tot - n_na}/{n_tot} populated ({100 - pct:.1f}% complete)"
        found_cols.append(col)
        print(status)
    elif candidates:
        print(f"⚠️  {col:<35} → Approximate match: {candidates}")
        missing_cols.append(col)
    else:
        print(f"❌ {col:<35} → NOT found")
        missing_cols.append(col)

print(f"\n✅ Priority columns retained : {len(found_cols)}")
print(f"❌ Priority columns absent   : {len(missing_cols)}")
if missing_cols:
    print(f"    Missing list: {missing_cols}")

# Export master integrated file
output_master_unificado = os.path.join(OUTPUT_DIR, "DATA_MASTER_CLUSTERS_Y_CLINICA.csv")
df_final_unificado.to_csv(output_master_unificado, index=False)

print(f"\n✅ Master dataset exported: {output_master_unificado}")
print(f"📊 Matrix dimensions: {df_final_unificado.shape[0]} samples × {df_final_unificado.shape[1]} columns")
print(f"\nFirst 10 columns : {df_final_unificado.columns[:10].tolist()} ...")
print(f"Last 10 columns  : {df_final_unificado.columns[-10:].tolist()} ...")
print("\n🎉 Clinical clustering pipeline executed successfully!")