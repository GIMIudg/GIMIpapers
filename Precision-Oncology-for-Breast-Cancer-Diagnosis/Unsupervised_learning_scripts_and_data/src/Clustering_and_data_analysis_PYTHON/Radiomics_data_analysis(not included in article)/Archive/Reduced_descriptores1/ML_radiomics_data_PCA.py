# ============================================================
# 🚀 Complete pipeline: Clinical/Demographic + Radiomics Clustering Analysis (PCA)
# INTEGRATED VERSION WITH RADIOMICS — Merge by patient_id (12 chars)
# ============================================================

# =========================
# 1. LIBRERÍAS Y CONFIGURACIÓN GLOBAL
# =========================
from sklearn.model_selection import ParameterGrid
import pandas as pd
import numpy as np
import os
import re
from sklearn.preprocessing import StandardScaler, OneHotEncoder, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering, Birch, DBSCAN, MeanShift, AffinityPropagation
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.impute import SimpleImputer
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

try:
    import hdbscan
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False
    print("⚠️ Module 'hdbscan' not found. This algorithm will be skipped.")

# PATH AND SEED CONFIGURATION
from pathlib import Path as _Path
_SCRIPT_DIR = _Path(__file__).resolve().parent
_DATA_ROOT  = _SCRIPT_DIR.parents[2] / 'Clinical_data_and_models_ids'
# Clinical data files (TCGA clinical, survival, metadata, model IDs)
PATH_BASE   = str(_DATA_ROOT / 'Clinical_Data') + '/'
PATH_FEATURES = str(_DATA_ROOT / 'Metabolic_Data' / 'FeatureMatrix_TumorPhenotype_All.csv')
PATH_CLINICAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.clinical.tsv')
PATH_SURVIVAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.survival.tsv.gz')
PATH_RADIOMICS = str(_DATA_ROOT / 'Radiomics_Data' / "TCGA-Run-2014_91cases_features_UChicago-V2010-MRI-Workstation.xls")
MODELS_NAME = str(_DATA_ROOT / 'GEMs_Data_for_construction' / "Model's_ids.txt")
# Note: the radiomics Excel file should also be placed in Clinical_data_and_models_ids/Clinical_Data/
COL_SAMPLE_ID = 'sample'
COL_SAMPLE_TYPE = 'sample_type.samples'
OUTPUT_DIR = "Results_clustering_PCA_radiomics"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SEEDS_TO_TEST = [42, 123, 100]
RANDOM_SEED = SEEDS_TO_TEST[0]
np.random.seed(RANDOM_SEED)

# =======================================================
# 2. DATA LOADING AND MASTER DATABASE CREATION
# =======================================================

def load_data(full_path):
    try:
        if full_path.endswith(('.xlsx', '.xls')):
            return pd.read_excel(full_path)
        elif full_path.endswith('.gz'):
            return pd.read_csv(full_path, sep="\t", compression='gzip')
        elif full_path.endswith('.txt'):
            # ajusta el separador según cómo esté guardado Model's_ids.txt
            return pd.read_csv(full_path, sep="\t")
        else:
            return pd.read_csv(full_path, sep="\t")
    except FileNotFoundError:
        print(f"❌ ERROR: File not found: {full_path}")
        return pd.DataFrame()
    except Exception as e:
        print(f"❌ ERROR loading '{full_path}': {e}")
        return pd.DataFrame()

df_clinical     = load_data(PATH_CLINICAL)
df_survival     = load_data(PATH_SURVIVAL)
df_metadata_raw = load_data(str(_DATA_ROOT / 'Clinical_Data' / 'MetaData.xlsx'))
df_model_names  = load_data(MODELS_NAME)
df_radiomics    = load_data(PATH_RADIOMICS)

for name, df in [("Clinical", df_clinical), ("Survival", df_survival), ("Model names", df_model_names)]:
    if df.empty:
        raise ValueError(f"❌ The file '{name}' is essential and could not be loaded. Check the path.")

print(f"Datasets loaded: Clinical ({len(df_clinical)}), Survival ({len(df_survival)}), "
      f"Metadata ({len(df_metadata_raw)}), Radiomics ({len(df_radiomics)})")


def extract_sample_id(filename):
    match = re.search(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}[A-Z0-9]?)', filename)
    if match:
        return match.group(0)[:16]
    return filename.split('_')[0].strip()[:16]


lista_modelos_unicos = df_model_names.iloc[:, 0].dropna().astype(str).tolist()
model_sample_ids = [extract_sample_id(model) for model in lista_modelos_unicos]
df_modelos_base = pd.DataFrame({
    COL_SAMPLE_ID: model_sample_ids,
    'model_full_path': lista_modelos_unicos,
})

# =======================================================
# 3. MAIN MERGE (clinical + survival + metadata)
# =======================================================

metadata_cols_to_keep = [
    'Menopausal Status', 'Cancer Type', 'ER', 'PR', 'HER2', 'Subtype',
    'Genetic Ancestry', 'Survival Status', 'Survival Time (years)', 'Sex'
]

if not df_metadata_raw.empty and 'hidden' in df_metadata_raw.columns:
    df_metadata_clean = df_metadata_raw.copy()
    df_metadata_clean['temp_id'] = (df_metadata_clean['hidden'].astype(str)
                                    .str.replace(r'\.', '-', regex=True).str.slice(0, 16))
    df_metadata_clean.rename(columns={'temp_id': COL_SAMPLE_ID}, inplace=True)
    df_metadata_clean = df_metadata_clean.drop_duplicates(subset=[COL_SAMPLE_ID], keep='first')
    available_meta_cols = [c for c in metadata_cols_to_keep if c in df_metadata_clean.columns]
    df_metadata_clean = df_metadata_clean[[COL_SAMPLE_ID] + available_meta_cols].copy()
    df_metadata_clean['patient_id'] = df_metadata_clean[COL_SAMPLE_ID].str.slice(0, 12)
    df_metadata_by_patient = (df_metadata_clean
                               .drop_duplicates(subset=['patient_id'], keep='first')
                               [['patient_id'] + available_meta_cols].copy())
    print(f"✅ Metadata loaded: {len(df_metadata_clean)} samples | "
          f"{df_metadata_clean['patient_id'].nunique()} unique patients")
else:
    print("❌ Metadata not loaded. Continuing without molecular metadata.")
    df_metadata_clean      = pd.DataFrame()
    df_metadata_by_patient = pd.DataFrame()
    available_meta_cols    = []

if 'sample' in df_clinical.columns:
    df_clinical['sample'] = df_clinical['sample'].apply(lambda x: extract_sample_id(str(x)))
if 'sample' in df_survival.columns:
    df_survival['sample'] = df_survival['sample'].apply(lambda x: extract_sample_id(str(x)))

df_survival_clean = (df_survival[['sample', 'OS.time', 'OS']]
                     .drop_duplicates(subset=['sample'], keep='first')
                     .rename(columns={'sample': COL_SAMPLE_ID}))

df_merged_clinical = pd.merge(
    df_modelos_base,
    df_clinical.drop(columns=['id', 'case_id'], errors='ignore'),
    on=COL_SAMPLE_ID, how='left'
)
df_final = pd.merge(df_merged_clinical, df_survival_clean, on=COL_SAMPLE_ID, how='left')

# Merge 1: por Sample ID exacto (16 chars)
if not df_metadata_clean.empty:
    cols_overlap = [col for col in available_meta_cols if col in df_final.columns]
    df_final = pd.merge(
        df_final.drop(columns=cols_overlap, errors='ignore'),
        df_metadata_clean.drop(columns=['patient_id']),
        on=COL_SAMPLE_ID, how='left'
    )

# Merge 2: relleno por Patient ID (12 chars)
if not df_metadata_by_patient.empty:
    df_final['patient_id'] = df_final[COL_SAMPLE_ID].str.slice(0, 12)
    missing_mask = df_final[available_meta_cols].isna().any(axis=1)
    n_missing_before = missing_mask.sum()
    df_fill = pd.merge(
        df_final.loc[missing_mask, ['patient_id']],
        df_metadata_by_patient, on='patient_id', how='left'
    )
    df_fill.index = df_final.index[missing_mask]
    for col in available_meta_cols:
        if col in df_fill.columns:
            df_final.loc[missing_mask, col] = df_final.loc[missing_mask, col].fillna(df_fill[col])
    df_final = df_final.drop(columns=['patient_id'])
    n_missing_after = df_final[available_meta_cols].isna().any(axis=1).sum()
    print(f"\n🔁 Patient ID fill-in:")
    print(f"   Rows with incomplete metadata BEFORE : {n_missing_before}")
    print(f"   Rows with incomplete metadata AFTER  : {n_missing_after}")
    print(f"   Rows recovered                       : {n_missing_before - n_missing_after}")
    print("\n📊 Remaining NAs per metadata column:")
    for col in available_meta_cols:
        if col in df_final.columns:
            n_na = df_final[col].isna().sum()
            print(f"   {col:<25} → {n_na} NAs ({n_na/len(df_final)*100:.1f}%)")

df = df_final.copy()
df_filtered = df.copy()

# =========================================================
# 3.1 Conversión de Tipos
# =========================================================
for c in ['is_ffpe.samples', 'oct_embedded.samples']:
    if c in df_filtered.columns:
        df_filtered[c] = df_filtered[c].replace({True: 1, False: 0})

for c in ['age_at_diagnosis.diagnoses', 'days_to_birth.demographic']:
    if c in df_filtered.columns and df_filtered[c].notna().any():
        if (df_filtered[c].dropna() > 1000).any():
            df_filtered[c] = df_filtered[c] / 365.25

# =======================================================
# 3.2 RADIOMIC DATA LOADING AND NORMALIZATION
#     Strategy: merge by patient_id (12 chars)
#     The .les files do not contain the TCGA sample suffix (01A/01B),
#     so the only reliable match is via the 12-char patient_id.
# =======================================================
# adjust as necessary
DESCRIPTORES_RADIOMICS = [
    'Maximum enhancement (K1)',
    'Time to peak (K2)',
    'Washout rate (K4)',
    'Maximum enhancement-variance (E1)',
    'Enhancement-Variance Time to Peak (E2)',
    'Entropy (T6)',
    'Variance (T14)',
    'Sphericity (G1)',
    'Margin Sharpness (M1)',
    'Variance of Margin Sharpness (M2)',
    'Size/Lesion volume (S1)',
    'Volume of most enhancing voxels (S4)',

]

""" lista completa
DESCRIPTORES_RADIOMICS = [
    'Maximum enhancement (K1)',
    'Time to peak (K2)',
    'Uptake rate (K3)',
    'Washout rate (K4)',
    'Curve shape index (K5)',
    'E1 (K6)',
    'Signal Enhancement Ratio (SER) (K7)',
    'Maximum enhancement-variance (E1)',
    'Enhancement-Variance Time to Peak (E2)',
    'Enhancement-variance Increasing Rate (E3)',
    'Enhancement-variance Decreasing Rate (E4)',
    'Contrast (T1)',
    'Correlation (T2)',
    'Difference Entropy (T3)',
    'Difference Variance (T4)',
    'Energy (T5)',
    'Entropy (T6)',
    'Homogeneity (T7)',
    'IMC1 (T8)',
    'IMC2 (T9)',
    'Maximum Correlation Coefficient (T10)',
    'Sum Average (T11)',
    'Sum Entropy (T12)',
    'Sum Variance (T13)',
    'Variance (T14)',
    'Sphericity (G1)',
    'Irregularity (G2)',
    'Margin Sharpness (M1)',
    'Variance of Margin Sharpness (M2)',
    'Variance of Radial Gradient Histogram (vRGH) (M3)',
    'Size/Lesion volume (S1)',
    'Effective Diameter (S2)',
    'Surface Area (S3)',
    'Surface Area to Volume ratio (G3)',
    'Volume of most enhancing voxels (S4)',
    'Maximum Diameter (S5)',
]
"""

def extract_patient_id_from_les(les_name):
    """
    Extrae los 12 chars del patient_id desde un nombre de archivo .les.
    Ejemplos:
      TCGA-AO-A03M-1.les     → TCGA-AO-A03M
      TCGA-AO-A0J9-S2-1.les  → TCGA-AO-A0J9  (ignora el sufijo S2 y el índice)
    """
    name = str(les_name).strip()
    match = re.search(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', name)
    if match:
        return match.group(1)   # siempre 12 chars exactos
    return name[:12]


def load_and_merge_radiomics(df_base, df_radiomics_raw, descriptores, col_sample_id='sample'):
    """
    Normaliza los IDs del archivo de radiomics a patient_id (12 chars)
    y hace el merge con df_base también por patient_id.
    Devuelve (df_merged, available_descriptors).
    """
    if df_radiomics_raw.empty:
        print("⚠️  df_radiomics está vacío. Se omite el merge de radiomics.")
        return df_base.copy(), []

    df_rad = df_radiomics_raw.copy()

    # ── Detectar columna de nombre de lesión ──────────────────────────────
    lesion_col_candidates = [c for c in df_rad.columns
                             if 'lesion' in c.lower() and 'name' in c.lower()]
    if lesion_col_candidates:
        lesion_col = lesion_col_candidates[0]
    else:
        lesion_col = df_rad.columns[0]
        print(f"⚠️  Columna 'Lesion Name' no encontrada. Usando primera columna: '{lesion_col}'")

    print(f"✅ Columna de ID en radiomics detectada: '{lesion_col}'")

    # ── Extraer patient_id (12 chars) desde el nombre .les ───────────────
    df_rad['patient_id'] = df_rad[lesion_col].astype(str).apply(extract_patient_id_from_les)
    df_rad['lesion_name_original'] = df_rad[lesion_col]

    # Verificar que todos tengan 12 chars
    bad_ids = df_rad[df_rad['patient_id'].str.len() != 12]['patient_id'].unique()
    if len(bad_ids) > 0:
        print(f"⚠️  {len(bad_ids)} IDs con longitud ≠ 12 chars: {bad_ids}")

    # ── Verificar descriptores disponibles ───────────────────────────────
    available_descriptors = [c for c in descriptores if c in df_rad.columns]
    missing_descriptors   = [c for c in descriptores if c not in df_rad.columns]

    if missing_descriptors:
        print(f"\n⚠️  Descriptores no encontrados en radiomics ({len(missing_descriptors)}):")
        for m in missing_descriptors:
            print(f"   - {m}")

    print(f"✅ Descriptores disponibles: {len(available_descriptors)}/{len(descriptores)}")

    if not available_descriptors:
        raise ValueError("❌ Ningún descriptor radiómico encontrado. Verifica los nombres de columnas.")

    # ── Mantener solo patient_id + descriptores; un registro por paciente ─
    df_rad_slim = (df_rad[['patient_id', 'lesion_name_original'] + available_descriptors]
                   .drop_duplicates(subset=['patient_id'], keep='first'))

    print(f"\nEjemplos de mapeo de IDs radiomics:")
    print(f"  {'Lesion Name':<35} → {'patient_id (12)'}")
    print(f"  {'-'*55}")
    for _, row in df_rad_slim.head(5).iterrows():
        print(f"  {str(row['lesion_name_original']):<35} → {row['patient_id']}")

    # ── Agregar patient_id al df_base para el merge ───────────────────────
    df_base_work = df_base.copy()
    df_base_work['patient_id'] = df_base_work[col_sample_id].str.slice(0, 12)

    # ── Merge por patient_id ──────────────────────────────────────────────
    df_merged = pd.merge(
        df_base_work,
        df_rad_slim,
        on='patient_id',
        how='inner'
    )

    # Limpiar columna auxiliar
    df_merged = df_merged.drop(columns=['patient_id'])

    n_original = len(df_base)
    n_merged   = len(df_merged)
    print(f"\n🔗 Merge radiomics completado:")
    print(f"   Muestras en df_filtered      : {n_original}")
    print(f"   Muestras con match radiomics : {n_merged}")
    print(f"   Muestras sin match           : {n_original - n_merged}")

    if n_merged == 0:
        print("\n🔍 DIAGNÓSTICO — patient_ids en df_filtered (primeros 10):")
        print(df_base_work['patient_id'].head(10).tolist())
        print("\n🔍 DIAGNÓSTICO — patient_ids en radiomics (primeros 10):")
        print(df_rad_slim['patient_id'].head(10).tolist())
        raise ValueError("❌ Merge resultó en 0 filas. Revisa el formato de los IDs.")

    return df_merged, available_descriptors


# Ejecutar merge de radiomics
df_filtered_rad, available_descriptors = load_and_merge_radiomics(
    df_base=df_filtered,
    df_radiomics_raw=df_radiomics,
    descriptores=DESCRIPTORES_RADIOMICS,
    col_sample_id=COL_SAMPLE_ID
)

# =========================
# 4. SELECCIÓN DE DESCRIPTORES
# =========================

# Usar SOLO los descriptores radiómicos como features de clustering
# (más las columnas de ingeniería de features que se agregarán en la sección 5)
descriptores_finales = available_descriptors   # los 36 features radiómicos confirmados

# Columnas de ID a preservar
id_cols_to_keep = [c for c in ['submitter_id', 'sample', 'modelo_path_completo']
                   if c in df_filtered_rad.columns]

final_cols = [c for c in descriptores_finales if c in df_filtered_rad.columns]

if not final_cols:
    raise ValueError("No se encontraron columnas de descriptores válidas en el DataFrame.")

df_aug = df_filtered_rad[final_cols + id_cols_to_keep +
                          [c for c in df_filtered_rad.columns
                           if c not in final_cols + id_cols_to_keep]].copy()

print(f"\n✅ df_aug construido: {df_aug.shape[0]} muestras × {df_aug.shape[1]} columnas")
print(f"   Features radiómicos incluidos: {len(final_cols)}")

# =========================
# 5. INGENIERÍA DE FEATURES
# =========================

def prior_treatment_flag(r):
    cols = ['prior_malignancy.diagnoses', 'prior_treatment.diagnoses',
            'progression_or_recurrence.diagnoses']
    for c in cols:
        if c in r.index and pd.notna(r[c]):
            val = str(r[c]).lower().strip()
            if val in ['yes', 'true', 'had prior treatment', 'recurrence', 'progression']:
                return 1
    return 0

df_aug['Prior_Treatment_Flag'] = df_aug.apply(prior_treatment_flag, axis=1)

if all(c in df_aug.columns for c in ['ajcc_pathologic_m.diagnoses', 'ajcc_pathologic_n.diagnoses']):
    df_aug['Metastasis_Flag'] = df_aug.apply(
        lambda r: 1 if ('m1' in str(r['ajcc_pathologic_m.diagnoses']).lower() or
                         'n2' in str(r['ajcc_pathologic_n.diagnoses']).lower() or
                         'n3' in str(r['ajcc_pathologic_n.diagnoses']).lower()) else 0, axis=1)
else:
    df_aug['Metastasis_Flag'] = 0

receptor_cols = ['ER', 'PR', 'HER2']
if all(c in df_aug.columns for c in receptor_cols):
    def classify_molecular_subtype(row):
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

# =========================
# 6. PREPROCESAMIENTO PARA ML
# =========================

label_cols_candidates = ['Molecular_Subtype', 'ER_PR_HER2_Combo', 'Subtype']
label_encoders = {}

for col in label_cols_candidates:
    if col in df_aug.columns and df_aug[col].dtype in ['object', 'category']:
        le = LabelEncoder()
        df_aug[f"{col}_encoded"] = le.fit_transform(df_aug[col].fillna('Unknown').astype(str))
        label_encoders[col] = le

id_and_meta_cols = set(id_cols_to_keep + ['submitter_id', 'sample', 'modelo_path_completo',
                                           'lesion_name_original'])

numeric_cols = [c for c in df_aug.select_dtypes(
                    include=['int64', 'float64', 'float32', 'int32', 'uint8']).columns
                if c not in id_and_meta_cols]

label_orig_cols = [col for col in label_cols_candidates
                   if col in df_aug.columns and f"{col}_encoded" in df_aug.columns]

categorical_cols = [c for c in df_aug.select_dtypes(include=['object', 'category']).columns
                    if c not in id_and_meta_cols and c not in label_orig_cols]

# Imputación numérica
cols_to_impute_zero = [c for c in numeric_cols if df_aug[c].isna().any()]
if cols_to_impute_zero:
    imputer_num = SimpleImputer(strategy='constant', fill_value=0)
    df_aug.loc[:, cols_to_impute_zero] = imputer_num.fit_transform(
        df_aug.loc[:, cols_to_impute_zero])

# Imputación categórica
for col in categorical_cols:
    df_aug[col] = df_aug[col].fillna('Missing').astype(str)

import sklearn
sklearn_version = tuple(int(x) for x in sklearn.__version__.split(".")[:2])
if sklearn_version >= (1, 2):
    ohe_kwargs = {'handle_unknown': 'ignore', 'sparse_output': False}
else:
    ohe_kwargs = {'handle_unknown': 'ignore', 'sparse': False}

transformers = []
if numeric_cols:
    transformers.append(('num', StandardScaler(), numeric_cols))
if categorical_cols:
    transformers.append(('cat', OneHotEncoder(**ohe_kwargs), categorical_cols))

if not transformers:
    raise ValueError("No hay features numéricas o categóricas válidas para el preprocesamiento.")

preprocessor = ColumnTransformer(transformers, remainder='drop')
X_scaled = preprocessor.fit_transform(df_aug)

final_columns = []
if 'num' in preprocessor.named_transformers_:
    final_columns.extend(numeric_cols)
if 'cat' in preprocessor.named_transformers_ and preprocessor.named_transformers_['cat'] is not None:
    cat_feature_names = preprocessor.named_transformers_['cat'].get_feature_names_out(categorical_cols)
    final_columns.extend(cat_feature_names)

final_columns = np.array(final_columns)
print(f"\n✅ ColumnTransformer completado: {len(final_columns)} features listos para ML.")
print(f"   (de los cuales {len(available_descriptors)} son features radiómicos)")

# =========================
# 7. EVALUACIÓN DE CLUSTERING
# =========================

K_VALORES = range(2, 10)

param_grids = {
    "KMeans":                  {'n_clusters': K_VALORES},
    "Agglomerative":           {'n_clusters': K_VALORES, 'linkage': ['ward', 'average', 'complete']},
    "Birch":                   {'n_clusters': K_VALORES},
    "GMM":                     {'n_components': K_VALORES},
    "BayesianGaussianMixture": {'n_components': K_VALORES},
    "DBSCAN":                  {'eps': [0.5, 1.0, 1.5], 'min_samples': [5, 10, 20]},
    "MeanShift":               {'bandwidth': [None]},
    "AffinityPropagation":     {'damping': [0.5, 0.9]},
}

if HDBSCAN_AVAILABLE:
    param_grids["HDBSCAN"] = {'min_cluster_size': [5, 10, 20], 'min_samples': [None, 5, 10]}

alg_classes = {
    'KMeans':                  KMeans,
    'Agglomerative':           AgglomerativeClustering,
    'Birch':                   Birch,
    'GMM':                     GaussianMixture,
    'DBSCAN':                  DBSCAN,
    'MeanShift':               MeanShift,
    'AffinityPropagation':     AffinityPropagation,
    'BayesianGaussianMixture': BayesianGaussianMixture,
}

if HDBSCAN_AVAILABLE:
    alg_classes['HDBSCAN'] = hdbscan.HDBSCAN

alg_classes  = {k: v for k, v in alg_classes.items() if v is not None}
param_grids  = {k: v for k, v in param_grids.items() if v}


def clustering_evaluation_table(X, algorithms, param_grids, seeds_list):
    results = {}
    summary_rows = []
    MAX_NOISE_PCT = 10.0

    for alg_name, alg_class in algorithms.items():
        best_score     = -1
        best_params    = None
        best_labels    = None
        best_n_clusters = 0
        best_noise_pct  = 0
        best_seed       = None

        current_grid_dict = param_grids.get(alg_name, [{}])
        if isinstance(current_grid_dict, dict) and any(isinstance(v, list) for v in current_grid_dict.values()):
            param_grid = ParameterGrid(current_grid_dict)
        else:
            param_grid = ParameterGrid([current_grid_dict])

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

                    valid_mask  = labels != -1
                    n_clusters  = len(set(labels[valid_mask]))
                    noise_pct   = (labels == -1).sum() / len(labels) * 100 if -1 in labels else 0.0

                    if noise_pct > MAX_NOISE_PCT:
                        continue
                    if n_clusters > 1 and valid_mask.sum() > 1:
                        score = silhouette_score(X[valid_mask], labels[valid_mask])
                    else:
                        score = -1

                    if score > best_score:
                        best_score     = score
                        best_params    = param_val
                        best_labels    = labels
                        best_n_clusters = n_clusters
                        best_noise_pct  = noise_pct
                        best_seed       = current_seed
                except Exception:
                    continue

        if best_labels is not None:
            try:
                if best_n_clusters > 1 and np.sum(best_labels != -1) > 1:
                    mask     = best_labels != -1
                    db_score = davies_bouldin_score(X[mask], best_labels[mask])
                    ch_score = calinski_harabasz_score(X[mask], best_labels[mask])
                else:
                    db_score = np.nan
                    ch_score = np.nan
            except Exception:
                db_score = np.nan
                ch_score = np.nan

            results[alg_name] = {
                "best_score": best_score, "best_params": best_params,
                "best_labels": best_labels, "n_clusters": best_n_clusters,
                "noise_pct": best_noise_pct, "best_seed": best_seed,
                "db_score": db_score, "ch_score": ch_score
            }
            summary_rows.append({
                "Algorithm":               alg_name,
                "Silhouette Score":        best_score,
                "Davies-Bouldin Score":    db_score,
                "Calinski-Harabasz Score": ch_score,
                "Best Params":             best_params,
                "Clusters":                best_n_clusters,
                "Noise %":                 best_noise_pct,
                "Best Seed":               best_seed,
                "Unique Labels":           np.unique(best_labels)
            })

    summary_df = pd.DataFrame(summary_rows).sort_values(by="Silhouette Score", ascending=False)
    print("\n=== Resumen de clustering ===")
    print(summary_df.to_string(index=False))
    return results, summary_df

# =======================================================
# 8. GENERAR EMBEDDINGS PCA
# =======================================================
X_input = X_scaled

n_samples, n_features = X_input.shape
max_components = min(n_samples, n_features)

PCA_N_COMPONENTS_LIST = [2, 3, 5, 10, 15, 20, 50]
PCA_N_COMPONENTS_LIST = [c for c in PCA_N_COMPONENTS_LIST if c <= max_components]
PCA_WHITEN_LIST       = [False, True]

total_emb = len(PCA_N_COMPONENTS_LIST) * len(PCA_WHITEN_LIST)
print(f"\n🔸 Ejecutando PCA en X_scaled con grid de hiperparámetros:")
print(f"   Features de entrada : {n_features} métricas secundarias")
print(f"   Modelos             : {n_samples}")
print(f"   Componentes máx.    : {max_components}")
print(f"   Configuraciones     : {len(PCA_N_COMPONENTS_LIST)} n_components × "
      f"{len(PCA_WHITEN_LIST)} whiten = {total_emb} embeddings")

embedding_matrices = {}

for n_comp in PCA_N_COMPONENTS_LIST:
    for whiten in PCA_WHITEN_LIST:
        emb_name = f"PCA_C{n_comp}_W{int(whiten)}"
        try:
            pca_model = PCA(n_components=n_comp, whiten=whiten, random_state=RANDOM_SEED)
            embedding_matrices[emb_name] = pca_model.fit_transform(X_input)
            var_exp = pca_model.explained_variance_ratio_.cumsum()[-1]
            print(f"   ✅ {emb_name} → varianza explicada acumulada: {var_exp:.3f}")
        except Exception as e:
            print(f"   ❌ Error en PCA {emb_name}: {e}")
            continue

print(f"\n✅ {len(embedding_matrices)} embeddings PCA generados.")

# ====================================
# 9. PIPELINE FINAL: CLUSTERING Y EXPORTACIÓN
# ====================================

df_ids = (df_aug[['sample']].copy().rename(columns={'sample': 'ModelName'})
          if 'sample' in df_aug.columns
          else df_aug[[id_cols_to_keep[0]]].copy().rename(columns={id_cols_to_keep[0]: 'ModelName'}))

df_clusters_master  = df_ids.copy()
all_clustering_results = []

for emb_name_base, X_emb in embedding_matrices.items():
    print(f"\n🔹 Clustering en Embedding: {emb_name_base}")
    results, summary_df = clustering_evaluation_table(
        X_emb, alg_classes, param_grids, seeds_list=SEEDS_TO_TEST)

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
            'Configuration':           full_name,
            'Algorithm':               row['Algorithm'],
            'Embedding':               emb_name_base,
            'Silhouette Score':        row['Silhouette Score'],
            'Davies-Bouldin Score':    db_val,
            'Calinski-Harabasz Score': ch_val,
            'Clusters':                row['Clusters'],
            'Noise_Pct':               noise_pct_val,
            'Labels':                  results[row['Algorithm']]['best_labels'],
            'Embedding_Key':           emb_name_base,
            'Best Seed':               row['Best Seed']
        })

        if row['Silhouette Score'] > 0:
            df_clusters_master[f"Cluster_{full_name}"] = results[row['Algorithm']]['best_labels']

output_filename_master = os.path.join(OUTPUT_DIR, "pacientes_clusterizados_todos_sinfiltro.csv")
df_clusters_master.to_csv(output_filename_master, index=False)
print(f"\n💾 CSV de clusters master guardado en: {output_filename_master}")

# ------------------------------------------------------------------------
# 10. FILTRADO Y EXPORTACIÓN DE LOS ALGORITMOS SELECCIONADOS
# ------------------------------------------------------------------------

df_all_results      = pd.DataFrame(all_clustering_results)
SILHOUETTE_THRESHOLD = 0.1

df_selected_models = df_all_results[
    df_all_results['Silhouette Score'] >= SILHOUETTE_THRESHOLD
].sort_values(
    by=['Silhouette Score', 'Calinski-Harabasz Score', 'Davies-Bouldin Score'],
    ascending=[False, False, True]
).copy().reset_index(drop=True)

num_selected_models = len(df_selected_models)
print(f"\n✨ Encontrados {num_selected_models} modelos con Silhouette Score >= {SILHOUETTE_THRESHOLD:.1f}.")

if not df_selected_models.empty:
    df_selected_labels = df_ids.copy()
    for index, row in df_selected_models.iterrows():
        config_name = row['Configuration']
        labels      = row['Labels']
        if len(labels) == df_ids.shape[0]:
            df_selected_labels[config_name] = labels
        else:
            print(f"⚠️ Etiquetas de {config_name} ({len(labels)}) no coinciden "
                  f"con muestras ({df_ids.shape[0]}). Saltando.")

    output_filename_selected = os.path.join(OUTPUT_DIR, "pacientes_clusterizados_seleccion_sinfiltro2.csv")
    df_selected_labels.to_csv(output_filename_selected, index=False)
    print(f"\n💾 CSV de {num_selected_models} modelos seleccionados guardado en: {output_filename_selected}")
else:
    print(f"\n⚠️ No se encontraron modelos con Silhouette Score >= {SILHOUETTE_THRESHOLD:.1f}.")

# ------------------------------------------------------------------------
# 11. GRAFICACIÓN DE LOS TOP N ALGORITMOS
# ------------------------------------------------------------------------

df_top_to_plot = df_selected_models.head(5).copy()


def generate_pca_plots(df_top_models, embedding_matrices, output_dir):
    if df_top_models.empty:
        return

    print(f"\n📈 Generando gráficas PCA para los {len(df_top_models)} modelos con mejor Score...")

    for i, row in df_top_models.reset_index(drop=True).iterrows():
        emb_key     = row['Embedding']
        labels      = row['Labels']
        config_name = row['Configuration']

        if emb_key not in embedding_matrices:
            print(f"Error: Matriz {emb_key} no encontrada.")
            continue

        X_pca  = embedding_matrices[emb_key]
        n_dims = X_pca.shape[1]

        title = (f"TOP {i+1}: {config_name}\n"
                 f"Silhouette: {row['Silhouette Score']:.3f} | K={int(row['Clusters'])} | "
                 f"Ruido: {row['Noise_Pct']:.1f}%")

        safe_config_name = re.sub(r'[\\/*?:"<>|]', '_', config_name)
        filename = os.path.join(output_dir, f"TOP_{i+1}_{safe_config_name}.png")

        df_plot = pd.DataFrame(X_pca, columns=[f'PC{j+1}' for j in range(n_dims)])
        df_plot['Cluster'] = labels.astype(str)

        cluster_list               = sorted(df_plot['Cluster'].unique())
        cluster_labels_for_palette = [c for c in cluster_list if c != '-1']
        palette   = sns.color_palette("tab10", n_colors=max(len(cluster_labels_for_palette), 1))
        color_map = {c: palette[j] for j, c in enumerate(cluster_labels_for_palette)}
        if '-1' in cluster_list:
            color_map['-1'] = 'gray'

        if n_dims == 2:
            plt.figure(figsize=(10, 8))
            sns.scatterplot(x='PC1', y='PC2', hue='Cluster', data=df_plot,
                            palette=color_map, legend="full", alpha=0.7, s=50)
            plt.xlabel('PC1'); plt.ylabel('PC2')
            plt.title(title, fontsize=12)
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"   > Gráfica 2D guardada: {filename}")

        elif n_dims == 3:
            fig = plt.figure(figsize=(12, 10))
            ax  = fig.add_subplot(111, projection='3d')
            for cluster_label in cluster_list:
                df_subset = df_plot[df_plot['Cluster'] == cluster_label]
                ax.scatter(df_subset['PC1'], df_subset['PC2'], df_subset['PC3'],
                           label=f'Cluster {cluster_label}',
                           color=color_map.get(cluster_label, 'black'), alpha=0.7, s=50)
            ax.set_title(title, fontsize=12)
            ax.set_xlabel('PC1'); ax.set_ylabel('PC2'); ax.set_zlabel('PC3')
            ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"   > Gráfica 3D guardada: {filename}")
        else:
            # Para n_dims > 3: scatter matrix de las primeras 4 PCs
            n_plot = min(n_dims, 4)
            cols_to_plot = [f'PC{j+1}' for j in range(n_plot)]
            fig, axes = plt.subplots(n_plot, n_plot, figsize=(4 * n_plot, 4 * n_plot))
            fig.suptitle(title, fontsize=11, y=1.01)
            for r in range(n_plot):
                for c in range(n_plot):
                    ax = axes[r][c]
                    if r == c:
                        ax.hist(df_plot[cols_to_plot[r]], bins=20, color='steelblue', alpha=0.7)
                        ax.set_xlabel(cols_to_plot[r])
                    else:
                        for cl in cluster_list:
                            sub = df_plot[df_plot['Cluster'] == cl]
                            ax.scatter(sub[cols_to_plot[c]], sub[cols_to_plot[r]],
                                       color=color_map.get(cl, 'black'),
                                       alpha=0.6, s=30, label=f'Cluster {cl}')
                        ax.set_xlabel(cols_to_plot[c])
                        ax.set_ylabel(cols_to_plot[r])
            handles = [plt.Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor=color_map.get(cl, 'black'), markersize=8,
                                   label=f'Cluster {cl}') for cl in cluster_list]
            fig.legend(handles=handles, title='Cluster',
                       bbox_to_anchor=(1.01, 0.5), loc='center left')
            plt.tight_layout()
            plt.savefig(filename, bbox_inches='tight')
            plt.close()
            print(f"   > Scatter matrix {n_plot}×{n_plot} PCs guardada: {filename}")


generate_pca_plots(df_top_to_plot, embedding_matrices, OUTPUT_DIR)

# ========================================================================
# 12. GENERACIÓN DE BASE DE DATOS INTEGRADA
# ========================================================================

print("\n🔗 Generando archivo maestro integrado...")

sample_col = 'sample' if 'sample' in df_aug.columns else id_cols_to_keep[0]

PRIORITY_COLS = [
    'Menopausal Status', 'Cancer Type', 'ER', 'PR', 'HER2', 'Subtype', 'Genetic Ancestry',
    'OS.time', 'OS',
    'gender.demographic', 'sex.demographic', 'race.demographic', 'ethnicity.demographic',
    'age_at_diagnosis.diagnoses', 'days_to_birth.demographic',
    'ajcc_pathologic_stage.diagnoses', 'ajcc_pathologic_t.diagnoses',
    'ajcc_pathologic_n.diagnoses',     'ajcc_pathologic_m.diagnoses',
    'ER_PR_HER2_Combo', 'ER_PR_HER2_Combo_encoded', 'Subtype_encoded',
    'Prior_Treatment_Flag', 'Metastasis_Flag',
    # Radiomics
    'Maximum enhancement (K1)', 'Signal Enhancement Ratio (SER) (K7)',
    'Sphericity (G1)', 'Irregularity (G2)', 'Size/Lesion volume (S1)',
    'Maximum Diameter (S5)', 'Margin Sharpness (M1)',
]

df_clinica_final = df_aug.drop_duplicates(subset=[sample_col]).copy()

survival_cols = [c for c in ['OS.time', 'OS'] if c in df_final.columns
                 and c not in df_clinica_final.columns]
if survival_cols:
    df_clinica_final = pd.merge(
        df_clinica_final,
        df_final[[sample_col] + survival_cols].drop_duplicates(subset=[sample_col]),
        on=sample_col, how='left'
    )

meta_cols_missing_in_aug = [c for c in available_meta_cols
                             if c not in df_clinica_final.columns
                             or df_clinica_final[c].isna().all()]
if meta_cols_missing_in_aug and not df_final.empty:
    df_clinica_final = df_clinica_final.drop(columns=meta_cols_missing_in_aug, errors='ignore')
    df_clinica_final = pd.merge(
        df_clinica_final,
        df_final[[sample_col] + meta_cols_missing_in_aug].drop_duplicates(subset=[sample_col]),
        on=sample_col, how='left'
    )

df_final_unificado = pd.merge(
    df_clusters_master,
    df_clinica_final,
    left_on='ModelName', right_on=sample_col, how='inner'
)

if sample_col in df_final_unificado.columns and sample_col != 'ModelName':
    df_final_unificado = df_final_unificado.drop(columns=[sample_col])

# Verificación de columnas prioritarias
print("\n📋 Verificación de columnas prioritarias en archivo final:")
found_cols   = []
missing_cols = []

for col in PRIORITY_COLS:
    candidates = [c for c in df_final_unificado.columns if col.lower() in c.lower()]
    if col in df_final_unificado.columns:
        n_na  = df_final_unificado[col].isna().sum()
        n_tot = len(df_final_unificado)
        pct   = n_na / n_tot * 100
        print(f"✅ {col:<45} → {n_tot - n_na}/{n_tot} con datos ({100 - pct:.1f}% completo)")
        found_cols.append(col)
    elif candidates:
        print(f"⚠️  {col:<45} → no exacto, posibles: {candidates}")
        missing_cols.append(col)
    else:
        print(f"❌  {col:<45} → NO encontrada")
        missing_cols.append(col)

print(f"\n✅ Columnas prioritarias presentes : {len(found_cols)}")
print(f"❌ Columnas prioritarias ausentes  : {len(missing_cols)}")
if missing_cols:
    print(f"   {missing_cols}")

output_master_unificado = os.path.join(OUTPUT_DIR, "DATA_MASTER_CLUSTERS_Y_CLINICA.csv")
df_final_unificado.to_csv(output_master_unificado, index=False)

print(f"\n✅ Archivo maestro guardado: {output_master_unificado}")
print(f"📊 Dimensiones: {df_final_unificado.shape[0]} muestras × {df_final_unificado.shape[1]} columnas")
print(f"\nPrimeras columnas: {df_final_unificado.columns[:10].tolist()} ...")
print(f"Últimas columnas:  {df_final_unificado.columns[-10:].tolist()} ...")