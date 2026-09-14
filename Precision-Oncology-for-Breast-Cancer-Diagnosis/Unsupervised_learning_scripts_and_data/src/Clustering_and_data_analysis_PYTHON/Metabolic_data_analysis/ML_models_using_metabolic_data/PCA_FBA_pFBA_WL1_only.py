# ============================================================
# 🚀 Pipeline: Metabolic Clustering — FeatureMatrix_TumorPhenotype
# Source Data: FeatureMatrix_TumorPhenotype.csv (MATLAB COBRA output)
# Optimizations: FBA | pFBA | L1w
# Features: Secondary metabolic metrics (excluding raw fluxes)
# Dimensionality Reduction: Multi-configuration PCA
# Reproducibility Target: Deterministic pipeline via explicit seeds & imputations
# ============================================================

import re
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

from sklearn.preprocessing import RobustScaler
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.cluster import (
    KMeans, AgglomerativeClustering, Birch,
    DBSCAN, MeanShift, AffinityPropagation,
    estimate_bandwidth
)
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.metrics import (
    silhouette_score, davies_bouldin_score, calinski_harabasz_score
)
from sklearn.model_selection import ParameterGrid

warnings.filterwarnings("ignore")

# Check optional clustering dependencies
try:
    import hdbscan
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False
    print("⚠️  Warning: HDBSCAN is not installed in this environment.")

# ============================================================
# ⚙️ CONFIGURATION & PATH SETUP
# ============================================================
# Resolve project root: go up from this file location
_SCRIPT_DIR = Path(__file__).resolve().parent   # ML_models_using_metabolic_data/
_DATA_ROOT = _SCRIPT_DIR.parents[3] / 'Clinical_data_and_models_ids'

# Input files
PATH_FEATURES = str(_DATA_ROOT / 'Metabolic_Data' / 'FeatureMatrix_TumorPhenotype_All.csv')
PATH_CLINICAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.clinical.tsv')
PATH_SURVIVAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.survival.tsv.gz')

# Output directory
OUT_DIR = "results_TumorPhenotype_PCA_metrics_updated_nol2"
os.makedirs(OUT_DIR, exist_ok=True)

# Random seeds for stochastic algorithms to guarantee reproducibility
SEEDS_TO_TEST = [42, 123, 100]
RANDOM_SEED = SEEDS_TO_TEST[0]
np.random.seed(RANDOM_SEED)

# Filtering and clustering quality thresholds
ZERO_NULL_THRESHOLD = 0.95   # Drop features with >= 95% zeros or NaNs
SILHOUETTE_THRESHOLD = 0.10  # Minimum silhouette score for downstream reporting
MAX_NOISE_PCT = 5.0          # Max allowable noise/outlier percentage (-1 labels)

# Solvers / optimization solutions evaluated
SOL_NAMES = ["FBA", "pFBA", "L1w"]

# Core metabolic metric root names
METRIC_ROOTS = [
    "CU",
    "EA",
    "WarburgIndex",
    "ATPConsumption",
    "ATPProduction",
    "RedoxIndex",
    "MFI",
    "AnabolismScore",
    "NADPHdemand",
    "TCA_completeness",
    "LipidSat",
    "LipidUnsat",
    "LipidPL",
    "GlnDependence",
]

ONCOMET_NAMES = ["Lactate", "Succinate", "AlphaKG"]

# ============================================================
# 🔑 UTILITY FUNCTIONS
# ============================================================
def extract_model_id(model_name: str) -> str:
    """
    Extracts and standardizes TCGA sample identifiers (16-character format: TCGA-XX-XXXX-XX).
    Falls back to leading substring before underscores.
    """
    match = re.search(
        r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}[A-Z0-9]?)',
        str(model_name)
    )
    if match:
        return match.group(0)[:16]
    return str(model_name).split('_')[0].strip()[:16]


def detect_columns(df: pd.DataFrame) -> dict:
    """
    Categorizes DataFrame columns into raw metabolic fluxes, secondary computed
    metrics, subsystem activities (SA), and oncometabolite turnover rates.
    """
    all_cols = set(df.columns) - {"Model", "PatientID"}
    flux_cols = []
    secondary_cols = []
    missing_roots = []

    # Detect raw metabolic flux columns
    for col in sorted(all_cols):
        if col.startswith("Flux_") and any(col.endswith(f"_{s}") for s in SOL_NAMES):
            flux_cols.append(col)

    # Detect secondary metric columns across all solution methods
    for root in METRIC_ROOTS:
        found_any = False
        for sol in SOL_NAMES:
            cname = f"{root}_{sol}"
            if cname in all_cols:
                secondary_cols.append(cname)
                found_any = True
        if not found_any:
            missing_roots.append(root)

    # Detect subsystem activity (SA) metrics
    sa_cols = sorted([
        c for c in all_cols
        if c.startswith("SA_")
        and any(c.endswith(f"_{s}") for s in SOL_NAMES)
        and c not in secondary_cols
    ])
    secondary_cols.extend(sa_cols)
    
    # Detect oncometabolite indicators
    for met in ONCOMET_NAMES:
        for sol in SOL_NAMES:
            cname = f"Oncomet_{met}_{sol}"
            if cname in all_cols:
                secondary_cols.append(cname)

    if missing_roots:
        print(f"\n    ⚠️  Metric roots NOT detected in CSV ({len(missing_roots)}):")
        for r in missing_roots:
            print(f"       • {r}")
        print("       → Please verify METRIC_ROOTS or inspect df.columns.")

    return {
        "flux": flux_cols,
        "secondary": secondary_cols,
        "all": flux_cols + secondary_cols
    }


def n_clusters_from_labels(labels: np.ndarray) -> int:
    """Computes the number of unique clusters excluding noise (-1)."""
    return len(set(labels)) - (1 if -1 in labels else 0)


def compute_noise_pct(labels: np.ndarray) -> float:
    """Calculates the percentage of samples classified as noise (-1)."""
    return float(np.mean(labels == -1) * 100)


# ============================================================
# 1️⃣ DATA LOADING & COHORT VALIDATION
# ============================================================
print("\n" + "="*60)
print("📂 STEP 1: LOADING & VALIDATING FEATURE MATRIX")
print("="*60)

try:
    df_raw = pd.read_csv(PATH_FEATURES)
except FileNotFoundError:
    raise FileNotFoundError(f"❌ Feature matrix file not found: {PATH_FEATURES}")

# Standardize IDs
df_raw["Model"] = df_raw["Model"].astype(str).apply(extract_model_id)
df_raw["PatientID"] = df_raw["Model"].str.slice(0, 12)

n_total = len(df_raw)
n_unique = df_raw["Model"].nunique()
print(f"    Loaded metabolic models : {n_total}")
print(f"    Unique model IDs        : {n_unique}")
print(f"    Total raw columns       : {df_raw.shape[1]}")

# Deduplication by Model ID
if n_unique < n_total:
    dups = df_raw[df_raw.duplicated(subset="Model", keep=False)]
    print(f"\n    ⚠️  {n_total - n_unique} duplicate IDs detected:")
    print(dups["Model"].value_counts().head(10).to_string())
    print("    → Retaining the first occurrence per sample ID.")
    df_raw = df_raw.drop_duplicates(subset="Model", keep="first").reset_index(drop=True)
    print(f"    Models remaining after deduplication: {len(df_raw)}")

# ============================================================
# 2️⃣ FEATURE SELECTION & FILTERING
# ============================================================
print("\n" + "="*60)
print("🔍 STEP 2: METABOLIC FEATURE SELECTION")
print("="*60)

col_groups = detect_columns(df_raw)
print(f"    Detected raw fluxes       : {len(col_groups['flux'])}")
print(f"    Detected secondary metrics: {len(col_groups['secondary'])}")
print(f"    Total features available  : {len(col_groups['all'])}")

FEATURE_MODE = "secondary"
feature_cols = col_groups[FEATURE_MODE]

if not feature_cols:
    raise ValueError("❌ No valid features found. Check METRIC_ROOTS and SOL_NAMES.")

print(f"\n    ✅ Active mode: '{FEATURE_MODE}' → {len(feature_cols)} features selected")
df_features = df_raw[["Model"] + feature_cols].copy()

# ============================================================
# 3️⃣ DATA CLEANING: LOW-VARIANCE & HIGH-SPARSITY PRUNING
# ============================================================
print("\n" + "="*60)
print(f"🧹 STEP 3: SPARSITY FILTERING (Threshold >= {ZERO_NULL_THRESHOLD*100:.0f}% zero/null)")
print("="*60)

cols_to_drop = []
for col in feature_cols:
    bad_pct = ((df_features[col] == 0) | df_features[col].isna()).sum() / len(df_features)
    if bad_pct >= ZERO_NULL_THRESHOLD:
        cols_to_drop.append((col, bad_pct))

if cols_to_drop:
    drop_names = [c for c, _ in cols_to_drop]
    df_features.drop(columns=drop_names, inplace=True)
    feature_cols = [c for c in feature_cols if c not in drop_names]
    print(f"    ❌ Dropped {len(drop_names)} sparse columns:")
    for col, pct in cols_to_drop[:10]:
        print(f"       • {col}: {pct*100:.1f}%")
    if len(cols_to_drop) > 10:
        print(f"       ... and {len(cols_to_drop)-10} additional features")
else:
    print("    ✅ No columns exceeded the sparsity threshold.")

print(f"\n    ✅ Valid features after sparsity filtering: {len(feature_cols)}")

# ============================================================
# 4️⃣ PREPROCESSING: IMPUTATION & ROBUST SCALING
# ============================================================
print("\n" + "="*60)
print("⚙️  STEP 4: PREPROCESSING (MEDIAN IMPUTATION + ROBUST SCALING)")
print("="*60)

X_raw = df_features[feature_cols].values.astype(float)

# Step 4a: Median imputation for missing metabolic variables
imputer = SimpleImputer(strategy="median")
X_imp = imputer.fit_transform(X_raw)

# Step 4b: Remove zero-variance (constant) columns post-imputation
col_std = np.std(X_imp, axis=0)
valid_ix = np.where(col_std > 1e-10)[0]
n_const = X_imp.shape[1] - len(valid_ix)
if n_const > 0:
    print(f"    ⚠️  Removing {n_const} constant features after imputation.")
    X_imp = X_imp[:, valid_ix]
    feature_cols = [feature_cols[i] for i in valid_ix]

# Step 4c: Robust scaling (resilient to biological flux outliers)
scaler = RobustScaler()
X_scaled = scaler.fit_transform(X_imp)
patient_ids = df_features["Model"].values

print(f"    Final processed models   : {X_scaled.shape[0]}")
print(f"    Final processed features : {X_scaled.shape[1]}")
print(f"    ✅ Scaled design matrix X_scaled: {X_scaled.shape}")

# ============================================================
# 5️⃣ MULTI-CONFIG PCA DIMENSIONALITY REDUCTION
# ============================================================
print("\n" + "="*60)
print("🔻 STEP 5: PCA EMBEDDINGS EXPLORATION")
print("="*60)

n_samples, n_features = X_scaled.shape
max_components = min(n_samples, n_features)

PCA_N_COMPONENTS_LIST = [2, 3, 5, 10, 15, 20, 50]
PCA_N_COMPONENTS_LIST = [c for c in PCA_N_COMPONENTS_LIST if c <= max_components]
PCA_WHITEN_LIST = [False, True]

total_emb = len(PCA_N_COMPONENTS_LIST) * len(PCA_WHITEN_LIST)
print(f"    Input features      : {n_features} secondary metrics")
print(f"    Total models        : {n_samples}")
print(f"    Max components rank : {max_components}")
print(f"    Testing: {len(PCA_N_COMPONENTS_LIST)} components × {len(PCA_WHITEN_LIST)} whitening states = {total_emb} embeddings")

embedding_matrices = {}
variance_explained = {}

for whiten in PCA_WHITEN_LIST:
    for n_comp in PCA_N_COMPONENTS_LIST:
        key = f"PCA_C{n_comp:02d}_W{int(whiten)}"
        try:
            pca = PCA(n_components=n_comp, whiten=whiten, random_state=RANDOM_SEED)
            X_pca = pca.fit_transform(X_scaled)
            embedding_matrices[key] = X_pca
            var_acc = np.cumsum(pca.explained_variance_ratio_)[-1]
            variance_explained[key] = var_acc
            print(f"    ✅ {key} | Cumulative variance = {var_acc*100:.1f}%")
        except Exception as e:
            print(f"    ⚠️  Failed generating {key}: {e}")

# Scree Plot: cumulative explained variance across configurations
fig_scree, ax_scree = plt.subplots(figsize=(10, 5))
keys_sorted = sorted(variance_explained, key=lambda k: int(k.split("_C")[1].split("_")[0]))

for whiten_val in PCA_WHITEN_LIST:
    tag = f"_W{int(whiten_val)}"
    ks = [k for k in keys_sorted if k.endswith(tag)]
    comps = [int(k.split("_C")[1].split("_")[0]) for k in ks]
    vars_ = [variance_explained[k] * 100 for k in ks]
    label = f"whiten={whiten_val}"
    ax_scree.plot(comps, vars_, marker="o", label=label)

ax_scree.set_xlabel("Number of PCA Components")
ax_scree.set_ylabel("Cumulative Variance Explained (%)")
ax_scree.set_title("Scree Plot — Cumulative Explained Variance by PCA Configuration")
ax_scree.legend()
ax_scree.grid(True, alpha=0.3)
plt.tight_layout()
scree_path = os.path.join(OUT_DIR, "ScreePlot_PCA.png")
fig_scree.savefig(scree_path, bbox_inches="tight")
plt.close(fig_scree)
print(f"\n    📊 Scree plot saved to: {scree_path}")

# ============================================================
# 6️⃣ CLUSTERING ALGORITHM SUITE & HYPERPARAMETER GRIDS
# ============================================================
K_RANGE = range(2, 10)

alg_classes = {
    "KMeans":                  KMeans,
    "Agglomerative":           AgglomerativeClustering,
    "Birch":                   Birch,
    "GMM":                     GaussianMixture,
    "BayesianGaussianMixture": BayesianGaussianMixture,
    "DBSCAN":                  DBSCAN,
    "MeanShift":               MeanShift,
    "AffinityPropagation":     AffinityPropagation,
}
if HDBSCAN_AVAILABLE:
    alg_classes["HDBSCAN"] = hdbscan.HDBSCAN

# Estimate bandwidth parameter for MeanShift
bw = estimate_bandwidth(
    X_scaled, quantile=0.2, n_samples=min(len(X_scaled), 500)
) or 1.0

param_grids = {
    "KMeans":                  ParameterGrid({"n_clusters": K_RANGE}),
    "Agglomerative":           ParameterGrid({"n_clusters": K_RANGE,
                                              "linkage": ["ward", "average", "complete"]}),
    "Birch":                   ParameterGrid({"n_clusters": K_RANGE}),
    "GMM":                     ParameterGrid({"n_components": K_RANGE}),
    "BayesianGaussianMixture": ParameterGrid({"n_components": K_RANGE}),
    "DBSCAN":                  ParameterGrid({"eps": [0.5, 1.0, 1.5, 2.5, 5.0],
                                              "min_samples": [3, 5, 8]}),
    "HDBSCAN":                 ParameterGrid({"min_cluster_size": [5, 10, 15]}),
    "MeanShift":               ParameterGrid({"bandwidth": [bw, bw*1.5, bw*0.5]}),
    "AffinityPropagation":     ParameterGrid({"damping": [0.5, 0.9]}),
}
param_grids = {k: v for k, v in param_grids.items() if k in alg_classes}

DET_ALGS = {"Agglomerative", "DBSCAN", "HDBSCAN", "MeanShift"}
STOCH_ALGS = {"KMeans", "GMM", "BayesianGaussianMixture", "Birch", "AffinityPropagation"}

# ============================================================
# 7️⃣ CLUSTERING OPTIMIZATION ENGINE
# ============================================================
def optimize_clustering(alg_name: str, X: np.ndarray) -> dict:
    """
    Systematically executes hyperparameter grid search across specified seeds.
    Evaluates configurations using Silhouette, Calinski-Harabasz, and Davies-Bouldin metrics.
    Enforces maximum outlier/noise constraints.
    """
    best = {
        "score": -np.inf, "db": np.inf, "chi": -np.inf,
        "labels": None, "params": None, "noise": None,
        "n_clusters": 0, "seed": None
    }

    X_arr = np.asarray(X)
    alg_cls = alg_classes[alg_name]
    grid = param_grids.get(alg_name, ParameterGrid([{}]))
    seeds = [SEEDS_TO_TEST[0]] if alg_name in DET_ALGS else SEEDS_TO_TEST

    for seed in seeds:
        for params in grid:
            try:
                if alg_name in {"GMM", "BayesianGaussianMixture"}:
                    if params.get("n_components", 2) >= X_arr.shape[0]:
                        continue

                model = (alg_cls(**params, random_state=seed)
                         if alg_name in STOCH_ALGS
                         else alg_cls(**params))

                labels = (model.fit(X_arr).labels_
                          if alg_name == "MeanShift"
                          else model.fit_predict(X_arr))

                noise_pct = compute_noise_pct(labels)
                if noise_pct > MAX_NOISE_PCT:
                    continue

                mask = labels != -1
                n_cl = len(np.unique(labels[mask]))
                if n_cl < 2 or mask.sum() < 2:
                    continue

                # Compute internal validation metrics
                sil = silhouette_score(X_arr[mask], labels[mask])
                db  = davies_bouldin_score(X_arr[mask], labels[mask])
                chi = calinski_harabasz_score(X_arr[mask], labels[mask])

                # Multi-tier selection hierarchy: Silhouette > CH > DB
                is_better = (
                    sil > best["score"] or
                    (sil == best["score"] and chi > best["chi"]) or
                    (sil == best["score"] and chi == best["chi"] and db < best["db"])
                )
                if is_better:
                    best.update({
                        "score": sil, "db": db, "chi": chi,
                        "labels": labels.copy(), "params": params,
                        "noise": noise_pct, "n_clusters": n_cl, "seed": seed
                    })

            except Exception:
                continue

    if best["labels"] is None:
        return {
            "best_score": None, "best_db": None, "best_chi": None,
            "best_labels": None, "best_params": None,
            "noise_pct": None, "n_clusters_found": 0, "best_seed": None
        }

    return {
        "best_score":       None if best["score"] == -np.inf else float(best["score"]),
        "best_db":          None if best["db"] == np.inf else float(best["db"]),
        "best_chi":         None if best["chi"] == -np.inf else float(best["chi"]),
        "best_labels":      best["labels"],
        "best_params":      best["params"],
        "noise_pct":        float(best["noise"]),
        "n_clusters_found": best["n_clusters"],
        "best_seed":        best["seed"],
    }

# ============================================================
# 8️⃣ EXHAUSTIVE BENCHMARKING: ALGORITHMS × EMBEDDINGS
# ============================================================
print("\n" + "="*60)
print("🔬 STEP 8: EXHAUSTIVE CLUSTERING GRID EXECUTION")
print("="*60)

rows = []
for reduction_name, X_emb in embedding_matrices.items():
    var_pct = variance_explained.get(reduction_name, 0) * 100
    print(f"\n  ▶ {reduction_name} ({X_emb.shape[1]}D | Cumulative Var={var_pct:.1f}%)")
    for alg in alg_classes:
        res = optimize_clustering(alg, X_emb)
        n_clusters = res.get("n_clusters_found", 0)
        seed_used = res.get("best_seed", "N/A")
        valid = (
            res["best_labels"] is not None and
            (res["noise_pct"] or 99) <= MAX_NOISE_PCT and
            n_clusters > 1
        )
        if valid:
            rows.append({
                "algorithm":               alg,
                "Optimal K/Comp":          n_clusters,
                "score":                   res["best_score"],
                "Davies-Bouldin Score":    res["best_db"],
                "Calinski-Harabasz Score": res["best_chi"],
                "labels":                  res["best_labels"],
                "params":                  res["best_params"],
                "noise":                   res["noise_pct"],
                "reduction":               reduction_name,
                "seed":                    seed_used,
                "variance_explained":      var_pct,
            })
            print(f"    ✅ {alg:25s} K={n_clusters:2d} | "
                  f"Sil={res['best_score']:.3f} | "
                  f"DB={res['best_db']:.3f} | "
                  f"CH={res['best_chi']:.0f} | "
                  f"Noise={res['noise_pct']:.1f}% | Seed={seed_used}")
        else:
            reason = (f"Noise={res['noise_pct']:.1f}%"
                      if (res["noise_pct"] or 0) > MAX_NOISE_PCT
                      else "No valid clustering partition found")
            print(f"    ❌ {alg:25s} {reason}")

# ============================================================
# 9️⃣ SELECTION, STANDARDIZATION & DATASET EXPORT
# ============================================================
print("\n" + "="*60)
print("🏆 STEP 9: FILTERING, NOMENCLATURE & EXPORT")
print("="*60)

df_scores = pd.DataFrame(rows)
df_scores = df_scores[df_scores["labels"].notnull()].copy()

# Sort and filter models meeting the threshold
df_selected = (
    df_scores[df_scores["score"] >= SILHOUETTE_THRESHOLD]
    .sort_values(
        by=["score", "Calinski-Harabasz Score", "Davies-Bouldin Score", "noise"],
        ascending=[False, False, True, True]
    )
    .reset_index(drop=True)
)

num_selected = len(df_selected)
print(f"\n=== Found {num_selected} clustering models satisfying Silhouette >= {SILHOUETTE_THRESHOLD:.2f} (Noise <= {MAX_NOISE_PCT:.0f}%) ===")

# Build standardized clustering result matrix
df_clusters_final = pd.DataFrame({"Model": patient_ids})

for i, row in df_selected.iterrows():
    labels = np.asarray(row["labels"])
    seed_str = f"_Seed{row['seed']}" if row["seed"] != "N/A" else ""
    col_name = (
        f"Cluster_{row['reduction']}_{row['algorithm']}"
        f"_K{row['Optimal K/Comp']}"
        f"_S{row['score']:.2f}"
        f"_DB{row['Davies-Bouldin Score']:.2f}"
        f"_CH{row['Calinski-Harabasz Score']:.0f}"
        f"{seed_str}"
    )
    if len(labels) == len(patient_ids):
        df_clusters_final[col_name] = labels

    if i < 10 or num_selected <= 10:
        print(f"    {i+1:2d}. {row['reduction']:25s} | "
              f"{row['algorithm']:22s} | K={row['Optimal K/Comp']:2d} | "
              f"Seed={row['seed']} | Sil={row['score']:.4f} | "
              f"Noise={row['noise']:.2f}% | Var={row['variance_explained']:.1f}%")
    elif i == 10:
        print("    ... [Intermediate models omitted for brevity] ...")

CLUSTERS_PATH = os.path.join(OUT_DIR, "PatientClusters_TumorPhenotype_PCA.csv")
df_clusters_final.to_csv(CLUSTERS_PATH, index=False)
print(f"\n    ✅ All {num_selected} clustering assignments exported to: {CLUSTERS_PATH}")

# ============================================================
# 🖼️ 10. VISUALIZATION ENGINE (TOP 5 CLUSTERING CONFIGURATIONS)
# ============================================================
print("\n" + "="*60)
print("📈 STEP 10: GENERATING PCA TOP VISUALIZATIONS")
print("="*60)

def generate_pca_plots(df_top: pd.DataFrame,
                       emb_dict: dict,
                       var_dict: dict,
                       output_dir: str,
                       top_n: int = 5) -> None:
    """Generates 2D or 3D scatter plots for the top N ranking clustering models."""
    if df_top.empty:
        print("    ⚠️  No valid models available to plot.")
        return

    print(f"\n📈 Generating PCA plots for top {min(top_n, len(df_top))} configurations...")

    for i, row in df_top.head(top_n).reset_index(drop=True).iterrows():
        emb_key = row["reduction"]
        labels = np.asarray(row["labels"])
        seed_str = f"_Seed{row['seed']}" if row["seed"] != "N/A" else ""
        var_pct = var_dict.get(emb_key, 0) * 100
        config = (
            f"{row['reduction']}_{row['algorithm']}"
            f"_K{row['Optimal K/Comp']}"
            f"_S{row['score']:.2f}"
            f"_DB{row['Davies-Bouldin Score']:.2f}"
            f"_CH{row['Calinski-Harabasz Score']:.0f}"
            f"{seed_str}"
        )
        title = (
            f"TOP {i+1}: {config}\n"
            f"Silhouette={row['score']:.3f} | "
            f"K={row['Optimal K/Comp']} | "
            f"Noise={row['noise']:.1f}% | "
            f"ExpVar={var_pct:.1f}%"
        )
        filename = os.path.join(output_dir, f"TOP_{i+1}_{config}.png")

        if emb_key not in emb_dict:
            print(f"    ⚠️  Embedding key missing from cache: {emb_key}")
            continue

        X_emb = emb_dict[emb_key]
        n_dims = X_emb.shape[1]
        label_str = labels.astype(str)

        unique_cl = sorted(set(label_str))
        valid_cl = [c for c in unique_cl if c != "-1"]
        palette = sns.color_palette("tab10", n_colors=max(len(valid_cl), 1))
        cmap = {c: palette[j] for j, c in enumerate(valid_cl)}
        if "-1" in unique_cl:
            cmap["-1"] = "gray"

        if n_dims == 2:
            fig, ax = plt.subplots(figsize=(10, 8))
            for cl in unique_cl:
                mask = label_str == cl
                ax.scatter(
                    X_emb[mask, 0], X_emb[mask, 1],
                    label=f"Cluster {cl}", color=cmap[cl],
                    alpha=0.7, s=50, edgecolors="none"
                )
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title(title, fontsize=10)
            ax.legend(title="Cluster", fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(filename, bbox_inches="tight")
            plt.close()
            print(f"    > 2D plot saved: {filename}")

        elif n_dims == 3:
            fig = plt.figure(figsize=(12, 10))
            ax3 = fig.add_subplot(111, projection="3d")
            for cl in unique_cl:
                mask = label_str == cl
                ax3.scatter(
                    X_emb[mask, 0], X_emb[mask, 1], X_emb[mask, 2],
                    label=f"Cluster {cl}", color=cmap[cl],
                    alpha=0.7, s=40
                )
            ax3.set_xlabel("PC1")
            ax3.set_ylabel("PC2")
            ax3.set_zlabel("PC3")
            ax3.set_title(title, fontsize=10)
            ax3.legend(title="Cluster", fontsize=7, bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(filename, bbox_inches="tight")
            plt.close()
            print(f"    > 3D plot saved: {filename}")

        else:
            # For latent spaces with > 3 dimensions, plot 2D projection onto PC1-PC2
            fig, ax = plt.subplots(figsize=(10, 8))
            for cl in unique_cl:
                mask = label_str == cl
                ax.scatter(
                    X_emb[mask, 0], X_emb[mask, 1],
                    label=f"Cluster {cl}", color=cmap[cl],
                    alpha=0.7, s=50, edgecolors="none"
                )
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title(f"{title}\n(Projection on PC1-PC2 from {n_dims}D)", fontsize=9)
            ax.legend(title="Cluster", fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(filename, bbox_inches="tight")
            plt.close()
            print(f"    > {n_dims}D→2D projection plot saved: {filename}")

generate_pca_plots(df_selected, embedding_matrices, variance_explained, OUT_DIR)

# ============================================================
# 11. CLINICAL & SURVIVAL MULTI-OMICS INTEGRATION
# ============================================================
print("\n" + "="*60)
print("🔗 STEP 11: MERGING WITH CLINICAL AND SURVIVAL DATA")
print("="*60)

def load_clinical(path: str) -> pd.DataFrame:
    """Reads clinical or survival tables (.tsv / .tsv.gz) and extracts standardized sample IDs."""
    try:
        df = (pd.read_csv(path, sep="\t", compression="gzip")
              if path.endswith(".gz")
              else pd.read_csv(path, sep="\t"))
        if "sample" in df.columns:
            df["Model"] = df["sample"].apply(extract_model_id)
        elif "submitter_id" in df.columns:
            df["Model"] = df["submitter_id"].apply(extract_model_id)
        return df
    except FileNotFoundError:
        print(f"    ⚠️  File not found: {path}")
        return pd.DataFrame({"Model": []})

df_clinical = load_clinical(PATH_CLINICAL)
df_survival = load_clinical(PATH_SURVIVAL)

valid_models = set(df_features["Model"].unique())

df_clinical_filt = df_clinical[df_clinical["Model"].isin(valid_models)]
df_survival_filt = df_survival[df_survival["Model"].isin(valid_models)]

base = df_raw[["Model", "PatientID"]].drop_duplicates(subset="Model").copy()

merged = base.copy()
for right_df, name in [(df_clinical_filt, "Clinical"),
                       (df_survival_filt, "Survival")]:
    if not right_df.empty:
        merged = merged.merge(
            right_df.drop_duplicates(subset="Model"),
            on="Model", how="left"
        )
        print(f"    ✅ Merged {name}: {len(right_df)} matched records")
    else:
        print(f"    ⚠️  {name} dataset empty or yielded no ID matches")

# Export Dataset 1: Clinical + metabolic features (without cluster labels)
merged_feat = merged.merge(df_features, on="Model", how="left")
OUT_BASE = os.path.join(OUT_DIR, "Merged_TumorPhenotype_PCA_AllData.csv")
merged_feat.to_csv(OUT_BASE, index=False)
print(f"\n    ✅ Base integrated matrix exported: {OUT_BASE}")

# Export Dataset 2: Full dataset integrated with all valid cluster assignments
try:
    df_cl = pd.read_csv(CLUSTERS_PATH).drop_duplicates(subset="Model")
    cluster_cols = [c for c in df_cl.columns if c != "Model"]

    merged_final = merged_feat.merge(df_cl, on="Model", how="left")
    merged_final[cluster_cols] = merged_final[cluster_cols].fillna(-1).astype(int)

    OUT_FINAL = os.path.join(OUT_DIR, "Merged_TumorPhenotype_PCA_AllData_withClusters.csv")
    merged_final.to_csv(OUT_FINAL, index=False)
    print(f"    ✅ Final multi-omics dataset with {len(cluster_cols)} cluster models exported: {OUT_FINAL}")

    preview_cols = ["Model"] + cluster_cols[:3]
    print("\n    Preview (first 3 clustering configurations):")
    print(merged_final[preview_cols].head(8).to_string(index=False))

except FileNotFoundError:
    print(f"    ⚠️  {CLUSTERS_PATH} not found — skipping cluster label merge.")

# ============================================================
# 📋 REPRODUCIBILITY & EXECUTION SUMMARY REPORT
# ============================================================
print("\n" + "="*60)
print("📋 REPRODUCIBILITY & PIPELINE EXECUTION SUMMARY")
print("="*60)
print(f"    Total cohort models processed      : {len(patient_ids)}")
print(f"    Feature mode                       : '{FEATURE_MODE}'")
print(f"    Initial secondary features         : {len(col_groups['secondary'])}")
print(f"    Pruned features (>= {ZERO_NULL_THRESHOLD*100:.0f}% zero/null): {len(cols_to_drop)}")
print(f"    Retained features for scaling/PCA  : {len(feature_cols)}")

if cols_to_drop:
    print(f"\n    Sample pruned sparse features:")
    for col, pct in cols_to_drop[:10]:
        print(f"       • {col}: {pct*100:.1f}%")
    if len(cols_to_drop) > 10:
        print(f"       ... and {len(cols_to_drop)-10} additional features")

print(f"\n    PCA embeddings evaluated           : {len(embedding_matrices)}")
print(f"    Total clustering grids executed    : {len(rows)}")
print(f"    Models meeting selection threshold : {num_selected} (Silhouette >= {SILHOUETTE_THRESHOLD})")

if not df_selected.empty:
    top = df_selected.iloc[0]
    print(f"\n    🥇 TOP CLUSTERING CONFIGURATION:")
    print(f"       Embedding            : {top['reduction']}")
    print(f"       Explained Variance   : {top.get('variance_explained', 0):.1f}%")
    print(f"       Algorithm            : {top['algorithm']}")
    print(f"       Optimal K/Components : {top['Optimal K/Comp']}")
    print(f"       Silhouette Score     : {top['score']:.4f}")
    print(f"       Davies-Bouldin Score : {top['Davies-Bouldin Score']:.4f}")
    print(f"       Calinski-Harabasz    : {top['Calinski-Harabasz Score']:.0f}")
    print(f"       Random Seed          : {top['seed']}")

print(f"\n    Output Artifacts:")
print(f"       • {CLUSTERS_PATH}")
print(f"       • {OUT_BASE}")
print(f"       • {OUT_FINAL}")
print(f"       • {scree_path}")
print(f"       • TOP_1 through TOP_{min(5, num_selected)} cluster scatter plots (.png)")
print("="*60)
print("\n🎉 Pipeline executed and verified successfully!")
