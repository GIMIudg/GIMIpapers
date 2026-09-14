# ============================================================
# 🚀 Pipeline: Metabolic Clustering — FeatureMatrix_TumorPhenotype
# Source Data: FeatureMatrix_TumorPhenotype.csv (MATLAB COBRA output)
# Optimizations: FBA | pFBA | L1w
# Features: Secondary metabolic metrics (excluding raw fluxes)
# Dimensionality Reduction: Multi-configuration UMAP Manifold Learning
# Reproducibility Target: Full seed tracking across manifold reduction & clustering
# ============================================================

import re
import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

from sklearn.preprocessing import RobustScaler
from sklearn.impute import SimpleImputer
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

# Verify UMAP installation
try:
    import umap
    UMAP_AVAILABLE = True
except ImportError:
    raise ImportError(
        "❌ UMAP is not installed in this environment. Run:\n"
        "   pip install umap-learn"
    )

# Verify optional HDBSCAN installation
try:
    import hdbscan
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False
    print("⚠️  Warning: HDBSCAN is not installed in this environment.")

# ============================================================
# ⚙️ CONFIGURATION & DIRECTORY SETUP
# ============================================================
_SCRIPT_DIR = Path(__file__).resolve().parent
_DATA_ROOT  = _SCRIPT_DIR.parents[3] / 'Clinical_data_and_models_ids'
PATH_FEATURES = str(_DATA_ROOT / 'Metabolic_Data' / 'FeatureMatrix_TumorPhenotype_agregado.csv')
PATH_CLINICAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.clinical.tsv')
PATH_SURVIVAL = str(_DATA_ROOT / 'Clinical_Data' / 'TCGA-BRCA.survival.tsv.gz')
OUT_DIR       = "results_TumorPhenotype_UMAP_metrics_updated_noL2"
os.makedirs(OUT_DIR, exist_ok=True)

# Random seeds for stochastic algorithms and manifold projections
SEEDS_TO_TEST = [42, 123, 100]
RANDOM_SEED   = SEEDS_TO_TEST[0]
np.random.seed(RANDOM_SEED)

# Quality cutoffs and validation thresholds
ZERO_NULL_THRESHOLD  = 0.95   # Drop features with >= 95% zero or missing values
SILHOUETTE_THRESHOLD = 0.10   # Minimum acceptable silhouette coefficient
MAX_NOISE_PCT        = 5.0    # Maximum allowable noise/unassigned sample percentage

# ============================================================
# 🔧 METABOLIC FEATURE IDENTIFIERS
# ============================================================
SOL_NAMES = ["FBA", "pFBA", "L1w"]

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
    Extracts standardized 16-character TCGA barcodes (TCGA-XX-XXXX-XX).
    Falls back to leading string slice if standard regex pattern is not matched.
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
    Scans matrix columns to identify raw fluxes, computed secondary biomarkers,
    subsystem pathway activities (SA), and oncometabolite turnover variables.
    """
    all_cols       = set(df.columns) - {"Model", "PatientID"}
    flux_cols      = []
    secondary_cols = []
    missing_roots  = []

    # Detect raw metabolic flux columns
    for col in sorted(all_cols):
        if col.startswith("Flux_") and any(col.endswith(f"_{s}") for s in SOL_NAMES):
            flux_cols.append(col)

    # Detect secondary computed metabolic features
    for root in METRIC_ROOTS:
        found_any = False
        for sol in SOL_NAMES:
            cname = f"{root}_{sol}"
            if cname in all_cols:
                secondary_cols.append(cname)
                found_any = True
        if not found_any:
            missing_roots.append(root)

    # Automatically detect Subsystem Activity (SA) metrics
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
        print(f"\n    ⚠️  Metric roots not found in CSV ({len(missing_roots)}):")
        for r in missing_roots:
            print(f"       • {r}")
        print("       → Please check METRIC_ROOTS definitions or source headers.")

    return {
        "flux": flux_cols,
        "secondary": secondary_cols,
        "all": flux_cols + secondary_cols
    }


def compute_noise_pct(labels: np.ndarray) -> float:
    """Computes percentage of samples classified as unassigned or noise (-1)."""
    return float(np.mean(labels == -1) * 100)


# ============================================================
# 1️⃣ DATA INGESTION & COHORT VALIDATION
# ============================================================
print("\n" + "="*60)
print("📂 STEP 1: INGESTING AND VALIDATING FEATURE MATRIX")
print("="*60)

try:
    df_raw = pd.read_csv(PATH_FEATURES)
except FileNotFoundError:
    raise FileNotFoundError(f"❌ Feature file not found at: {PATH_FEATURES}")

# Normalize patient and sample IDs
df_raw["Model"]     = df_raw["Model"].astype(str).apply(extract_model_id)
df_raw["PatientID"] = df_raw["Model"].str.slice(0, 12)

n_total  = len(df_raw)
n_unique = df_raw["Model"].nunique()
print(f"    Loaded metabolic models : {n_total}")
print(f"    Unique model IDs        : {n_unique}")
print(f"    Total input columns     : {df_raw.shape[1]}")

if n_unique < n_total:
    dups = df_raw[df_raw.duplicated(subset="Model", keep=False)]
    print(f"\n    ⚠️  {n_total - n_unique} duplicate model IDs detected:")
    print(dups["Model"].value_counts().head(10).to_string())
    print("    → Retaining the first occurrence per unique identifier.")
    df_raw = df_raw.drop_duplicates(subset="Model", keep="first").reset_index(drop=True)
    print(f"    Cohort size post-deduplication: {len(df_raw)}")

flux_example = [c for c in df_raw.columns if c.startswith("Flux_")][:3]
sec_example  = [c for c in df_raw.columns
                if not c.startswith("Flux_") and c not in ("Model", "PatientID")][:6]
print(f"\n    Sample flux features     : {flux_example}")
print(f"    Sample secondary metrics : {sec_example}")

# ============================================================
# 2️⃣ FEATURE SELECTION
# ============================================================
print("\n" + "="*60)
print("🔍 STEP 2: METABOLIC FEATURE SELECTION")
print("="*60)

col_groups = detect_columns(df_raw)
print(f"    Detected raw fluxes       : {len(col_groups['flux'])}")
print(f"    Detected secondary metrics: {len(col_groups['secondary'])}")
print(f"    Total available features  : {len(col_groups['all'])}")

FEATURE_MODE = "secondary"
feature_cols = col_groups[FEATURE_MODE]

if not feature_cols:
    raise ValueError("❌ No valid features found. Check METRIC_ROOTS and SOL_NAMES.")

print(f"\n    ✅ Selected mode: '{FEATURE_MODE}' → {len(feature_cols)} features retained")
print(f"    Head sample: {feature_cols[:4]} ...")

df_features = df_raw[["Model"] + feature_cols].copy()

# ============================================================
# 3️⃣ QUALITY CONTROL: HIGH-SPARSITY PRUNING
# ============================================================
print("\n" + "="*60)
print(f"🧹 STEP 3: SPARSITY FILTERING (Threshold >= {ZERO_NULL_THRESHOLD*100:.0f}% zeros/nulls)")
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
    print(f"    ❌ Pruned {len(drop_names)} high-sparsity columns:")
    for col, pct in cols_to_drop[:10]:
        print(f"       • {col}: {pct*100:.1f}%")
    if len(cols_to_drop) > 10:
        print(f"       ... and {len(cols_to_drop)-10} additional columns")
else:
    print("    ✅ Zero columns exceeded the sparsity threshold.")

print(f"\n    ✅ Features remaining after sparsity cleaning: {len(feature_cols)}")

# ============================================================
# 4️⃣ PREPROCESSING: IMPUTATION & ROBUST SCALING
# ============================================================
print("\n" + "="*60)
print("⚙️  STEP 4: DATA PREPROCESSING (MEDIAN IMPUTATION + ROBUST SCALING)")
print("="*60)

X_raw = df_features[feature_cols].values.astype(float)

# Impute missing values with column-wise medians
imputer = SimpleImputer(strategy="median")
X_imp   = imputer.fit_transform(X_raw)

# Drop invariant (zero-variance) columns post-imputation
col_std  = np.std(X_imp, axis=0)
valid_ix = np.where(col_std > 1e-10)[0]
n_const  = X_imp.shape[1] - len(valid_ix)
if n_const > 0:
    print(f"    ⚠️  Removing {n_const} constant columns post-imputation")
    X_imp        = X_imp[:, valid_ix]
    feature_cols = [feature_cols[i] for i in valid_ix]

# Scale features using median and interquartile range (IQR) to reduce outlier bias
scaler   = RobustScaler()
X_scaled = scaler.fit_transform(X_imp)

patient_ids = df_features["Model"].values

print(f"    Final processed models   : {X_scaled.shape[0]}")
print(f"    Final processed features : {X_scaled.shape[1]}")
print(f"    ✅ Scaled design matrix X_scaled: {X_scaled.shape}")

# ============================================================
# 5️⃣ MULTI-CONFIGURATION UMAP MANIFOLD LEARNING
#
#    UMAP (Uniform Manifold Approximation and Projection) recovers
#    non-linear biological manifolds while balancing local neighborhood
#    connectivity and global cluster topology.
#
#    Explored hyperparameter combinations:
#      • n_neighbors  : 5, 15, 30, 50
#        Balances local detail (low) vs broader continuum structure (high).
#      • min_dist     : 0.0, 0.1, 0.5
#        Controls packing density (0.0 optimizes clustering compactness).
#      • n_components : 2, 3
#        Latent manifold projection dimension.
#      • metric       : euclidean, cosine
#        Distance metric in scaled metabolic metric space.
#      • seeds        : SEEDS_TO_TEST (42, 123, 100)
#        Evaluates manifold projection stability across stochastic runs.
#
#    Total configurations = 4 * 3 * 2 * 2 * 3 seeds = 144 embeddings
# ============================================================
print("\n" + "="*60)
print("🔻 STEP 5: UMAP MANIFOLD PROJECTIONS (SECONDARY METRICS)")
print("="*60)

UMAP_N_NEIGHBORS  = [5, 15, 30, 50]
UMAP_MIN_DIST     = [0.0, 0.1, 0.5]
UMAP_N_COMPONENTS = [2, 3]
UMAP_METRICS      = ["euclidean", "cosine"]

# Ensure n_neighbors does not exceed cohort sample count
UMAP_N_NEIGHBORS = [n for n in UMAP_N_NEIGHBORS if n < X_scaled.shape[0]]

total_emb = (len(UMAP_N_NEIGHBORS) * len(UMAP_MIN_DIST) *
             len(UMAP_N_COMPONENTS) * len(UMAP_METRICS) * len(SEEDS_TO_TEST))

print(f"    Input features      : {X_scaled.shape[1]} secondary metrics")
print(f"    Cohort samples      : {X_scaled.shape[0]}")
print(f"    Grid specifications : {len(UMAP_N_NEIGHBORS)} n_neighbors × "
      f"{len(UMAP_MIN_DIST)} min_dist × "
      f"{len(UMAP_N_COMPONENTS)} dims × "
      f"{len(UMAP_METRICS)} metrics × "
      f"{len(SEEDS_TO_TEST)} seeds = {total_emb} total embeddings")

embedding_matrices = {}

for seed in SEEDS_TO_TEST:
    np.random.seed(seed)
    for n_comp in UMAP_N_COMPONENTS:
        for n_neigh in UMAP_N_NEIGHBORS:
            for min_d in UMAP_MIN_DIST:
                for metric in UMAP_METRICS:
                    key = (f"UMAP_N{n_neigh:02d}_D{str(min_d).replace('.','')}"
                           f"_C{n_comp}_M{metric}_S{seed}")
                    try:
                        reducer = umap.UMAP(
                            n_neighbors  = n_neigh,
                            min_dist     = min_d,
                            n_components = n_comp,
                            metric       = metric,
                            random_state = seed,
                            n_jobs       = -1,
                            verbose      = False,
                        )
                        X_umap = reducer.fit_transform(X_scaled)
                        embedding_matrices[key] = X_umap
                    except Exception as e:
                        print(f"    ⚠️  Failed generating {key}: {e}")

    n_done = len([k for k in embedding_matrices if f"_S{seed}" in k])
    print(f"    ✅ Seed {seed} completed — {n_done} total embeddings compiled")

print(f"\n    Total generated UMAP embeddings : {len(embedding_matrices)}")
print(f"    Sample embedding keys           : {list(embedding_matrices.keys())[:4]} ...")

# ============================================================
# 6️⃣ CLUSTERING ALGORITHM SUITE & PARAMETER GRIDS
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
    "DBSCAN":                  ParameterGrid({"eps":         [0.5, 1.0, 1.5, 2.5, 5.0],
                                              "min_samples": [3, 5, 8]}),
    "HDBSCAN":                 ParameterGrid({"min_cluster_size": [5, 10, 15]}),
    "MeanShift":               ParameterGrid({"bandwidth": [bw, bw*1.5, bw*0.5]}),
    "AffinityPropagation":     ParameterGrid({"damping": [0.5, 0.9]}),
}
param_grids = {k: v for k, v in param_grids.items() if k in alg_classes}

DET_ALGS   = {"Agglomerative", "DBSCAN", "HDBSCAN", "MeanShift"}
STOCH_ALGS = {"KMeans", "GMM", "BayesianGaussianMixture", "Birch", "AffinityPropagation"}

# ============================================================
# 7️⃣ CLUSTERING OPTIMIZATION & INTERNAL VALIDATION
# ============================================================
def optimize_clustering(alg_name: str, X: np.ndarray) -> dict:
    """
    Executes a parameter search for a specified clustering algorithm on embedding matrix X.
    Evaluates partitions using Silhouette, Calinski-Harabasz, and Davies-Bouldin metrics.
    Enforces noise/outlier constraints (< MAX_NOISE_PCT).
    """
    best = {"score": -np.inf, "db": np.inf, "chi": -np.inf,
            "labels": None, "params": None, "noise": None,
            "n_clusters": 0, "seed": None}

    X_arr   = np.asarray(X)
    alg_cls = alg_classes[alg_name]
    grid    = param_grids.get(alg_name, ParameterGrid([{}]))
    seeds   = [SEEDS_TO_TEST[0]] if alg_name in DET_ALGS else SEEDS_TO_TEST

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

                sil = silhouette_score(X_arr[mask], labels[mask])
                db  = davies_bouldin_score(X_arr[mask], labels[mask])
                chi = calinski_harabasz_score(X_arr[mask], labels[mask])

                # Multi-objective criteria: Silhouette > Calinski-Harabasz > Davies-Bouldin
                is_better = (
                    sil > best["score"] or
                    (sil == best["score"] and chi > best["chi"]) or
                    (sil == best["score"] and chi == best["chi"] and db < best["db"])
                )
                if is_better:
                    best.update({"score": sil, "db": db, "chi": chi,
                                 "labels": labels.copy(), "params": params,
                                 "noise": noise_pct, "n_clusters": n_cl, "seed": seed})

            except Exception:
                continue

    if best["labels"] is None:
        return {"best_score": None, "best_db": None, "best_chi": None,
                "best_labels": None, "best_params": None,
                "noise_pct": None, "n_clusters_found": 0, "best_seed": None}

    return {
        "best_score":       None if best["score"] == -np.inf else float(best["score"]),
        "best_db":          None if best["db"]    == np.inf  else float(best["db"]),
        "best_chi":         None if best["chi"]   == -np.inf else float(best["chi"]),
        "best_labels":      best["labels"],
        "best_params":      best["params"],
        "noise_pct":        float(best["noise"]),
        "n_clusters_found": best["n_clusters"],
        "best_seed":        best["seed"],
    }

# ============================================================
# 8️⃣ EXHAUSTIVE BENCHMARKING ACROSS EMBEDDINGS
# ============================================================
print("\n" + "="*60)
print("🔬 STEP 8: SYSTEMATIC CLUSTERING GRID EXECUTION")
print("="*60)

rows = []
for reduction_name, X_emb in embedding_matrices.items():
    print(f"\n  ▶ Testing: {reduction_name} ({X_emb.shape[1]}D)")
    for alg in alg_classes:
        res        = optimize_clustering(alg, X_emb)
        n_clusters = res.get("n_clusters_found", 0)
        seed_used  = res.get("best_seed", "N/A")
        valid      = (res["best_labels"] is not None and
                      (res["noise_pct"] or 99) <= MAX_NOISE_PCT and
                      n_clusters > 1)
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
            })
            print(f"    ✅ {alg:25s} K={n_clusters:2d} | "
                  f"Sil={res['best_score']:.3f} | "
                  f"DB={res['best_db']:.3f} | "
                  f"CH={res['best_chi']:.0f} | "
                  f"Noise={res['noise_pct']:.1f}% | Seed={seed_used}")
        else:
            reason = (f"Noise={res['noise_pct']:.1f}%"
                      if (res["noise_pct"] or 0) > MAX_NOISE_PCT
                      else "No valid cluster assignment")
            print(f"    ❌ {alg:25s} {reason}")

# ============================================================
# 9️⃣ MODEL FILTERING & EXPORT OF CLUSTER PARTITIONS
# ============================================================
print("\n" + "="*60)
print("🏆 STEP 9: SELECTION, STANDARDIZED NAMING, AND EXPORT")
print("="*60)

df_scores = pd.DataFrame(rows)
df_scores = df_scores[df_scores["labels"].notnull()].copy()

df_selected = (
    df_scores[df_scores["score"] >= SILHOUETTE_THRESHOLD]
    .sort_values(
        by=["score", "Calinski-Harabasz Score", "Davies-Bouldin Score", "noise"],
        ascending=[False, False, True, True]
    )
    .reset_index(drop=True)
)

num_selected = len(df_selected)
print(f"\n=== 🏆 Found {num_selected} models meeting Silhouette >= {SILHOUETTE_THRESHOLD:.2f} (Max noise {MAX_NOISE_PCT:.0f}%) ===")

# Build final cluster matrix: one standardized column per configuration
df_clusters_final = pd.DataFrame({"Model": patient_ids})

for i, row in df_selected.iterrows():
    labels   = np.asarray(row["labels"])
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
        print(f"    {i+1:2d}. {row['reduction']:45s} | "
              f"{row['algorithm']:22s} | K={row['Optimal K/Comp']:2d} | "
              f"Seed={row['seed']} | Sil={row['score']:.4f} | "
              f"Noise={row['noise']:.2f}%")
    elif i == 10:
        print("    ... [Intermediate configurations omitted for brevity] ...")

CLUSTERS_PATH = os.path.join(OUT_DIR, "PatientClusters_TumorPhenotype_UMAP.csv")
df_clusters_final.to_csv(CLUSTERS_PATH, index=False)
print(f"\n    ✅ All {num_selected} clustering assignments exported to: {CLUSTERS_PATH}")

# ============================================================
# 🖼️ 10. VISUALIZATION (TOP 5 CLUSTERING SOLUTIONS)
# ============================================================
print("\n" + "="*60)
print("📈 STEP 10: GENERATING TOP-RANKED UMAP PROJECTION PLOTS")
print("="*60)

def generate_umap_plots(df_top: pd.DataFrame,
                        emb_dict: dict,
                        output_dir: str,
                        top_n: int = 5) -> None:
    """Generates 2D or 3D scatter plots for the top N ranking clustering models."""
    if df_top.empty:
        print("    ⚠️  No valid models available to plot.")
        return

    print(f"\n📈 Generating UMAP projection plots for top {min(top_n, len(df_top))} configurations...")

    for i, row in df_top.head(top_n).reset_index(drop=True).iterrows():
        emb_key   = row["reduction"]
        labels    = np.asarray(row["labels"])
        seed_str  = f"_Seed{row['seed']}" if row["seed"] != "N/A" else ""
        config    = (f"{row['reduction']}_{row['algorithm']}"
                     f"_K{row['Optimal K/Comp']}"
                     f"_S{row['score']:.2f}"
                     f"_DB{row['Davies-Bouldin Score']:.2f}"
                     f"_CH{row['Calinski-Harabasz Score']:.0f}"
                     f"{seed_str}")
        title     = (f"TOP {i+1}: {config}\n"
                     f"Silhouette={row['score']:.3f} | "
                     f"K={row['Optimal K/Comp']} | Noise={row['noise']:.1f}%")
        filename  = os.path.join(output_dir, f"TOP_{i+1}_{config}.png")

        if emb_key not in emb_dict:
            print(f"    ⚠️  Embedding key missing: {emb_key}")
            continue

        X_emb     = emb_dict[emb_key]
        n_dims    = X_emb.shape[1]
        label_str = labels.astype(str)

        unique_cl = sorted(set(label_str))
        valid_cl  = [c for c in unique_cl if c != "-1"]
        palette   = sns.color_palette("tab10", n_colors=max(len(valid_cl), 1))
        cmap      = {c: palette[j] for j, c in enumerate(valid_cl)}
        if "-1" in unique_cl:
            cmap["-1"] = "gray"

        if n_dims == 2:
            fig, ax = plt.subplots(figsize=(10, 8))
            for cl in unique_cl:
                mask = label_str == cl
                ax.scatter(X_emb[mask, 0], X_emb[mask, 1],
                           label=f"Cluster {cl}", color=cmap[cl],
                           alpha=0.7, s=50, edgecolors="none")
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.set_title(title, fontsize=10)
            ax.legend(title="Cluster", fontsize=8,
                      bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(filename, bbox_inches="tight")
            plt.close()
            print(f"    > 2D plot saved: {filename}")

        elif n_dims == 3:
            fig = plt.figure(figsize=(12, 10))
            ax3 = fig.add_subplot(111, projection="3d")
            for cl in unique_cl:
                mask = label_str == cl
                ax3.scatter(X_emb[mask, 0], X_emb[mask, 1], X_emb[mask, 2],
                            label=f"Cluster {cl}", color=cmap[cl],
                            alpha=0.7, s=40)
            ax3.set_xlabel("UMAP 1")
            ax3.set_ylabel("UMAP 2")
            ax3.set_zlabel("UMAP 3")
            ax3.set_title(title, fontsize=10)
            ax3.legend(title="Cluster", fontsize=7,
                       bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(filename, bbox_inches="tight")
            plt.close()
            print(f"    > 3D plot saved: {filename}")

generate_umap_plots(df_selected, embedding_matrices, OUT_DIR)

# ============================================================
# 11. MULTI-OMICS INTEGRATION: CLINICAL & SURVIVAL MERGE
# ============================================================
print("\n" + "="*60)
print("🔗 STEP 11: MERGING WITH CLINICAL AND SURVIVAL DATASETS")
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
        print(f"    ⚠️  File not found at: {path}")
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
        merged = merged.merge(right_df.drop_duplicates(subset="Model"),
                              on="Model", how="left")
        print(f"    ✅ Merged {name}: {len(right_df)} matching patient records")
    else:
        print(f"    ⚠️  {name} file yielded no matches or is empty")

# Dataset 1: Clinical + metabolic features (without cluster labels)
merged_feat = merged.merge(df_features, on="Model", how="left")

OUT_BASE = os.path.join(OUT_DIR, "Merged_TumorPhenotype_UMAP_AllData.csv")
merged_feat.to_csv(OUT_BASE, index=False)
print(f"\n    ✅ Base multi-omics dataset saved: {OUT_BASE}")

# Dataset 2: Integrated table containing all valid clusterings
try:
    df_cl        = pd.read_csv(CLUSTERS_PATH).drop_duplicates(subset="Model")
    cluster_cols = [c for c in df_cl.columns if c != "Model"]

    merged_final = merged_feat.merge(df_cl, on="Model", how="left")
    merged_final[cluster_cols] = merged_final[cluster_cols].fillna(-1).astype(int)

    OUT_FINAL = os.path.join(OUT_DIR,
                             "Merged_TumorPhenotype_UMAP_AllData_withClusters.csv")
    merged_final.to_csv(OUT_FINAL, index=False)
    print(f"    ✅ Final integrated matrix ({len(cluster_cols)} cluster models): {OUT_FINAL}")

    preview_cols = ["Model"] + cluster_cols[:3]
    print("\n    Preview (first 3 cluster configurations):")
    print(merged_final[preview_cols].head(8).to_string(index=False))

except FileNotFoundError:
    print(f"    ⚠️  {CLUSTERS_PATH} not found — skipping cluster merge")

# ============================================================
# 📋 FINAL REPORT
# ============================================================
print("\n" + "="*60)
print("📋 PIPELINE REPRODUCIBILITY & EXECUTION REPORT")
print("="*60)
print(f"    Total cohort models processed      : {len(patient_ids)}")
print(f"    Active feature extraction mode     : '{FEATURE_MODE}'")
print(f"    Initial secondary feature count    : {len(col_groups['secondary'])}")
print(f"    Pruned high-sparsity columns       : {len(cols_to_drop)}")
print(f"    Final metabolic features scaled    : {len(feature_cols)}")
if cols_to_drop:
    print(f"\n    Sample pruned features (>= {ZERO_NULL_THRESHOLD*100:.0f}% zero/null):")
    for col, pct in cols_to_drop[:10]:
        print(f"       • {col}: {pct*100:.1f}%")
    if len(cols_to_drop) > 10:
        print(f"       ... and {len(cols_to_drop)-10} more")
print(f"\n    Generated UMAP manifold embeddings : {len(embedding_matrices)}")
print(f"    Total clustering iterations        : {len(rows)}")
print(f"    Models satisfying quality cutoff   : {num_selected} (Silhouette >= {SILHOUETTE_THRESHOLD})")

if not df_selected.empty:
    top = df_selected.iloc[0]
    print(f"\n    🥇 OPTIMAL CLUSTERING SOLUTION:")
    print(f"       Embedding   : {top['reduction']}")
    print(f"       Algorithm   : {top['algorithm']}")
    print(f"       K / Comps   : {top['Optimal K/Comp']}")
    print(f"       Silhouette  : {top['score']:.4f}")
    print(f"       DB Score    : {top['Davies-Bouldin Score']:.4f}")
    print(f"       CH Score    : {top['Calinski-Harabasz Score']:.0f}")
    print(f"       Random Seed : {top['seed']}")

print(f"\n    Generated Output Files:")
print(f"       • {CLUSTERS_PATH}")
print(f"       • {OUT_BASE}")
print(f"       • {OUT_FINAL}")
print(f"       • Top-ranked projection plots (.png)")
print("="*60)
print("\n🎉 Pipeline finished successfully!")