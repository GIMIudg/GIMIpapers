# ============================================================
# 🔬 EXPLORATORY AND DESCRIPTIVE ANALYSIS
#     Radiomics — TCGA-Run-2014 (91 cases, UChicago V2010 MRI Workstation)
# ============================================================

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from scipy import stats
import re
import warnings
import logging
import sys

warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# ============================================================
# PATH AND SEED CONFIGURATION
# ============================================================
from pathlib import Path as _Path
_SCRIPT_DIR = _Path(__file__).resolve().parent
_DATA_ROOT  = _SCRIPT_DIR.parents[2] / 'Clinical_data_and_models_ids'
# Clinical data files (TCGA clinical, survival, metadata, model IDs)
PATH_BASE   = str(_DATA_ROOT / 'Radiomics_Data') + '/'
CLINICAL_DATA = str(_DATA_ROOT / 'MetaData') + '/'
# Note: the radiomics Excel file must also be located in PATH_BASE.

RADIOMICS_FILE = 'TCGA-Run-2014_91cases_features_UChicago-V2010-MRI-Workstation.xls'
RADIOMICS_PATH = PATH_BASE + RADIOMICS_FILE

SEED = 42
np.random.seed(SEED)

# Output folder (next to the script, not inside Clinical_Data)
OUTPUT_DIR = _SCRIPT_DIR / 'EDA_Radiomics_Output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def outpath(filename):
    return str(OUTPUT_DIR / filename)

# ------------------------------------------------------------
# OPTIONAL: GROUP COMPARISON SETTINGS (Section 10)
# Mann-Whitney U and Cliff's delta need a binary grouping variable
# (e.g., ER status, molecular subtype, tumor grade binarized...).
# Leave GROUP_FILE / GROUP_COL as None to skip that section.
# ------------------------------------------------------------
GROUP_FILE   = None            # e.g., PATH_BASE + 'TCGA_clinical.xlsx'  (.xlsx / .xls / .csv)
GROUP_ID_COL = 'ModelName'     # column in that file containing the TCGA patient barcode
GROUP_COL    = None            # e.g., 'ER_status'
N_BOOT       = 2000            # bootstrap resamples for the Cliff's delta 95% CI
ALPHA        = 0.05

# ============================================================
# GRAPHIC SPECIFICATIONS
# ============================================================
DPI        = 300
FONT_SIZE  = 10
FONT_TITLE = 13
FONT_AXIS  = 11

_available  = {f.name for f in fm.fontManager.ttflist}
FONT_FAMILY = 'Arial' if 'Arial' in _available else 'DejaVu Sans'
print(f"  Font: {FONT_FAMILY}" + ("  ✔" if FONT_FAMILY == 'Arial' else "  ⚠ falling back to DejaVu Sans"))

plt.rcParams.update({
    'font.family':     FONT_FAMILY,
    'font.size':       FONT_SIZE,
    'axes.titlesize':  FONT_TITLE,
    'axes.labelsize':  FONT_AXIS,
    'xtick.labelsize': FONT_SIZE - 1,
    'ytick.labelsize': FONT_SIZE - 1,
    'legend.fontsize': FONT_SIZE,
    'figure.dpi':      DPI,
    'savefig.dpi':     DPI,
})

def save_fig(fig, filename):
    path = outpath(filename)
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✅  {filename}")

# ============================================================
# 1️⃣  DATA LOADING
# ============================================================

print(f"\n📂 Loading: {RADIOMICS_FILE}")
try:
    df = pd.read_excel(RADIOMICS_PATH, engine='xlrd')
except ImportError:
    sys.exit("❌ The 'xlrd' package is required to read .xls files. Install it with: pip install xlrd")
except FileNotFoundError:
    sys.exit(f"❌ File not found at: {RADIOMICS_PATH}\n"
             f"   Check that it is located in: {PATH_BASE}")
except Exception as e:
    sys.exit(f"❌ Error while reading the file: {e}")

ID_COL = 'Lesion Name'
if ID_COL not in df.columns:
    sys.exit(f"❌ Expected ID column '{ID_COL}' was not found. "
             f"Available columns: {list(df.columns)}")

print(f"✅ Data loaded: {df.shape[0]} lesions × {df.shape[1]-1} features\n")

# ── Normalize ID to patient level (12-character TCGA barcode) ───────────────
def extract_tcga_barcode(lesion_name):
    m = re.search(r'(TCGA-\w{2}-\w{4})', str(lesion_name))
    return m.group(1) if m else None

df['ModelName']       = df[ID_COL].apply(extract_tcga_barcode)
n_no_match            = df['ModelName'].isna().sum()
if n_no_match > 0:
    print(f"⚠️  {n_no_match} 'Lesion Name' value(s) did not match the TCGA-XX-XXXX pattern:")
    print(df.loc[df['ModelName'].isna(), ID_COL].tolist())

# Reorder: ModelName first, then original Lesion Name, then features
feature_cols = [c for c in df.columns if c not in (ID_COL, 'ModelName')]
df = df[['ModelName', ID_COL] + feature_cols]

# ── Detect multiple lesions per patient (same ModelName) ────────────────────
dup_patients = df['ModelName'][df['ModelName'].duplicated(keep=False)].dropna().unique()
if len(dup_patients) > 0:
    print(f"\n⚠️  {len(dup_patients)} patient(s) with more than one recorded lesion:")
    print(df[df['ModelName'].isin(dup_patients)][['ModelName', ID_COL]].to_string(index=False))
    print("   (decide whether to average them, keep the index lesion, or keep them separate)")
else:
    print("\n✅ Each patient has exactly one recorded lesion (no duplicates).")

df.to_csv(outpath('radiomics_with_ModelName.csv'), index=False)
print("✅ 'radiomics_with_ModelName.csv' generated (with normalized patient ID).")

# ============================================================
# 2️⃣  FEATURE CATEGORIZATION (based on the code in parentheses)
# ============================================================

CATEGORY_NAMES = {
    'K': 'Enhancement kinetics (Kinetic curve)',
    'E': 'Enhancement variance (Enhancement-variance)',
    'T': 'Texture (GLCM Texture)',
    'G': 'Geometry (Geometry/Shape)',
    'M': 'Margin (Margin)',
    'S': 'Size (Size)',
}

def get_feature_code(col):
    m = re.search(r'\(([A-Z])(\d+)\)\s*$', col)
    return (m.group(1), f"{m.group(1)}{m.group(2)}") if m else (None, None)

feature_groups = {}
for col in feature_cols:
    letter, code = get_feature_code(col)
    if letter:
        feature_groups.setdefault(letter, []).append(col)

print(f"\n📂 Feature categories detected:")
for letter, cols in feature_groups.items():
    name = CATEGORY_NAMES.get(letter, letter)
    print(f"   {letter} — {name}: {len(cols)} variables")

# ============================================================
# 3️⃣  DESCRIPTIVE STATISTICS
# ============================================================

print("\n⏳ Computing descriptive statistics...")

desc_rows = []
for col in feature_cols:
    vals = df[col].dropna()
    letter, code = get_feature_code(col)
    desc_rows.append({
        'Feature':       col,
        'Category':      CATEGORY_NAMES.get(letter, 'Other'),
        'Code':          code,
        'N':             len(vals),
        'N_missing':     df[col].isna().sum(),
        'Mean':          vals.mean(),
        'Median':        vals.median(),
        'SD':            vals.std(),
        'Min':           vals.min(),
        'Max':           vals.max(),
        'Q1':            vals.quantile(0.25),
        'Q3':            vals.quantile(0.75),
        'IQR':           vals.quantile(0.75) - vals.quantile(0.25),
        'CV_%':          (vals.std() / vals.mean() * 100) if vals.mean() != 0 else np.nan,
        'Skewness':      stats.skew(vals),
        'Kurtosis':      stats.kurtosis(vals),
        'Shapiro_p':     stats.shapiro(vals)[1] if len(vals) >= 3 else np.nan,
    })

df_desc = pd.DataFrame(desc_rows)
df_desc['Normal_p05'] = np.where(df_desc['Shapiro_p'] >= 0.05, 'Yes (p≥0.05)', 'No (p<0.05)')
df_desc.to_csv(outpath('descriptive_statistics.csv'), index=False)
print("✅ 'descriptive_statistics.csv' generated.")

n_normal = (df_desc['Shapiro_p'] >= 0.05).sum()
print(f"\n📊 Normality summary (Shapiro-Wilk, α=0.05):")
print(f"   {n_normal}/{len(df_desc)} features consistent with a normal distribution.")
print(f"   {len(df_desc) - n_normal}/{len(df_desc)} show significant deviation from normality")
print(f"   → consider non-parametric tests / transformations for the latter.")

# ============================================================
# 4️⃣  MISSING VALUES
# ============================================================

missing_summary = df[feature_cols].isna().sum()
missing_summary = missing_summary[missing_summary > 0]

if len(missing_summary) > 0:
    print(f"\n⚠️  Features with missing values:")
    print(missing_summary.to_string())
else:
    print(f"\n✅ No missing values in any of the {len(feature_cols)} features.")

# ============================================================
# 5️⃣  OUTLIER DETECTION (IQR method)
# ============================================================

print("\n⏳ Detecting outliers (1.5×IQR rule)...")

outlier_rows = []
for col in feature_cols:
    vals = df[col].dropna()
    q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
    iqr = q3 - q1
    lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    mask = (vals < lower) | (vals > upper)
    outlier_ids = df.loc[vals[mask].index, 'ModelName'].tolist()
    outlier_rows.append({
        'Feature':      col,
        'N_outliers':   mask.sum(),
        'Pct_outliers': round(mask.sum() / len(vals) * 100, 1),
        'Lower_limit':  lower,
        'Upper_limit':  upper,
        'Outlier_patients': ', '.join(outlier_ids) if len(outlier_ids) <= 8 else f"{len(outlier_ids)} patients (see detail)"
    })

df_outliers = pd.DataFrame(outlier_rows).sort_values('N_outliers', ascending=False)
df_outliers.to_csv(outpath('outlier_detection.csv'), index=False)
print("✅ 'outlier_detection.csv' generated.")
print(f"\n📊 Top 5 features with the most outliers:")
print(df_outliers.head(5)[['Feature', 'N_outliers', 'Pct_outliers']].to_string(index=False))

# ============================================================
# 6️⃣  DISTRIBUTIONS BY CATEGORY (histogram + KDE)
# ============================================================

print("\n🚀 Generating distributions by category...")

for letter, cols in feature_groups.items():
    n_feat = len(cols)
    ncols_grid = 4
    nrows_grid = int(np.ceil(n_feat / ncols_grid))

    fig, axes = plt.subplots(nrows_grid, ncols_grid,
                              figsize=(ncols_grid * 3.2, nrows_grid * 2.8))
    axes = np.array(axes).reshape(-1)

    for i, col in enumerate(cols):
        ax = axes[i]
        vals = df[col].dropna()
        sns.histplot(vals, kde=True, ax=ax, color=sns.color_palette('Set2')[list(feature_groups.keys()).index(letter) % 8],
                     edgecolor='white', linewidth=0.5)
        short_name = re.sub(r'\s*\([A-Z]\d+\)\s*$', '', col)
        ax.set_title(short_name[:32], fontsize=FONT_SIZE - 1, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Frequency', fontsize=FONT_SIZE - 2)
        ax.tick_params(labelsize=FONT_SIZE - 2)

    for j in range(n_feat, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Distribution — {CATEGORY_NAMES.get(letter, letter)}",
                 fontsize=FONT_TITLE, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_fig(fig, f'Distribution_Category_{letter}.png')

# ============================================================
# 7️⃣  STANDARDIZED BOXPLOTS (z-score) — ALL FEATURES
#     Allows comparing dispersion/outliers across variables on
#     very different scales in a single plot per category.
# ============================================================

print("\n🚀 Generating standardized boxplots by category...")

for letter, cols in feature_groups.items():
    df_z = df[cols].apply(lambda x: (x - x.mean()) / x.std())
    df_z_long = df_z.melt(var_name='Feature', value_name='Z-score')
    df_z_long['Feature'] = df_z_long['Feature'].apply(lambda c: re.sub(r'\s*\([A-Z]\d+\)\s*$', '', c))

    fig, ax = plt.subplots(figsize=(max(6, len(cols) * 0.9), 5))
    sns.boxplot(data=df_z_long, x='Feature', y='Z-score', hue='Feature', ax=ax,
                palette='Set2', width=0.5, linewidth=1.1, legend=False)
    sns.stripplot(data=df_z_long, x='Feature', y='Z-score', ax=ax,
                  color='#333333', alpha=0.35, size=3, jitter=True)
    ax.axhline(0, color='red', lw=0.8, ls='--', alpha=0.5)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([re.sub(r'\s*\([A-Z]\d+\)\s*$', '', c) for c in cols],
                        rotation=35, ha='right', fontsize=FONT_SIZE - 1)
    ax.set_title(f"Standardized boxplot (z-score) — {CATEGORY_NAMES.get(letter, letter)}",
                 fontsize=FONT_TITLE - 1, fontweight='bold')
    ax.set_xlabel('')
    plt.tight_layout()
    save_fig(fig, f'Boxplot_Zscore_Category_{letter}.png')

# ============================================================
# 8️⃣  CORRELATION MATRIX (all features)
# ============================================================

print("\n🚀 Generating correlation matrix...")

corr_matrix = df[feature_cols].corr(method='spearman')
short_labels = [re.sub(r'\s*\(([A-Z]\d+)\)\s*$', r' [\1]', c) for c in feature_cols]

fig, ax = plt.subplots(figsize=(14, 13))
sns.heatmap(
    corr_matrix, annot=False, cmap='RdBu_r', center=0,
    vmin=-1, vmax=1, square=True, linewidths=0.3, linecolor='white',
    cbar_kws={'label': "Spearman correlation", 'shrink': 0.7},
    xticklabels=short_labels, yticklabels=short_labels, ax=ax
)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
ax.set_title(f"Correlation Matrix (Spearman) — {len(feature_cols)} Radiomic Features",
             fontsize=FONT_TITLE, fontweight='bold', pad=12)
plt.tight_layout()
save_fig(fig, 'Correlation_Matrix_Spearman.png')

# ── Highly correlated pairs (possible redundancy) ───────────────────────────
corr_pairs = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
corr_pairs = corr_pairs.stack().reset_index()
corr_pairs.columns = ['Feature_A', 'Feature_B', 'Correlation']
corr_pairs['Correlation_abs'] = corr_pairs['Correlation'].abs()
corr_pairs_high = corr_pairs[corr_pairs['Correlation_abs'] >= 0.90].sort_values(
    'Correlation_abs', ascending=False)

corr_pairs_high.to_csv(outpath('high_correlation_pairs_r090.csv'), index=False)
print(f"✅ 'high_correlation_pairs_r090.csv' generated ({len(corr_pairs_high)} pairs with |r|≥0.90).")
if len(corr_pairs_high) > 0:
    print("   Top 10 most correlated pairs (possible redundancy for modeling):")
    print(corr_pairs_high.head(10)[['Feature_A', 'Feature_B', 'Correlation']].to_string(index=False))

# ============================================================
# 9️⃣  EXPLORATORY PCA (quick multivariate view)
# ============================================================

print("\n🚀 Running exploratory PCA...")

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

X = df[feature_cols].values
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=min(10, len(feature_cols)), random_state=SEED)
pcs = pca.fit_transform(X_scaled)
var_exp = pca.explained_variance_ratio_ * 100
var_cum = np.cumsum(var_exp)

df_pca_summary = pd.DataFrame({
    'Component':           [f'PC{i+1}' for i in range(len(var_exp))],
    'Explained_variance_%': var_exp.round(2),
    'Cumulative_variance_%': var_cum.round(2),
})
df_pca_summary.to_csv(outpath('pca_explained_variance.csv'), index=False)
print("✅ 'pca_explained_variance.csv' generated.")
print(df_pca_summary.head(5).to_string(index=False))

n_pc_80 = int(np.argmax(var_cum >= 80) + 1)
print(f"\n   → {n_pc_80} components are needed to explain ≥80% of the total variance.")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].bar(range(1, len(var_exp) + 1), var_exp, color='#4C72B0', alpha=0.85, label='Individual')
axes[0].plot(range(1, len(var_exp) + 1), var_cum, color='#D55E00', marker='o', label='Cumulative')
axes[0].axhline(80, color='grey', ls='--', lw=0.8)
axes[0].set_xlabel('Principal Component')
axes[0].set_ylabel('Explained variance (%)')
axes[0].set_title('Scree Plot — PCA', fontsize=FONT_TITLE, fontweight='bold')
axes[0].legend(fontsize=FONT_SIZE - 1)

scatter = axes[1].scatter(pcs[:, 0], pcs[:, 1], c=df[feature_cols[0]], cmap='viridis',
                           alpha=0.8, edgecolors='white', linewidths=0.4, s=50)
axes[1].set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
axes[1].set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
axes[1].set_title('PC1 vs PC2 Projection\n(color = ' + re.sub(r'\s*\([A-Z]\d+\)$', '', feature_cols[0]) + ')',
                   fontsize=FONT_TITLE - 1, fontweight='bold')
cb = fig.colorbar(scatter, ax=axes[1], shrink=0.8)
cb.set_label(re.sub(r'\s*\([A-Z]\d+\)$', '', feature_cols[0]), fontsize=FONT_SIZE - 1)

plt.tight_layout()
save_fig(fig, 'PCA_scree_and_projection.png')

# ============================================================
# 🔟  NON-PARAMETRIC GROUP COMPARISON
#     Mann-Whitney U test + Cliff's delta (effect size)
#
#     Cliff's delta  δ = P(X > Y) − P(X < Y),  range [−1, 1]
#       δ = 0  → complete overlap between groups
#       δ = ±1 → no overlap (all X > all Y, or vice versa)
#     It is directly related to Mann-Whitney U:  δ = 2U/(n1·n2) − 1
#     Magnitude thresholds (Romano et al., 2006):
#       |δ| < 0.147 negligible | < 0.33 small | < 0.474 medium | ≥ 0.474 large
# ============================================================

def cliffs_delta(x, y):
    """Cliff's delta between two samples (positive → x tends to be larger than y)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    diff = x[:, None] - y[None, :]
    return (np.sum(diff > 0) - np.sum(diff < 0)) / (len(x) * len(y))

def cliffs_delta_ci(x, y, n_boot=2000, alpha=0.05, rng=None):
    """Percentile bootstrap confidence interval for Cliff's delta."""
    rng = rng if rng is not None else np.random.default_rng(SEED)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        xb = rng.choice(x, size=len(x), replace=True)
        yb = rng.choice(y, size=len(y), replace=True)
        boots[b] = cliffs_delta(xb, yb)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return lo, hi

def cliffs_magnitude(d):
    a = abs(d)
    if a < 0.147: return 'negligible'
    if a < 0.330: return 'small'
    if a < 0.474: return 'medium'
    return 'large'

def fdr_bh(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out

def load_group_labels(path, id_col, group_col):
    """Read the clinical file and return a Series (index = TCGA patient barcode, values = group)."""
    ext = _Path(path).suffix.lower()
    if ext == '.csv':
        clin = pd.read_csv(path)
    elif ext == '.xls':
        clin = pd.read_excel(path, engine='xlrd')
    else:
        clin = pd.read_excel(path)
    for c in (id_col, group_col):
        if c not in clin.columns:
            raise KeyError(f"Column '{c}' not found. Available columns: {list(clin.columns)}")
    clin = clin[[id_col, group_col]].copy()
    clin['ModelName'] = clin[id_col].apply(extract_tcga_barcode)
    clin = clin.dropna(subset=['ModelName', group_col]).drop_duplicates('ModelName')
    return clin.set_index('ModelName')[group_col]

print("\n⏳ Non-parametric group comparison (Mann-Whitney U + Cliff's delta)...")

if GROUP_FILE is None or GROUP_COL is None:
    print("⚠️  Skipped: set GROUP_FILE and GROUP_COL in the configuration block to enable this section.")
else:
    try:
        group_labels = load_group_labels(GROUP_FILE, GROUP_ID_COL, GROUP_COL)
    except Exception as e:
        group_labels = None
        print(f"❌ Could not load group labels: {e}")

    if group_labels is not None:
        groups = df['ModelName'].map(group_labels)
        levels = sorted(groups.dropna().unique(), key=str)
        print(f"   Grouping variable '{GROUP_COL}': levels = {levels} "
              f"({groups.notna().sum()}/{len(df)} lesions with a label)")

        if len(dup_patients) > 0:
            print("   ⚠️  Some patients have several lesions → observations are not fully independent.")

        if len(levels) != 2:
            print(f"❌ Exactly 2 groups are required (found {len(levels)}). "
                  f"Binarize '{GROUP_COL}' (e.g., ER+ vs ER−) or subset two levels.")
        else:
            g1, g2 = levels
            mask1 = (groups == g1).values
            mask2 = (groups == g2).values
            print(f"   Group 1 = '{g1}' (n={mask1.sum()}) | Group 2 = '{g2}' (n={mask2.sum()})")
            print(f"   Sign convention: δ > 0 → '{g1}' tends to have HIGHER values than '{g2}'.")
            if min(mask1.sum(), mask2.sum()) < 10:
                print("   ⚠️  A group has fewer than 10 lesions: interpret results with caution.")

            rng = np.random.default_rng(SEED)
            rows = []
            for col in feature_cols:
                x = df.loc[mask1, col].dropna().values
                y = df.loc[mask2, col].dropna().values
                if len(x) < 2 or len(y) < 2:
                    continue
                u_stat, p_val = stats.mannwhitneyu(x, y, alternative='two-sided')
                delta = cliffs_delta(x, y)
                ci_lo, ci_hi = cliffs_delta_ci(x, y, n_boot=N_BOOT, alpha=ALPHA, rng=rng)
                letter, code = get_feature_code(col)
                rows.append({
                    'Feature':          col,
                    'Category':         CATEGORY_NAMES.get(letter, 'Other'),
                    'Letter':           letter,
                    f'n_{g1}':          len(x),
                    f'n_{g2}':          len(y),
                    f'Median_{g1}':     np.median(x),
                    f'Median_{g2}':     np.median(y),
                    'U_statistic':      u_stat,
                    'MannWhitney_p':    p_val,
                    'Cliffs_delta':     delta,
                    'CI95_low':         ci_lo,
                    'CI95_high':        ci_hi,
                    'Magnitude':        cliffs_magnitude(delta),
                    'AUC_equivalent':   (delta + 1) / 2,   # P(X>Y) + 0.5·P(X=Y)
                })

            df_group = pd.DataFrame(rows)
            df_group['MannWhitney_p_FDR'] = fdr_bh(df_group['MannWhitney_p'].values)
            df_group['Significant_FDR05'] = np.where(df_group['MannWhitney_p_FDR'] < ALPHA, 'Yes', 'No')
            df_group = df_group.sort_values('Cliffs_delta', key=lambda s: s.abs(), ascending=False)
            df_group.to_csv(outpath('group_comparison_MannWhitney_CliffsDelta.csv'), index=False)
            print("✅ 'group_comparison_MannWhitney_CliffsDelta.csv' generated.")

            n_sig_raw = (df_group['MannWhitney_p'] < ALPHA).sum()
            n_sig_fdr = (df_group['MannWhitney_p_FDR'] < ALPHA).sum()
            print(f"\n📊 {n_sig_raw}/{len(df_group)} features with p<{ALPHA} (uncorrected); "
                  f"{n_sig_fdr}/{len(df_group)} after Benjamini-Hochberg FDR.")
            print("   Top 10 features by |Cliff's delta|:")
            print(df_group.head(10)[['Feature', 'Cliffs_delta', 'Magnitude',
                                     'MannWhitney_p', 'MannWhitney_p_FDR']].to_string(index=False))

            # ── Forest-style plot of Cliff's delta with bootstrap 95% CI ─────
            plot_df = df_group.iloc[::-1].reset_index(drop=True)   # largest |δ| at the top
            letters = list(feature_groups.keys())
            palette = sns.color_palette('Set2', n_colors=max(len(letters), 3))
            colors  = [palette[letters.index(l) % len(palette)] for l in plot_df['Letter']]

            fig, ax = plt.subplots(figsize=(9, max(5, len(plot_df) * 0.28 + 1.5)))
            y_pos = np.arange(len(plot_df))
            xerr = np.vstack([
                np.clip(plot_df['Cliffs_delta'] - plot_df['CI95_low'], 0, None),
                np.clip(plot_df['CI95_high'] - plot_df['Cliffs_delta'], 0, None),
            ])
            ax.barh(y_pos, plot_df['Cliffs_delta'], color=colors, edgecolor='white',
                    xerr=xerr, error_kw=dict(ecolor='#333333', lw=0.8, capsize=2))
            ax.set_yticks(y_pos)
            ax.set_yticklabels([re.sub(r'\s*\(([A-Z]\d+)\)\s*$', r' [\1]', f) +
                                (' *' if s == 'Yes' else '')
                                for f, s in zip(plot_df['Feature'], plot_df['Significant_FDR05'])],
                               fontsize=FONT_SIZE - 2)
            ax.axvline(0, color='black', lw=0.8)
            for t in (0.147, 0.330, 0.474):
                ax.axvline(t,  color='grey', ls='--', lw=0.6, alpha=0.6)
                ax.axvline(-t, color='grey', ls='--', lw=0.6, alpha=0.6)
            ax.set_xlim(-1, 1)
            ax.set_xlabel(f"Cliff's delta  (δ>0: '{g1}' > '{g2}')")
            ax.set_title(f"Cliff's delta by feature — {GROUP_COL}\n"
                         f"(* = Mann-Whitney U, FDR-adjusted p<{ALPHA}; dashed lines = 0.147 / 0.33 / 0.474)",
                         fontsize=FONT_TITLE - 2, fontweight='bold')
            handles = [plt.Rectangle((0, 0), 1, 1, color=palette[i % len(palette)])
                       for i in range(len(letters))]
            ax.legend(handles, [CATEGORY_NAMES.get(l, l) for l in letters],
                      fontsize=FONT_SIZE - 3, loc='lower right', frameon=True)
            plt.tight_layout()
            save_fig(fig, 'CliffsDelta_by_feature.png')

# ============================================================
# 1️⃣1️⃣  FINAL SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("  EXPLORATORY ANALYSIS COMPLETED")
print("=" * 60)
print(f"  Patients/lesions   : {df.shape[0]}")
print(f"  Features           : {len(feature_cols)} ({len(feature_groups)} categories: {','.join(feature_groups.keys())})")
print(f"  Missing values     : {df[feature_cols].isna().sum().sum()}")
print(f"  Output folder      : {OUTPUT_DIR}")
print(f"\n  CSV files:")
print(f"    radiomics_with_ModelName.csv                    (dataset with normalized patient ID)")
print(f"    descriptive_statistics.csv                      (mean, SD, skew, kurtosis, normality per feature)")
print(f"    outlier_detection.csv                           (IQR-rule outliers, with IDs)")
print(f"    high_correlation_pairs_r090.csv                 (redundant pairs, |r|≥0.90)")
print(f"    pca_explained_variance.csv                      (variance explained per component)")
print(f"    group_comparison_MannWhitney_CliffsDelta.csv    (only if GROUP_FILE / GROUP_COL are set)")
print(f"\n  Figures (.png, 300 DPI):")
print(f"    Distribution_Category_{{K,E,T,G,M,S}}.png")
print(f"    Boxplot_Zscore_Category_{{K,E,T,G,M,S}}.png")
print(f"    Correlation_Matrix_Spearman.png")
print(f"    PCA_scree_and_projection.png")
print(f"    CliffsDelta_by_feature.png                      (only if GROUP_FILE / GROUP_COL are set)")
print("=" * 60)