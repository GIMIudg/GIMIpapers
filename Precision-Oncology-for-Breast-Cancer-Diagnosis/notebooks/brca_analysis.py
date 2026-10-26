"""
brca_analysis.py - helpers for BRCA_Metabolic_Flux_Pipeline_02.ipynb

Self-contained (numpy / scipy / scikit-learn / matplotlib only) so that every number in the
manuscript can be reproduced without optional packages:
  * figure style and the colour palette shared by all figures
  * UMAP (McInnes et al. 2018) implemented in numpy  -> `umap_embed`
  * statistics: Benjamini-Hochberg, Cliff's delta, Cramer's V, Kaplan-Meier, log-rank, Cox PH
  * MWUSelector: Mann-Whitney U + BH feature selection usable inside scikit-learn pipelines
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import curve_fit, minimize
from scipy.stats import mannwhitneyu, chi2_contingency, chi2
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import SpectralEmbedding
import matplotlib as mpl
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ palette (validated, see README)
COL = {
    # tissue
    "Normal": "#1baf7a", "Tumor": "#eb6834", "Metastatic": "#6250d6",
    # PAM50 (Normal-like shares the color of normal tissue)
    "LumA": "#2a78d6", "LumB": "#eda100", "Her2": "#e87ba4", "Basal": "#e34948", "Normal-like": "#1baf7a",
    # neutral
    "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#9e9d98", "grid": "#e6e5e0",
}
ALGO_COLORS = {  # fixed order, never cycled
    "Logistic Regression": "#2a78d6", "SVM": "#eb6834", "KNN": "#1baf7a", "Random Forest": "#eda100",
    "Gradient Boosting": "#e87ba4", "MLP": "#008300", "Decision Tree": "#6250d6", "Naive Bayes": "#e34948",
    "Dummy (prior)": "#9e9d98",
}
CLUSTER_COLORS = ["#6250d6", "#eda100", "#1baf7a", "#e87ba4", "#2a78d6", "#008300"]
SEQ = "Blues"          # magnitude
DIV = "RdBu_r"         # polarity (blue = lower, red = higher)


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": COL["ink2"], "axes.linewidth": 0.6, "xtick.color": COL["ink2"], "ytick.color": COL["ink2"],
        "axes.labelcolor": COL["ink"], "text.color": COL["ink"], "axes.grid": False,
        "figure.dpi": 110, "savefig.dpi": 300, "lines.linewidth": 1.5, "pdf.fonttype": 42,
    })


def panel_label(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def save(fig, path_noext):
    """PLOS: TIFF 300 dpi (LZW) + PNG preview."""
    fig.savefig(f"{path_noext}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{path_noext}.tif", bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    plos_tiff(f"{path_noext}.tif")


def plos_tiff(path, max_w=2250, max_h=2625):
    """PLOS ONE figure check: RGB (no alpha), 300 dpi, width <= 7.5 in, height <= 8.75 in, LZW."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    s = min(1.0, max_w / im.size[0], max_h / im.size[1])
    if s < 1:
        im = im.resize((int(im.size[0] * s), int(im.size[1] * s)), Image.LANCZOS)
    im.save(path, compression="tiff_lzw", dpi=(300, 300))


# ------------------------------------------------------------------ statistics
def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return out
    o = np.argsort(p[ok])
    r = p[ok][o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    t = np.empty(n); t[o] = r
    out[ok] = t
    return out


def cliffs_delta(x, y):
    u = mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1, u.pvalue


def cliffs_matrix(A, B):
    """vectorised Cliff's delta and MWU p for every column (A, B: samples x features)."""
    r = mannwhitneyu(A, B, axis=0, alternative="two-sided")
    return 2 * r.statistic / (A.shape[0] * B.shape[0]) - 1, r.pvalue


def cramers_v(a, b):
    t = pd.crosstab(a, b)
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan, np.nan
    c, p, _, _ = chi2_contingency(t)
    return np.sqrt(c / (t.values.sum() * (min(t.shape) - 1))), p


def kaplan_meier(t, e):
    t, e = np.asarray(t, float), np.asarray(e, int)
    times = np.unique(t[e == 1])
    s, S, lo, hi, var = 1.0, [], [], [], 0.0
    for u in times:
        n = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / n
        var += d / (n * (n - d)) if n > d else 0
        se = s * np.sqrt(var)
        S.append(s); lo.append(max(0, s - 1.96 * se)); hi.append(min(1, s + 1.96 * se))
    return np.r_[0, times], np.r_[1, S], np.r_[1, lo], np.r_[1, hi]


def logrank(t, e, g):
    """multi-group log-rank test; returns chi2, df, p."""
    t, e, g = np.asarray(t, float), np.asarray(e, int), np.asarray(g)
    groups = np.unique(g); k = len(groups)
    O = np.zeros(k); E = np.zeros(k); V = np.zeros((k, k))
    for u in np.unique(t[e == 1]):
        at = t >= u
        n = at.sum(); d = ((t == u) & (e == 1)).sum()
        nj = np.array([(at & (g == x)).sum() for x in groups])
        dj = np.array([((t == u) & (e == 1) & (g == x)).sum() for x in groups])
        O += dj; E += d * nj / n
        if n > 1:
            V += (d * (n - d) / (n - 1)) * (np.diag(nj / n) - np.outer(nj, nj) / n ** 2)
    z = (O - E)[:-1]; Vr = V[:-1, :-1]
    stat = float(z @ np.linalg.pinv(Vr) @ z)
    return stat, k - 1, float(chi2.sf(stat, k - 1))


def cox_ph(t, e, X, names):
    """Cox proportional hazards (Breslow ties), Newton via scipy; returns HR table."""
    t, e, X = np.asarray(t, float), np.asarray(e, int), np.asarray(X, float)
    o = np.argsort(-t); t, e, X = t[o], e[o], X[o]

    def nll(b):
        eta = X @ b
        ll = 0.0; g = np.zeros_like(b)
        m = eta.max()
        cs = np.cumsum(np.exp(eta - m)); csx = np.cumsum(np.exp(eta - m)[:, None] * X, axis=0)
        # risk set of i = all j with t_j >= t_i -> positions up to the last tie of t_i
        last = np.searchsorted(-t, -t, side="right") - 1
        ev = e == 1
        ll = (eta[ev] - (np.log(cs[last[ev]]) + m)).sum()
        g = (X[ev] - csx[last[ev]] / cs[last[ev]][:, None]).sum(0)
        return -ll, -g

    r = minimize(nll, np.zeros(X.shape[1]), jac=True, method="BFGS")
    b = r.x
    # observed information by finite differences of the gradient
    h = 1e-5; H = np.zeros((len(b), len(b)))
    for j in range(len(b)):
        d = np.zeros_like(b); d[j] = h
        H[:, j] = (nll(b + d)[1] - nll(b - d)[1]) / (2 * h)
    se = np.sqrt(np.diag(np.linalg.pinv(H)))
    from scipy.stats import norm
    zv = b / se
    return pd.DataFrame({"covariate": names, "HR": np.exp(b), "HR_low95": np.exp(b - 1.96 * se),
                         "HR_high95": np.exp(b + 1.96 * se), "p": 2 * norm.sf(np.abs(zv))})


# ------------------------------------------------------------------ feature selection
class MWUSelector(BaseEstimator, TransformerMixin):
    """Keep features with BH-adjusted Mann-Whitney p <= alpha (fitted on training data only).
    If fewer than k_min pass, the k_min smallest p-values are kept."""

    def __init__(self, alpha=0.05, k_min=10):
        self.alpha, self.k_min = alpha, k_min

    def fit(self, X, y):
        X, y = np.asarray(X, float), np.asarray(y)
        p = mannwhitneyu(X[y == 1], X[y == 0], axis=0, alternative="two-sided").pvalue
        p = np.where(np.isnan(p), 1.0, p)
        keep = bh(p) <= self.alpha
        if keep.sum() < self.k_min:
            keep = np.zeros_like(keep); keep[np.argsort(p)[: self.k_min]] = True
        self.support_, self.pvalues_ = keep, p
        return self

    def transform(self, X):
        return np.asarray(X, float)[:, self.support_]


# ------------------------------------------------------------------ UMAP (numpy implementation)
def _find_ab(spread=1.0, min_dist=0.1):
    def f(x, a, b):
        return 1.0 / (1.0 + a * x ** (2 * b))
    xv = np.linspace(0, spread * 3, 300)
    yv = np.where(xv < min_dist, 1.0, np.exp(-(xv - min_dist) / spread))
    (a, b), _ = curve_fit(f, xv, yv)
    return a, b


def _fuzzy_graph(X, k, metric):
    nn = NearestNeighbors(n_neighbors=k, metric=metric).fit(X)
    d, idx = nn.kneighbors(X)
    n = X.shape[0]
    rho = np.array([r[r > 0][0] if (r > 0).any() else 0.0 for r in d])
    target = np.log2(k)
    sig = np.ones(n)
    for i in range(n):                                   # binary search per point
        lo, hi, mid = 0.0, np.inf, 1.0
        di = np.maximum(d[i, 1:] - rho[i], 0)
        for _ in range(64):
            s = np.exp(-di / mid).sum()
            if abs(s - target) < 1e-5:
                break
            if s > target:
                hi = mid; mid = (lo + hi) / 2
            else:
                lo = mid; mid = mid * 2 if hi == np.inf else (lo + hi) / 2
        sig[i] = max(mid, 1e-3 * np.mean(d))
    w = np.exp(-np.maximum(d[:, 1:] - rho[:, None], 0) / sig[:, None])
    rows = np.repeat(np.arange(n), k - 1)
    P = sp.coo_matrix((w.ravel(), (rows, idx[:, 1:].ravel())), shape=(n, n)).tocsr()
    P = P + P.T - P.multiply(P.T)
    return P.tocoo()


def umap_embed(X, n_neighbors=15, min_dist=0.1, n_components=2, metric="euclidean",
               n_epochs=500, negative_sample_rate=5, random_state=42):
    """UMAP embedding (same algorithm and defaults as umap-learn 0.5: fuzzy simplicial set,
    spectral initialisation, SGD with negative sampling)."""
    rng = np.random.default_rng(random_state)
    X = np.asarray(X, float)
    P = _fuzzy_graph(X, n_neighbors, metric)
    keep = P.data >= P.data.max() / n_epochs
    head, tail, w = P.row[keep], P.col[keep], P.data[keep]
    a, b = _find_ab(1.0, min_dist)
    try:
        Y = SpectralEmbedding(n_components=n_components, affinity="precomputed",
                              random_state=random_state).fit_transform(P.tocsr())
        Y = 10 * (Y - Y.min(0)) / (Y.max(0) - Y.min(0) + 1e-12)
    except Exception:
        Y = rng.uniform(-10, 10, (X.shape[0], n_components))
    Y = Y + rng.normal(0, 1e-4, Y.shape)
    eps = w.max() / w                                    # epochs per sample
    nxt = eps.copy()
    n = Y.shape[0]
    for ep in range(n_epochs):
        alpha = 1.0 - ep / n_epochs
        m = nxt <= ep + 1
        if not m.any():
            continue
        h, t = head[m], tail[m]
        diff = Y[h] - Y[t]
        d2 = (diff ** 2).sum(1)
        coef = np.where(d2 > 0, -2 * a * b * d2 ** (b - 1) / (1 + a * d2 ** b), 0.0)
        g = np.clip(coef[:, None] * diff, -4, 4) * alpha
        np.add.at(Y, h, g)
        np.add.at(Y, t, -g)
        nxt[m] += eps[m]
        hn = np.repeat(h, negative_sample_rate)
        tn = rng.integers(0, n, hn.size)
        diff = Y[hn] - Y[tn]
        d2 = (diff ** 2).sum(1)
        coef = np.where(d2 > 0, 2 * b / ((0.001 + d2) * (1 + a * d2 ** b)), 0.0)
        g = np.where(d2[:, None] > 0, np.clip(coef[:, None] * diff, -4, 4), 4.0) * alpha
        g[tn == hn] = 0
        np.add.at(Y, hn, g)
    return Y


# ------------------------------------------------------------------ Fig 1 (workflow schematic)
def workflow_figure(path_noext):
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    set_style()
    plt.rcParams["savefig.dpi"] = 300
    fig, ax = plt.subplots(figsize=(7.5, 3.1))
    ax.set_xlim(0, 100); ax.set_ylim(0, 42); ax.axis("off")

    DATA, MODEL, ML, CLIN = "#eef3fb", "#f3f1fb", "#eaf7f1", "#fdf1ea"
    EDGE = {DATA: "#2a78d6", MODEL: "#6250d6", ML: "#1baf7a", CLIN: "#eb6834"}

    def box(x, y, w, h, title, lines, fill):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fill, ec=EDGE[fill], lw=0.9))
        ax.text(x + w / 2, y + h - 1.8, title, ha="center", va="top", fontsize=7.4, fontweight="bold", color=COL["ink"],
                linespacing=1.2)
        ax.text(x + w / 2, y + h - 3.0 - 2.9 * (title.count("\n") + 1), "\n".join(lines), ha="center", va="top", fontsize=6.3, color=COL["ink2"], linespacing=1.35)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=0.9, color=COL["ink2"],
                                     shrinkA=0, shrinkB=0))

    box(0.5, 23, 15, 18, "a) TCGA-BRCA\ndata", ["RNA-seq (FPKM)", "clinical data", "1,226 samples", "1,095 patients"], DATA)
    box(0.5, 1, 15, 18, "b) Manual\ncuration", ["breast cancer genes,", "reactions and", "metabolites"], DATA)
    box(20, 8, 16, 26, "c) Patient-specific\nGEMs", ["Recon3D", "XomicsToModel", "(thermoKernel)", "active gene:", "FPKM ≥ e²"], MODEL)
    box(40.5, 8, 18, 26, "d) Flux\nestimation", ["uniform medium", "(HPLM + lipids)", "FBA · pFBA · L1w", "multi-objective", "derived metrics", "model QC"], MODEL)
    box(63, 23, 16, 18, "e) Classification", ["tumor vs normal", "8 algorithms", "nested, patient-", "grouped CV"], ML)
    box(63, 1, 16, 18, "f) Clustering", ["PCA / UMAP", "7 algorithm families", "artifact checks", "clinical vs metabolic"], ML)
    box(83.5, 8, 16, 26, "g) Clinical\ninterpretation", ["receptor status", "PAM50 subtype", "clinical variables", "predicted from fluxes", "survival"], CLIN)

    arrow(16.3, 32, 19.6, 25); arrow(16.3, 10, 19.6, 17)
    arrow(36.8, 21, 40.1, 21)
    arrow(59.3, 25, 62.6, 31); arrow(59.3, 17, 62.6, 11)
    arrow(79.8, 31, 83.1, 25); arrow(79.8, 11, 83.1, 17)
    save(fig, path_noext)
