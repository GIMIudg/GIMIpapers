# BRCA Metabolic Flux Analysis — Unified Pipeline

Reproducible Google Colab notebook for the paper:

> **Patient-Specific Metabolic Fluxes Reveal Functional Organization and Heterogeneity in Breast Cancer**  
> Ruiz Robles E., Rincón-Ballesteros R., Chacón Méndez S.A., Alvarez-Padilla F.J., Preciat G.  
> University of Guadalajara, Mexico. 2025.

---

## Project Structure

```
BRCA_Metabolic_Flux_Analysis/
├── notebooks/
│   └── BRCA_Metabolic_Flux_Pipeline.ipynb   ← Main Google Colab notebook
├── data/
│   ├── MetaData.xlsx                          ← Molecular metadata (ER, PR, HER2, Subtype, Ancestry)
│   ├── Model's_ids.txt                        ← List of 1,226 patient-specific GEM IDs
│   └── TCGA-BRCA.survival.tsv.gz             ← Overall survival data (OS.time, OS)
└── results/                                   ← (Created at runtime in Colab)
```

### Large Files — Automatically Downloaded in Colab

| File | Source | Size |
|---|---|---|
| `FeatureMatrix_TumorPhenotype_All.csv` | GitHub LFS (GIMIpapers) | ~179 MB |
| `TCGA-BRCA.clinical.tsv` | TCGA GDC | ~5 MB |
| `ParetoSurface_*.csv` | GitHub LFS (GIMIpapers) | ~10–100 MB |

---

## How to Run

### Option A — Google Colab (recommended)

1. Upload this folder to Google Drive (or clone the repo)
2. Open `notebooks/BRCA_Metabolic_Flux_Pipeline.ipynb` in Google Colab
3. Set runtime: **Runtime → Change runtime type → T4 GPU** (optional, speeds up UMAP)
4. Run all cells: **Runtime → Run all**
5. Expected runtime: ~60–90 min

### Option B — Local execution

```bash
pip install umap-learn hdbscan lifelines openpyxl pingouin statsmodels
jupyter notebook notebooks/BRCA_Metabolic_Flux_Pipeline.ipynb
```

---

## What the Notebook Produces

### Supervised Learning (Sections 4)
| Figure | Description |
|---|---|
| `Fig04_UMAP_before_after_feature_selection.png` | UMAP projection before/after feature selection |
| `Fig05_confusion_matrices.png` | KNN vs Decision Tree confusion matrices |
| `Fig06_ROC_AUC_curves.png` | ROC-AUC curves for all 5 classifiers |
| `Fig07_decision_matrix.png` | Multi-criteria classifier comparison heatmap |

### Unsupervised Clustering (Sections 5–8)
| Figure | Description |
|---|---|
| `Fig09_clinical_cluster_3D.png` | Clinical patient clustering — 3D view |
| `Fig10_metabolic_cluster_3D.png` | Metabolic patient clustering — 3D view |
| `Fig08_top10_ARI_heatmap.png` | Top 10 concordances (ARI) clinical vs metabolic |
| `Fig02_contingency_matrix.png` | Patient distribution clinical vs metabolic clusters |
| `Fig03_cohort_pie.png` | Cohort composition: Core vs Divergent |

### Divergent Subgroup Analysis (Sections 9–13)
| Figure | Description |
|---|---|
| `Fig11_CliffsDelta_forest_plot.png` | Top 40 metabolic features — Cliff's delta effect sizes |
| `Fig12_KaplanMeier_survival.png` | Kaplan–Meier overall survival curves |
| `Fig13_heatmap_metabolic_signatures_subtypes.png` | Metabolic signatures across molecular subtypes |
| `Fig03_metabolic_pathways_bubble.png` | Significantly altered metabolic pathways |
| `Fig_Pareto_panel.png` | Pareto front multi-objective trade-offs |

---

## Key Results Reproduced

| Result | Paper | Notebook |
|---|---|---|
| Best classifier accuracy | KNN 0.988 | Computed |
| Best ROC-AUC | KNN/SVM 1.000 | Computed |
| Metabolic clustering Silhouette | ~0.98 | Computed |
| Clinical clustering Silhouette | ~0.87 | Computed |
| Divergent group size | ~3.3% (n≈39) | Computed |
| Divergent subgroup: TNBC enrichment | χ², p=1.99×10⁻¹² | Computed |

---

## Data Availability

- **Repository:** https://github.com/GIMIudg/GIMIpapers/tree/main/Precision-Oncology-for-Breast-Cancer-Diagnosis
- **Zenodo archive:** https://doi.org/10.5281/zenodo.19339596
- **TCGA-BRCA:** https://portal.gdc.cancer.gov/projects/TCGA-BRCA

---

## Citation

If you use this pipeline, please cite:

```
Ruiz Robles E., Rincón-Ballesteros R., Chacón Méndez S.A., Alvarez-Padilla F.J., Preciat G. (2025).
Patient-specific metabolic fluxes reveal functional organization and heterogeneity in breast cancer.
DOI: 10.5281/zenodo.19339596
```
