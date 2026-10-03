# BRCA Metabolic Flux Analysis — Unified Pipeline

Reproducible Google Colab notebook for the paper:

> **Patient-Specific Metabolic Fluxes Reveal Functional Organization and Heterogeneity in Breast Cancer**  
> Ruiz Robles E., Rincón-Ballesteros R., Chacón Méndez S.A., Alvarez-Padilla F.J., Preciat G.  
> University of Guadalajara, Mexico. 2025.

---

## Project Structure

```
Precision-Oncology-for-Breast-Cancer-Diagnosis/
├── Unsupervised_learning_scripts_and_data/
│   ├── Clinical_data_and_models_ids/      ← Clinical metadata and GEM files
│   │   ├── All_models_created/            ← Contains all 1,226 generated models
│   │   ├── Clinical_Data/                 ← Downloaded clinical and survival datasets. Metadata was downloaded from [UCSC Xena](https://xena.ucsc.edu/), while the clinical and survival data were obtained from [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA).
│   │   ├── GEMs_Data_for_construction/    ← Data used during GEMs construction
│   │   └── Metabolic_Data/                ← Flux results and Pareto solutions
│   └── src/                               ← Analysis scripts (MATLAB, Python, R)
├── main_notebook/
│   ├── BRCA_Paired_Patient_Classification.ipynb ← Supervised learning notebook
│   └── Figures/                           ← Output figures and plots
├── requirements.txt                       ← Python dependencies
└── README.md                              ← Project overview & reproducibility guide
```

- **`Clinical_Data/`**: Metadata was downloaded from [UCSC Xena](https://xena.ucsc.edu/), while the clinical and survival data were obtained from [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA).

## How to Run

### Environment Setup (Python)

It is highly recommended to use a virtual environment to avoid dependency conflicts:

```bash 
# Create a virtual environment
python3 -m venv env
source env/bin/activate  # On Windows use: env\Scripts\activate

# Install the required dependencies
pip install -r requirements.txt

# Launch Jupyter Notebook
jupyter notebook
```
Then, you can open `main_notebook/BRCA_Paired_Patient_Classification.ipynb` to execute the supervised learning pipeline. For the rest of the analyses, follow the **Detailed Reproducibility Guide** below.

---

## Detailed Reproducibility Guide

Follow the steps below to fully reproduce the results of this project across the various scripts. It is highly recommended to install the dependencies listed in `requirements.txt` to ensure the Python algorithms run exactly as intended.

### 1. Creation of Metabolic Models

To construct the metabolic models, we utilized "Recon3D_301.mat" as the base model, which contains all known human metabolic reactions. Using the `xomicsToModel` algorithm implemented in the COBRAToolbox, 1,226 patient-specific breast cancer models were created.

#### Data Processing
1. **Gene Expression Data:** Gene expression profiles for 1,226 breast cancer patients were obtained from the TCGA (The Cancer Genome Atlas) project.
2. **Structuring:** The data was structured into an individual `.txt` file per patient.
3. **Preprocessing:** During preprocessing, we found that multiple Ensembl identifiers mapped to the same Entrez ID. To resolve this and prevent gene duplication, we calculated the average of their corresponding expression levels. Genes with null (zero) expression values were preserved.
4. **Model Initialization:** The resulting matrix was processed using the `xomicsToModel` function. The scripts responsible for this data processing are `creator.py` and `utils.py`.
5. **Bibliomic Data:** Simultaneously, an Excel file named `Bibliomics_Data` was curated. It contains a list of genes associated with breast cancer, the mandatory inclusion of a biomass reaction, and specific constraints required to generate consistent personalized models.

#### Model Construction
- **Script:** A MATLAB script named `Construction_Gems` was used. The `xomicsToModel` tool utilizes the gene expression and bibliomics data to find all possible reactions given the selected genes, ensuring the models are thermodynamically consistent. This procedure is based on the methodology described in:
  > Preciat, G., Wegrzyn, A.B., Luo, X. et al. (2026). XomicsToModel: omics data integration and generation of thermodynamically consistent metabolic models. *Nat Protoc*, 21, 2665–2706. [https://doi.org/10.1038/s41596-025-01288-9](https://doi.org/10.1038/s41596-025-01288-9)
- **Data directory:** `Unsupervised_learning_scripts_and_data/Clinical_data_and_models_ids/GEMs_Data_for_construction`
- **Code directory:** `Unsupervised_learning_scripts_and_data/src/GEMS_Code_for_Construction/GEMs_construction`

### 2. Obtaining Metabolic Fluxes

Once the models were built, various general and secondary metabolic fluxes were calculated. These fluxes act as descriptors that reflect the heterogeneity across the different metabolic models.

- **Data directory:** `Unsupervised_learning_scripts_and_data/Clinical_data_and_models_ids/GEMs_Data_for_construction`
- **Code directory:** `Unsupervised_learning_scripts_and_data/src/Metabolic_fluxes_calculation_MATLAB`
- **Main Script:** `Fluxes_calculation.m` runs the principal flux optimizations.
- **Pareto Optimization Script:** `Pareto_multiobjective_1226_models.m` calculates multi-objective optimality in the metabolic models. This Pareto procedure follows the same methodology described in: 
  > Dai, Z., Yang, S., Xu, L., et al. (2019). Identification of cancer-associated metabolic vulnerabilities by modeling multi-objective optimality in metabolism. *Cell Communication and Signaling*, 17, 124. [https://doi.org/10.1186/s12964-019-0439-y](https://doi.org/10.1186/s12964-019-0439-y)
- **Model Names Extraction:** `Models_names_saved.m` is a utility script that exclusively extracts and saves the names of the generated models.
- **Auxiliary functions:** `PrepareModel.m` and `countcarbons.m` ensure that the biomass reaction is correctly present in the model before running the flux optimizations.

### 3. Supervised Learning Pipeline

A pipeline was designed for processing metabolic data and incorporating it into supervised learning models to predict whether samples are tumor or normal tissue.

- **Directory:** `main_notebook/`
- **Dependencies:** Python algorithms require exact library versions to ensure reproducibility, as updates to libraries may alter random number generation or specific metric calculations. Please install them using the `requirements.txt` file located in the root of the repository.

### 4. Unsupervised Learning and Clinical/Metabolic Correlations

We then searched for correlations between clinical and metabolic data using unsupervised learning models and statistical metrics.

#### Clinical Data Analysis
Clinical data was processed by performing a refined selection of variables across different datasets. Two dimensionality reduction methods (PCA and UMAP) and various unsupervised learning models were applied.
- **Code directory:** `Unsupervised_learning_scripts_and_data/src/Clustering_and_data_analysis_PYTHON/Clinical_data_analysis/ML_models_using_clinical_data`
- This folder contains `.py` scripts for PCA and UMAP that run the models using the refined variable selection.
- **Results:** The outputs are saved in the corresponding `results/` folders.

#### Metabolic Data Analysis
The exact same dimensionality reduction and unsupervised learning procedures were applied to the metabolic data.
- **Code directory:** `Unsupervised_learning_scripts_and_data/src/Clustering_and_data_analysis_PYTHON/Metabolic/ML_models_using_metabolic_data`

#### Cluster Comparisons & Concordance
Once clusters were generated from both clinical and metabolic data, the algorithms were compared using the Adjusted Rand Index (ARI).
- **Code directory:** `Unsupervised_learning_scripts_and_data/src/Clustering_and_data_analysis_PYTHON/Cluster_correlations`
- **Script:** `correlations_main.py` is responsible for finding the pair of clinical and metabolic algorithms with the highest number of concordant patients. This identified the divergent group with a characteristic quiescent signature.
- **Generated Results:**
  - `results/divergent_patients_study_pyn.csv`: Contains the patient codes and their divergent classification.
  - `results/core_patients_correlation_ParetoAndNorms`: Contains the patients belonging to the "core" group.
  - Additional outputs include figures of the selected clusters, a heatmap of the algorithms that grouped similarly, and a correlation matrix of the groups.

#### Statistical Analysis
Finally, significant differences and effect sizes between groups were analyzed.
- **Metrics:** Mann-Whitney U test, Cliff's Delta, and Benjamini-Hochberg correction.
- **Script:** `Clustering_Correlation_Analysis.ipynb` (written in R due to its superior capabilities for statistical graphics and biological analysis).
- **Results:** All generated figures and plots from this analysis are stored in the `Figures_PLOS_v6` folder.

---

## What the Notebook Produces

### Supervised Learning 
| Figure | Description |
|---|---|
| `Fig04_UMAP_before_after_feature_selection.png` | UMAP projection before/after feature selection |
| `Fig05_confusion_matrices.png` | KNN vs Decision Tree confusion matrices |
| `Fig06_ROC_AUC_curves.png` | ROC-AUC curves for all 5 classifiers |
| `Fig07_decision_matrix.png` | Multi-criteria classifier comparison heatmap |

### Unsupervised Clustering 
| Figure | Description |
|---|---|
| `Fig09_clinical_cluster_3D.png` | Clinical patient clustering — 3D view |
| `Fig10_metabolic_cluster_3D.png` | Metabolic patient clustering — 3D view |
| `Fig08_top10_ARI_heatmap.png` | Top 10 concordances (ARI) clinical vs metabolic |
| `Fig02_contingency_matrix.png` | Patient distribution clinical vs metabolic clusters |
| `Fig03_cohort_pie.png` | Cohort composition: Core vs Divergent |

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
