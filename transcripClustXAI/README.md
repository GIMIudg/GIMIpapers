# 🧬 T2D Metabolic Subgroups Discovery via Liver Transcriptomics

 🧠 This repository contains methodology and image results from the study with the objective of identifying **latent metabolic subpopulations in Type 2 Diabetes (T2D)** by combining unsupervised clustering on liver transcriptomic data with post-hoc clinical interpretability.

-----

## 🔍 Overview

Type 2 Diabetes is a highly heterogeneous disorder where traditional monolithic classifications often fail to capture patient-specific variability. This project leverages a dataset of 910 subjects at risk of T2D to uncover hidden metabolic structures.

We applied **machine learning** and **bioinformatics** techniques to:

  - Map and filter 70,523 microarray probes into biologically functional metabolic families.
  - Cluster binary gene activity profiles to find patient subgroups.
  - Extract interpretable clinical rules governing these groups using decision trees.
  - Bridge the gap between high-dimensional molecular signatures and observable physiological phenotypes.

-----

## ⚙️ Methodology

```plaintext
Raw Microarray Data → ID Mapping → KEGG Family Stratification → Binarization & PCA → DBSCAN Clustering → Decision Tree Explainer
```

1.  **Schema Extraction & Processing** Gene expression data (GSE130991) is mapped from HTA 2.0 probes to Entrez IDs, filtering out non-coding sequences. Features are then stratified into 10 functional families using the KEGG database (e.g., Lyases, Ligases, Translocases).
2.  **Feature Engineering** Continuous variables are $\log_2$-transformed and binarized (threshold = 1.5) to capture functional activation patterns, reducing experimental noise. PCA is applied to retain 85% variance.
3.  **Unsupervised Clustering** DBSCAN, K-Means, and BIRCH are evaluated. DBSCAN is selected for its superior ability to handle irregular cluster shapes and complex densities in biological data.
4.  **Supervised Interpretation** A post-hoc analysis using Decision Trees acts as an explainer, mapping the high-dimensional transcriptomic subgroups back to 17 standard clinical variables (like BMI, HbA1c, and Triglycerides).

-----

## 🚀 Key Results

  - **Dimensionality Mitigation:** Grouping genes into functional metabolic families acts as biological feature selection, successfully bypassing the curse of dimensionality ($p \gg n$) and improving computational stability.
  - **The Lipid-Driven Dichotomy:** In the *Lyases* family, subgroups are overwhelmingly driven by systemic **Triglycerides** (87.15%), showing a distinct threshold-like transcriptomic response.
  - **Multifactorial Reprogramming:** In the *Ligases* family, the stratification is shaped by a combination of cumulative stressors: **BMI** (33.95%), **HbA1c** (25.44%), and **LDL Cholesterol** (16.38%).
  - **Insulin Sensitivity:** *Translocases* exhibit highly fragmented, individualized groupings strongly associated with Fasting Insulin (42.43%).

-----

## 📁 Repository Structure

```plaintext
project/
├── data/
│   ├── raw/                  # GSE130991 SOFT files and clinical metadata
│   └── processed/            # Binary activity matrices by KEGG family
├── src/
│   ├── preprocess.py         # ID mapping and binarization pipeline
│   └── clustering.py         # K-Means, BIRCH, and DBSCAN implementations
├── notebooks/
│   ├── 01_feature_extraction.ipynb
│   ├── 02_clustering_evaluation.ipynb
│   └── 03_decision_tree_explainers.ipynb
└── README.md
```