# GEMs Code for Construction

This folder contains the complete computational pipeline for constructing and analyzing personalized Genome-Scale Metabolic Models (GEMs) from transcriptomic data.

## 📂 Directory Structure

### `GEMs_construction/`
This folder contains the scripts needed to process the raw RNA-seq data and construct the metabolic models:
- **`creator.py`**: A Python script that maps Ensembl IDs to Entrez IDs, processes patient transcriptomic profiles, handles duplicated genes by averaging their expression, preserves unexpressed genes (zero values), and exports the data into TXT (`Xomics`) files required for model building.
- **`utils.py`**: Helper functions and utilities for `creator.py`.
- **`GEMs_creation.m`**: A MATLAB script that takes the generated `Xomics` files and uses thermodynamic and metabolic network constraints (e.g., tINIT algorithm) to construct the final personalized genome-scale metabolic models.

### `GEMS_Exploratory_Analysis/`
This folder contains scripts for the quality control, metrics extraction, and initial statistical exploration of the constructed metabolic models:
- **`GEMs_metrics.m`**: MATLAB script to evaluate the topological and functional metrics of the constructed GEMs (e.g., number of reactions, metabolites, constraints).
- **`GEMs_exploratory_plots_1226.R`**: R script used to generate statistical plots and visualizations (boxplots, distributions) of the extracted model metrics.
- **`model_metrics_1226.ipynb`**: A Jupyter Notebook for interactive data exploration and analysis of the metabolic properties of the 1,226 patient models.

## 🚀 Workflow Summary
1. **Preprocessing (`creator.py`)**: Convert raw RNA-seq data (FPKM) into sample-specific `Xomics` expression files, taking care to map IDs properly and resolve duplicates.
2. **Model Building (`GEMs_creation.m`)**: Run the MATLAB pipeline using the `Xomics` files to generate the context-specific thermodynamic GEMs.
3. **Exploration & Quality Control (`GEMS_Exploratory_Analysis/`)**: Extract metrics from the resulting models and generate exploratory plots in R and Python to validate their consistency.
