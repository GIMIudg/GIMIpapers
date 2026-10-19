# Precision Oncology for Breast Cancer Diagnosis
Master’s Thesis Research Repository

## Overview

This directory contains the datasets and resources used for the construction of genome-scale metabolic models (GEMs), metabolic flux analysis, and clinical data integration for precision oncology research in breast cancer.

---

## `/Clinical_Data`

This directory contains the clinical datasets used for unsupervised learning analysis. Two datasets were obtained directly from the National Cancer Institute (TCGA-BRCA), while `MetaData.xlsx` was downloaded from the UCSC Xena Browser. The metadata file contains additional clinical and phenotypic variables associated with the same patients.

---

## `/GEMS_Data_for_Construction`

This directory contains the bibliomic and transcriptomic data required for the reconstruction of genome-scale metabolic models using the `XomicsToModel` pipeline.

It includes:
- Transcriptomic data from 1,226 patients
- Input data required for model reconstruction
- A list of all generated metabolic models

---

## `/Metabolic_Data`

This directory contains the final datasets used for unsupervised learning and metabolic analysis.

The files include:
- Metabolic flux distributions
- Pareto surface calculations
