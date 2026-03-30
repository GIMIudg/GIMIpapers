
# 🧬 Metabolic Flux Feature Engineering (RBC-GEM)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Status: In Progress](https://img.shields.io/badge/Status-In%20Progress-orange)]()  
[![MATLAB](https://img.shields.io/badge/MATLAB-required-blue)](https://www.mathworks.com/)  
[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)  

---

> 🧠 Framework for extracting **interpretable features from metabolic flux spaces** using genome-scale metabolic models (GEMs) and MCMC sampling.

---

## 🔍 Overview

This project develops a pipeline to:

- Compare **healthy vs pathological metabolic states**
- Sample the flux space
- Define a **healthy reference (FBA)**
- Extract **biologically meaningful features** from flux vectors
- Analyze states using **interpretable machine learning (decision trees)**

---

## 🧬 Model

- Base model: RBC-GEM  
- System: Red blood cell (RBC)  
- Case study: **G6PD deficiency**

---

## ⚙️ Pipeline

```

GEM → perturbation → sampling → fluxes → feature engineering → analysis

```

---

## 🧩 Features

- **Intrinsic**: norms, sparsity, entropy  
- **Comparative**: distances to healthy state  
- **Biochemical**: ATP, NADPH, redox balance  
- **Pathway-level**: glycolysis, PPP usage  
- **Advanced**: Mahalanobis distance, PCA space  

---

## 📁 Structure

```

project/
├── modelsAndSamples/   # Models + sampled fluxes
├── matlab/             # Model modification + sampling
├── python/
│   ├── scripts/        # Feature engineering
│   └── notebooks/      # Analysis
└── README.md

```

---

## 🚀 Workflow

1. **MATLAB**
   - Modify model (simulate G6PD)
   - Run sampling

2. **Python**
   - Compute features
   - Analyze + visualize
   - Train decision trees

---

## 🎯 Goal

Define a **feature-based representation of metabolic states** to quantify proximity to health and enable interpretable analysis of perturbations.
