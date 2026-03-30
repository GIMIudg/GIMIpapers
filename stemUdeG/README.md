# 📊 STEM Admission Trends at Universidad de Guadalajara

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status: Completed](https://img.shields.io/badge/Status-Completed-blue)](https://github.com/Gpreciat/dataTRICKS)
[![Made with Python](https://img.shields.io/badge/Made%20with-Python%203.10-blue.svg)](https://www.python.org/)
[![Paper](https://img.shields.io/badge/View-Paper-9cf?logo=Overleaf)](https://example.com/paper-pdf) <!-- Reemplaza si tienes un link final -->

---

> 🧠 This repository contains the full codebase, data, and supplementary materials for the project **"AI-Based Analysis of STEM Admission Trends in the University of Guadalajara"**, presented at the The Mexican International Conference on Artificial Intelligence (MICAI) 2025.

---

## 🔍 Project Overview

This study analyzes longitudinal admission data (2003A–2025B) from the Universidad de Guadalajara (UdeG), focusing on STEM-related careers—especially engineering programs at the CUCEI campus.

We applied **machine learning** and **statistical modeling** techniques to:

- Detect long-term stagnation or decline in engineering interest
- Evaluate the impact of COVID-19 on student behavior
- Cluster programs by competitiveness using PCA
- Identify anomalous trajectories via Isolation Forest and LOF
- Forecast future applicant volume using SARIMA

All analyses are reproducible and open source.

---

## 📁 Repository Structure

```plaintext
projects/stemUdeG/
├── data/
│   └── PM TOTALES CU.xlsx                # Preprocessed dataset
├── notebooks/
│   ├── stemAdmissionPredictor.ipynb      # Formal analysis
└── README.md
