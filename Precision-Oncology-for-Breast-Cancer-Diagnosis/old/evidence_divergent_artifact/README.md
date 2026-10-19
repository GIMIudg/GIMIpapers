# Evidence: the "Divergent" group of the submitted manuscript is an artefact

Kept on purpose; do not delete. Files restored from commit `e6676451b71789b1853f1e77e6dca941f97f4a72` (main, before the 2026-10 reorganisation), with identical content (same SHA-256 as the LFS objects).

| File | Original path (commit e667645) | Why it matters |
|---|---|---|
| `pacientes_divergentes_para_estudio_pyn.csv` | `old/src/Clustering and data analysis PYTHON/Cluster_correlations/RESULTADOS_MEJORADOS/` | The 40 "Divergent" tumours of the paper. 24 have no MetaData (42% of them from tissue source site EW, vs 2.7% in the core; p≈1e-9) and 16 are models that XomicsToModel had to relax (biomass 0). |
| `pacientes_core_correlacion_Paretoynormas.csv`, `Estadisticas_Metabolicas_CoreVsDivergent.csv`, `Resumen_Significancia_PorGrupo.csv` | same folder | Core list and Core vs Divergent statistics reported in the paper. |
| `correlation_main_original.py` | `…/Cluster_correlations/correlation_main.py` | Lines ~128–175: "Core" = largest cell of the clinical × metabolic contingency table; every other patient is labelled divergent. |
| `ML_clinical_data_UMAP_seleccion_reducida_2_MODIFICACIONVARIABLES.py` (and PCA) | `…/Clinical data analysis/ML models using clinical data/` | Lines ~325–332: missing clinical values imputed as 0 / "Missing", so patients without MetaData cluster together in the clinical space. |
| `FeatureMatrix_TumorPhenotype_All_previous.csv` | `old/Clinical_data_and_models_ids/Metabolic_Data/` | Feature matrix of the 2025 models (LFS pointer). |

The new pipeline guards against both causes: MATLAB marks relaxed/zero-biomass models as `QC_pass = 0`, and notebook 02 tests every cluster against missing MetaData, tissue source site and QC warnings before interpreting it.
