# old

Code and inputs of the previous analysis (2025 models, submitted manuscript). Not used by the current pipeline; kept for traceability.
`inventory_2026-10-05.csv` lists every file of the repository before the reorganisation: original path, new path, what it is, which script generated it, what used it and the decision (old / toBeDeleted).

| Folder | Content | Replaced by |
|---|---|---|
| `evidence_divergent_artifact/` | files documenting the "Divergent" artefact — **keep** | — |
| `gem_construction/` | `creator.py` (expression → genes,expVal), `GEMs_creation.m`, `inputs/` (bibliomicData, Recon3D_301.mat, model IDs) | notebook 01, `matlab/modelGen.m`, `bibliomicData_new.xlsx` |
| `flux_calculation_matlab/` | `Fluxes_calculation.m`, `prepareModel.m`, Pareto scripts | `matlab/BRCA_fluxPipeline.m` |
| `gem_exploratory_analysis/` | model size / similarity scripts (R, MATLAB, notebook) | notebook 02 §7 |
| `clinical_ml/` | clinical clustering scripts (`previous_versions/`: earlier variants) | notebook 02 §5 |
| `metabolic_ml/` | metabolic clustering scripts | notebook 02 §4 |
| `cluster_correlations/` | clinical–metabolic concordance (Core/Divergent) | notebook 02 §5–6 |
| `pareto_exploration/` | Pareto notebook | notebook 02 §7 |
| `paired_classification_notebooks/` | student notebooks of Oct 2026 (paired classification, script explanation) | notebook 02 §3 |
| `previous_readmes/` | previous READMEs, `.gitignore`, `.gitattributes`, `requirements.txt` | — |
