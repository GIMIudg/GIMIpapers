# MATLAB

COBRA Toolbox v3 + Gurobi. Run from the folder that contains `bibliomicData.xlsx`, `matlabResults/` and `pythonResults/` (Drive: `MCBCI_analisisCancer/results`).

| Script | Input | Output |
|---|---|---|
| `modelGen.m` | `bibliomicData.xlsx`, `pythonResults/processedData/xomics_files/*_Xomics.txt` | `matlabResults/models/<sample>_Xomics/Model.mat` (XomicsToModel) |
| `BRCA_fluxAnalysis.m` | `matlabResults/models/` | `matlabResults/FeatureMatrix_TumorPhenotype_All.csv`, `ModelQC.csv`, `ParetoSurface_100sol.csv`, `FeatureDictionary.csv`, `MaxBiomass_perModel.csv` — everything `notebooks/BRCA_Metabolic_Flux_Pipeline_02.ipynb` reads |

`BRCA_fluxAnalysis.m` replaces `Fluxes_calculation.m`, `Norms_1_2_fluxes-nurvo.mlx`, `prepareModel.m` and `countCarbons.m` (now in `old/flux_calculation_matlab/2026-09_versions/`): same five flux solutions (FBA, pFBA, L1w, L2, L2w), same metrics and column names, Recon3D reaction IDs, model QC (PASS/WARN/FAIL) and resumable runs (`matlabResults/perModel/`). Its header lists every difference with the old scripts.
