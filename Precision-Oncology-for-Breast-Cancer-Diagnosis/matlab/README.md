# MATLAB

Model generation and flux calculation (COBRA Toolbox v3, Gurobi). Scripts to be added:

- `modelGen.m` — builds one model per sample with XomicsToModel from `results/processedData/xomics_files` and `bibliomicData_new.xlsx` (`param.transcriptomicThreshold = 0`, `param.biomassRxnAbbr = 'biomass_ob'`). Output: `matlabResults/models`.
- `BRCA_fluxPipeline.m` — validates every model (PASS/WARN/FAIL), computes FBA, pFBA, L1w, L2, L2w and the metrics, and writes `FeatureMatrix_TumorPhenotype_All.csv`, `ModelQC.csv`, `FeatureDictionary.csv`, `ParetoSurface_100sol.csv` to `matlabResults/`, which notebook 02 reads.
