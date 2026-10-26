# Patient-specific metabolic fluxes in breast cancer (TCGA-BRCA)

Code for the manuscript *Patient-Specific Metabolic Fluxes: Functional Organization and Heterogeneity in Breast Cancer* (PLOS ONE, PONE-D-26-33254).
Genome-scale metabolic models (Recon3D + XomicsToModel) are built for 1,226 TCGA-BRCA samples (1,106 primary tumours, 113 normal tissues, 7 metastases), fluxes are computed with several objective functions, and flux-derived features are analysed with supervised, unsupervised and clinical methods.

## Workflow

```
data/rawData (Drive)                          raw TCGA-BRCA downloads
   │
   ▼  notebooks/BRCA_Metabolic_Flux_Pipeline_01.ipynb   data processing
results/pythonResults/processedData           xomics_files/*_Xomics.txt, clinical_processed.csv, sample_inventory.csv
   │
   ▼  matlab/modelGen.m                        models (XomicsToModel)
results/matlabResults/models
   │
   ▼  matlab/BRCA_fluxAnalysis.m               fluxes, metrics, model QC, Pareto
results/matlabResults                         FeatureMatrix_TumorPhenotype_All.csv, ModelQC.csv, ParetoSurface_100sol.csv
   │
   ▼  notebooks/BRCA_Metabolic_Flux_Pipeline_02.ipynb   supervised, unsupervised and clinical analysis
results/pythonResults                         figures/ and tables/
```

Data and results are **not** stored in this repository; they live in the project Drive folder (`MCBCI_analisisCancer/`):

| Folder | Content |
|---|---|
| `data/rawData/` | `TCGA-BRCA_log_FPKM.tsv`, `TCGA-BRCA.clinical.tsv` (GDC), `TCGA-BRCA.survival.tsv.gz`, `MetaData.xlsx` (UCSC Xena) |
| `data/reference/` | `Recon3D_genes.xlsx` (optional QC in notebook 01) |
| `results/bibliomicData.xlsx` | bibliomic data for `modelGen.m` |
| `results/pythonResults/processedData/` | outputs of notebook 01 |
| `results/matlabResults/` | models and outputs of `BRCA_fluxAnalysis.m` |
| `results/pythonResults/` | outputs of notebook 02 |

## How to run

1. `pip install -r requirements.txt`
2. Set the Drive path: `export BRCA_DRIVE_DIR=/path/to/MCBCI_analisisCancer` (or edit `DRIVE_DIR` in the first cell).
3. Run `notebooks/BRCA_Metabolic_Flux_Pipeline_01.ipynb`. The first run queries MyGene.info (internet) and caches the Ensembl→Entrez map.
4. Build the models with `matlab/modelGen.m`, then run `matlab/BRCA_fluxAnalysis.m` (COBRA Toolbox v3, Gurobi). Both ask for the `results` folder.
5. Run `notebooks/BRCA_Metabolic_Flux_Pipeline_02.ipynb`.

## Repository layout

| Path | Content |
|---|---|
| `notebooks/BRCA_Metabolic_Flux_Pipeline_01.ipynb` | data processing, expression → `genes,expVal` |
| `notebooks/BRCA_Metabolic_Flux_Pipeline_02.ipynb` | supervised, unsupervised and clinical analysis; writes every number, table and figure of the manuscript |
| `notebooks/brca_analysis.py` | helpers for notebook 02: color palette, figure export (PLOS TIFF), statistics, UMAP, Kaplan–Meier/log-rank/Cox, workflow figure |
| `matlab/modelGen.m` | model generation (XomicsToModel) |
| `matlab/BRCA_fluxAnalysis.m` | fluxes, metrics, model QC and Pareto surfaces for notebook 02 |
| `old/` | previous code, kept for traceability (see `old/README.md`, `old/inventory_2026-10-05.csv`) |
| `toBeDeleted/` | copies, session files and generated results waiting for review before deletion (see `toBeDeleted/README.md`) |

## Authors

Ruiz Robles E., Rincón-Ballesteros R., Chacón Méndez S.A., Alvarez-Padilla F.J., Preciat G. — Universidad de Guadalajara.
