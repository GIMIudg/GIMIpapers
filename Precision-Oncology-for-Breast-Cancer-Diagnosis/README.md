# Patient-specific metabolic fluxes in breast cancer (TCGA-BRCA)

Code for the manuscript *Patient-Specific Metabolic Fluxes: Functional Organization and Heterogeneity in Breast Cancer* (PLOS ONE, PONE-D-26-33254).
Genome-scale metabolic models (Recon3D + XomicsToModel) are built for 1,226 TCGA-BRCA samples (1,106 primary tumours, 113 normal tissues, 7 metastases), fluxes are computed with several objective functions, and flux-derived features are analysed with supervised, unsupervised and clinical methods.

## Workflow

```
Data/rawData (Drive)                         raw TCGA-BRCA downloads
   │
   ▼  BRCA_Metabolic_Flux_Pipeline_01.ipynb   data processing
results/processedData (Drive)                xomics_files/*_Xomics.txt, clinical_processed.csv, sample_inventory.csv
   │
   ▼  matlab/ (modelGen.m, BRCA_fluxPipeline.m)   models (XomicsToModel) and fluxes
matlabResults (Drive)                        models, FeatureMatrix_TumorPhenotype_All.csv, ModelQC.csv, ParetoSurface_100sol.csv
   │
   ▼  BRCA_Metabolic_Flux_Pipeline_02.ipynb   supervised, unsupervised and clinical analysis
pythonResults (Drive)                        figures/ and tables/
```

Data and results are **not** stored in this repository; they live in the project Drive folder (`MCBCI_analisisCancer/`):

| Folder | Content |
|---|---|
| `Data/rawData/` | `TCGA-BRCA_log_FPKM.tsv`, `TCGA-BRCA.clinical.tsv` (GDC), `TCGA-BRCA.survival.tsv.gz`, `MetaData.xlsx` (UCSC Xena) |
| `Data/reference/` | `Recon3D_genes.xlsx` (optional QC in notebook 01) |
| `results/processedData/` | outputs of notebook 01 |
| `matlabResults/` | outputs of the MATLAB scripts |
| `pythonResults/` | outputs of notebook 02 |

## How to run

1. `pip install -r requirements.txt`
2. Set the Drive path: `export BRCA_DRIVE_DIR=/path/to/MCBCI_analisisCancer` (or edit `DRIVE_DIR` in the first cell).
3. Run `BRCA_Metabolic_Flux_Pipeline_01.ipynb`. The first run queries MyGene.info (internet) and caches the Ensembl→Entrez map.
4. Build the models with `matlab/modelGen.m` using `results/processedData/xomics_files` and `param.transcriptomicThreshold = 0`, then compute fluxes with `matlab/BRCA_fluxPipeline.m` (COBRA Toolbox v3, Gurobi).
5. Run `BRCA_Metabolic_Flux_Pipeline_02.ipynb`.

## Repository layout

| Path | Content |
|---|---|
| `BRCA_Metabolic_Flux_Pipeline_01.ipynb` | data processing, expression → `genes,expVal` |
| `BRCA_Metabolic_Flux_Pipeline_02.ipynb` | supervised, unsupervised and clinical analysis |
| `matlab/` | model generation and flux calculation |
| `old/` | previous code, kept for traceability (see `old/README.md`, `old/inventory_2026-10-05.csv`) |
| `toBeDeleted/` | copies, session files and generated results waiting for review before deletion (see `toBeDeleted/README.md`) |

## Authors

Ruiz Robles E., Rincón-Ballesteros R., Chacón Méndez S.A., Alvarez-Padilla F.J., Preciat G. — Universidad de Guadalajara.
