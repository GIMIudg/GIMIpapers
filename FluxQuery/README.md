# FluxQuery code

Code of the FluxQuery project (axonal growth in dopaminergic, cholinergic and glutamatergic neurons).
Data and results are kept in a separate project folder with `data/` and `results/` (supplementary material);
every script asks for that folder when it starts.

Only code lives here. Inputs are read from `<project>/data` and `<project>/results`, outputs are written to
`<project>/results` (and, for step 1, `<project>/data/omics`). Each script asks for the project folder with a
dialog (Python: also the environment variable `FQ_PROJECT`, or set `PROJECT_DIR` in the first cell).

```
code/
├── python/
│   ├── preprocessingFluxQuery.ipynb      step 1: omics -> XomicsToModel input files
│   └── samplingAnalysisFluxQuery.ipynb   step 4: autoencoders, FluxQuery and medium design (3 neurons)
└── matlab/
    ├── modelGenFluxQuery.m               step 2: cholinergic and glutamatergic models (XomicsToModel)
    └── samplesGenFluxQuery.m             step 3: CHRR sampling of the three models
```

Run order: 1 → 2 → 3 → 4. Steps 2 and 3 can be run as scripts or opened as Live Scripts
(right click > *Open as Live Script*, then save as `.mlx`).

## Step 1 – `preprocessingFluxQuery.ipynb`

* Transcriptomics (FPKM, `genes,expVal`, Entrez IDs, mean of the control replicates):
  * cholinergic: GSE308143 (hiPSC-derived cholinergic neurons, NT1_C–NT5_C; Salmon counts → FPKM)
  * glutamatergic: GSE112732 (hiPSC-derived glutamatergic neurons, day 49, 6 samples; authors' FPKM)
  * dopaminergic reference (same study as the cholinergic one, only for comparison): GSE308143
* Exometabolomics (glutamatergic only, Metabolomics Workbench ST004403, iPSC-derived neurons,
  spent vs fresh medium). Converted to exchange rates as in the iDopaNeuro paper:
  `v = ΔC·V / (Δt·gDW)` with ΔC = (spent/fresh − 1)·C_medium, V = 2 mL, Δt = 48 h and
  gDW from total protein; mean ± SD in µmol/gDW/h. Metabolites that are released but absent from
  the medium only get a direction (`exometabolomicDirection_glutamatergicNeuron.txt`).
  Glutamine and alanine are excluded (the medium contains the alanyl-glutamine dipeptide).
* No exometabolomic dataset exists for human cholinergic neurons (checked October 2026).

Downloads are cached in `data/omics/<neuron>/raw`, so the notebook runs offline once they exist.
Requirements: Python ≥ 3.10, numpy, pandas, requests.

## Step 2 – `modelGenFluxQuery.m`

XomicsToModel (COBRA Toolbox) with the iDopaNeuro options, from Recon3DModel_301_xomics_input.mat
(`CBTDIR/papers/2023_iDopaNeuro`, or chosen with a dialog), the bibliomic workbook
`results/bibliomic/bibliomicData_<neuron>Neuron.xlsx` and the omics files of step 1.
The glutamatergic model also uses the exometabolomic rates and adds the secreted metabolites to the
core set. Output: `results/matlab/models/<neuron>/iCholiNeuro.mat`, `iGlutaNeuro.mat` and the log.

## Step 3 – `samplesGenFluxQuery.m`

CHRR (`sampleCbModel`, nFiles = 50, nPointsReturned = 8000) for iDopaNeuroC, iCholiNeuro and
iGlutaNeuro. Three independent sets per model (the rounding of set 1 is reused), a random pilot set
of 8,000 states (`smallSamples.mat`) and `samplingReport.txt` (steady state, bounds, agreement
between sets). Output: `results/matlab/sampling/<neuron>/`.

Requirements (steps 2–3): MATLAB, COBRA Toolbox, Gurobi (or another LP/MILP solver).

## Step 4 – `samplingAnalysisFluxQuery.ipynb`

For each neuron with a model and a sampling: conversion of the `.mat` samples to NumPy (once),
removal of thermodynamically infeasible cycles, translation of the axonal-growth query
(`results/axonalGrowthQuery/axonalGrowthQuery_v4.docx`) into flux thresholds, training of an
autoencoder (AE) and a physics-informed autoencoder (PI-AE, penalised by ‖S·v‖² and bound violations), FluxQuery latent search, exact LP
repair, the growth ensemble, comparison with classical objectives (pFBA), and the LP design of a
growth medium. Section 7 compares the neurons. `SAMPLING_SET = "pilot"` uses `smallSamples.mat`
for a quick run; `"full"` uses `samples1-3.mat`.
Requirements: Python ≥ 3.10, numpy, scipy (≥ 1.9, HiGHS LP solver), pandas, matplotlib, scikit-learn,
python-docx, openpyxl. The autoencoders are implemented in NumPy (no deep-learning framework needed).

> Note: the MATLAB scripts were written against the COBRA Toolbox API but could not be executed
> where they were prepared (no MATLAB). Run them first on a small test and report any error.
