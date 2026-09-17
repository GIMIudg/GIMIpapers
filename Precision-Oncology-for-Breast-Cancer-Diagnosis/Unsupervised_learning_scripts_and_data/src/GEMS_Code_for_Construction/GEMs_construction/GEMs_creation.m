% ============================================================
% GEMs_creation.m
% Description:
%   Generates patient-specific genome-scale metabolic models (GEMs) by
%   integrating TCGA-BRCA transcriptomic data with bibliomic information
%   and the generic Recon3D model using the thermoKernel algorithm.
%
% Inputs (resolved via Clinical_data_and_models_ids/):
%   - Recon3D_301.mat              : Generic GEM reference model
%   - TCGA-BRCA_log_FPKM.tsv      : RNA-seq expression matrix
%   - Bibliomic_Gems.xlsx          : Manually curated bibliomic constraints
%   - Xomics_files/*.txt           : Per-sample Xomics files (from creator.py)
%
% Output:
%   - All_models_created/*.mat     : 1226 patient-specific GEMs
%
% Requirements: COBRA Toolbox, Gurobi solver
% ============================================================

% Initialize CobraToolbox and solver
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');
%% 

% ── Resolve paths relative to this script's location ────────────────────────
% Hierarchy: GEMs_construction/ -> GEMS_Code_for_Construction/ -> src/ -> Unsupervised_learning_scripts_and_data/
scriptDir    = fileparts(mfilename('fullpath'));
dataRoot     = fullfile(scriptDir, '..', '..', '..', 'Clinical_data_and_models_ids');

recon3dPath       = fullfile(dataRoot, 'GEMs_Data_for_construction', 'Recon3D_301.mat');
XomicsFilesPath   = fullfile(dataRoot, 'GEMs_Data_for_construction', 'Xomics_files');
bibliomicDataPath = fullfile(dataRoot, 'GEMs_Data_for_construction', 'Bibliomic_Gems.xlsx');
outputDir         = fullfile(dataRoot, 'All_models_created');

% ── Load generic reference model ────────────────────────────────────────────
model = readCbModel(recon3dPath);

% Strip version suffix from gene IDs (e.g. "1234.1" → "1234")
for i = 1:numel(model.genes)
    parts = split(model.genes{i}, '.');
    model.genes{i} = parts{1};
end

% ── Load Xomics files ────────────────────────────────────────────────────────
% Get all per-sample Xomics files
XomicsFilesInfo  = dir(fullfile(XomicsFilesPath, '*.txt'));

% Build full path list
XomicsFilesNames = {};
for i = 1:length(XomicsFilesInfo)
    XomicsFileName   = XomicsFilesInfo(i).name;
    fullPath         = fullfile(XomicsFilesPath, XomicsFileName);
    XomicsFilesNames{i,1} = fullPath;

    % Read table to verify readability
    T = readtable(fullPath);
end

% Get file names without .txt extension (used as model identifiers)
filesNames = {};
for i = 1:length(XomicsFilesInfo)
    [~, XomicsFileName, ~] = fileparts(XomicsFilesInfo(i).name);
    filesNames{i,1} = XomicsFileName;
end


%% Load bibliomic data ───────────────────────────────────────────────────────
specificData = preprocessingOmicsModel(bibliomicDataPath);

% Create output directory if it does not exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Transcriptomic Integration ────────────────────────────────────────────────

% Clean gene IDs: extract numeric part only
genes   = cell(numel(model.genes), 1);
pattern = '\d+';                          % Pattern to extract numeric gene ID part
for i = 1:numel(model.genes)
    match     = regexp(model.genes{i}, pattern, 'match', 'once');
    genes{i}  = match;
end
model.genes = genes;

% Clean grRules: remove version suffixes and special characters
grRules = cell(numel(model.grRules), 1);
pattern = '[_,.#/-]\w*';
for i = 1:numel(model.grRules)
    grRules{i} = regexprep(model.grRules{i}, pattern, '');
end
model.grRules = grRules;

model.biomassRxnAbbr = 'biomass_ob';

% XomicsToModel integration parameters
param.inactiveGenesTranscriptomics = true;
param.curationOverOmics            = false;
param.metabolomicWeights           = 'mean';
param.transcriptomicThreshold      = 0.05;
param.weightsFromOmics             = true;
param.fluxCCmethod                 = 'fastcc';
param.printLevel                   = 3;
param.debug                        = false;
param.diaryFilename                = 0;
param.verbose                      = false;

getCobraSolverParams('LP', 'feasTol', 1e-6);
feasTol = getCobraSolverParams('LP', 'feasTol');

param.allowRxnRelaxation  = true;
param.TolMinBoundary      = -1e4;
param.TolMaxBoundary      =  1e4;
param.boundPrecisionLimit = feasTol * 100;
param.fluxEpsilon         = feasTol * 100;

param.closeIons              = false;
param.closeUptakes           = false;
param.nonCoreSinksDemands    = 'closeAll';
param.sinkDMinactive         = false;
param.allowNetFlux           = true;
param.protectEssentialRxns   = true;
param.protectRxns            = {'biomass_ob','biomass_ex'};
param.activeGenesApproach    = 'oneRxnPerActiveGene';
param.modelExtractionAlgorithm = 'thermoKernel';
param.coreFromActiveGenes    = true;
param.core                   = {'biomass_ob','biomass_ex'};
params.specificRxns          = {'biomass_ob', 'biomass_ex'};

%% Generate patient-specific GEMs ────────────────────────────────────────────
for i = 1:length(XomicsFilesNames)
    specificData.transcriptomicData = readtable(XomicsFilesNames{i});
    [specificModel, ~]  = XomicsToModel(model, specificData, param);
    specificModelPath   = fullfile(outputDir, [filesNames{i}, '_specificModel.mat']);
    delete(fullfile(pwd, '*thermoKernel*'));  % Clean up temporary thermoKernel files
    save(specificModelPath, "specificModel");
    disp(['✅ Model ', filesNames{i}, ' saved successfully']);
end

disp("✅ All patient-specific GEMs generated successfully.");
