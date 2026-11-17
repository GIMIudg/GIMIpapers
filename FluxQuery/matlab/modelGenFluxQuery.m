%% FluxQuery - Step 2: generation of the cholinergic and glutamatergic models
% Builds context-specific genome-scale models of human cholinergic and
% glutamatergic neurons with XomicsToModel (COBRA Toolbox), following the
% iDopaNeuro pipeline. The dopaminergic neuron is NOT rebuilt: the published
% iDopaNeuroC model (data/models/iDopaNeuroC.mat) is used in steps 3 and 4.
%
% Inputs (project folder, chosen with a dialog; nothing is hard-coded):
%   results/bibliomic/bibliomicData_<neuron>Neuron.xlsx      bibliomic data (step 0, manual curation)
%   data/omics/<neuron>/transcriptomicData_<neuron>Neuron.txt  genes,expVal (FPKM)   (step 1)
%   data/omics/glutamatergic/exometabolomicData_glutamatergicNeuron.txt       rxns,...,mean,SD (umol/gDW/h)
%   data/omics/glutamatergic/exometabolomicDirection_glutamatergicNeuron.txt  secreted/uptaken, no rate
%   generic model: Recon3DModel_301_xomics_input.mat (COBRA Toolbox, papers/2023_iDopaNeuro)
%
% Outputs:
%   results/matlab/models/<neuron>/i<Name>Neuro.mat     model (variable 'model') and report
%   results/matlab/models/<neuron>/modelGenerationReport.txt
%
% Requirements: MATLAB, COBRA Toolbox (with XomicsToModel), an LP/MILP solver
% (Gurobi recommended). Run as a script, or open it as a Live Script
% (right click > Open as Live Script) and save as .mlx.
%
% Next step: samplesGenFluxQuery.m

%% Settings
clear
initCobraToolbox(false)
global CBTDIR

neurons = {'cholinergic', 'glutamatergic'};     % neurons to build
modelNames = containers.Map(neurons, {'iCholiNeuro', 'iGlutaNeuro'});
solver = 'gurobi';

% Exometabolomic data are used when present (glutamatergic). Metabolites that
% are clearly secreted but cannot be given an absolute rate are added to the
% core set (their exchange must be able to carry flux) when this is true:
useSecretedAsActive = true;

%% Project folder and generic model
projectDir = uigetdir(pwd, 'Select the FluxQuery project folder (contains data/ and results/)');
if isequal(projectDir, 0), error('No project folder selected'); end
dataDir = fullfile(projectDir, 'data');
resultsDir = fullfile(projectDir, 'results');

genericModelFile = fullfile(CBTDIR, 'papers', '2023_iDopaNeuro', 'Recon3DModel_301_xomics_input.mat');
if ~isfile(genericModelFile)
    [f, p] = uigetfile('*.mat', 'Select the generic model (Recon3DModel_301_xomics_input.mat)');
    if isequal(f, 0), error('No generic model selected'); end
    genericModelFile = fullfile(p, f);
end
generic = load(genericModelFile);
fn = fieldnames(generic);
genericModel = generic.(fn{1});
if isfield(genericModel, 'expressionRxns')
    genericModel = rmfield(genericModel, 'expressionRxns');
end
fprintf('Generic model: %s (%d reactions)\n', genericModelFile, numel(genericModel.rxns));

changeCobraSolver(solver, 'all', 0);

% Use the COBRA Toolbox XomicsToModel: other copies on the path (e.g. in COBRA.papers)
% may be older versions that fail on the transcriptomic data (model.geneExpVal)
xomicsDir = fullfile(CBTDIR, 'src', 'dataIntegration', 'XomicsToModel');
addpath(xomicsDir, '-begin');
allCopies = which('XomicsToModel', '-all');
if numel(allCopies) > 1
    fprintf('Other copies of XomicsToModel on the path (ignored):\n');
    fprintf('   %s\n', allCopies{2:end});
end
fprintf('Using %s\n', which('XomicsToModel'));

%% XomicsToModel options (as for iDopaNeuro)
feasTol = getCobraSolverParams('LP', 'feasTol');
param = struct();
param.setObjective = '';                        % no objective function
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';              % weights from the exometabolomic mean
param.sinkDMinactive = 1;                       % non-core sinks and demands inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true;                      % only uptakes of medium components (mediaData)
param.findThermoConsistentFluxSubset = 1;
param.tissueSpecificSolver = 'thermoKernel';
param.activeGenesApproach = 'oneRxnPerActiveGene';
param.transcriptomicThreshold = 0;              % log2(FPKM): genes with FPKM > 1 are expressed
param.printLevel = 1;

%% Build each model
for i = 1:numel(neurons)
    neuron = neurons{i};
    name = modelNames(neuron);
    fprintf('\n===== %s neuron (%s) =====\n', neuron, name);
    outDir = fullfile(resultsDir, 'matlab', 'models', neuron);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % bibliomic data
    bibliomicFile = fullfile(resultsDir, 'bibliomic', ['bibliomicData_' neuron 'Neuron.xlsx']);
    specificData = preprocessingOmicsModel(bibliomicFile, 1, 1);

    % transcriptomic data (genes, expVal)
    omicsDir = fullfile(dataDir, 'omics', neuron);
    transcriptomic = readtable(fullfile(omicsDir, ['transcriptomicData_' neuron 'Neuron.txt']));
    transcriptomic.genes = string(transcriptomic.genes);
    specificData.transcriptomicData = transcriptomic;
    fprintf('transcriptomic data: %d genes\n', height(transcriptomic));

    % exometabolomic data (rxns, mets, metName, platform, mean, SD), when available
    exoFile = fullfile(omicsDir, ['exometabolomicData_' neuron 'Neuron.txt']);
    if isfile(exoFile)
        specificData.exoMet = readtable(exoFile, 'TextType', 'string');
        for v = {'rxns', 'mets', 'metName', 'platform'}                % text columns as cell arrays
            specificData.exoMet.(v{1}) = cellstr(specificData.exoMet.(v{1}));
        end
        fprintf('exometabolomic data: %d exchange rates\n', height(specificData.exoMet));
    else
        fprintf('no exometabolomic data for this neuron\n');
    end

    % secreted metabolites without an absolute rate -> core (active) reactions
    dirFile = fullfile(omicsDir, ['exometabolomicDirection_' neuron 'Neuron.txt']);
    if useSecretedAsActive && isfile(dirFile)
        direction = readtable(dirFile, 'TextType', 'string');
        use = direction.direction == "secreted";
        if ismember('useInModel', direction.Properties.VariableNames)   % false: not producible by Recon3D (see note)
            use = use & (string(direction.useInModel) == "True" | string(direction.useInModel) == "true" | string(direction.useInModel) == "1");
        end
        secreted = cellstr(direction.rxns(use));
        secreted = secreted(ismember(secreted, genericModel.rxns));
        if ~isfield(specificData, 'activeReactions') || isempty(specificData.activeReactions)
            specificData.activeReactions = secreted(:);
        elseif istable(specificData.activeReactions)              % table with a 'rxns' column
            specificData.activeReactions = unique([cellstr(specificData.activeReactions.rxns); secreted(:)]);
        else                                                          % cell array of reaction IDs
            specificData.activeReactions = unique([specificData.activeReactions(:); secreted(:)]);
        end
        fprintf('%d secreted exchanges added to the core set\n', numel(secreted));
    end

    % uptake only of medium components: core exchange reactions keep the open bounds of the
    % generic model, so exchanges that are core but not medium components (e.g. secretions
    % observed in the exometabolome or inherited from iDopaNeuro) are made secretion-only.
    % Gases, water, protons and inorganic ions stay reversible; uptake hypotheses
    % (rxnsHypothesis) and explicit constraints (rxns2constrain) are kept as given.
    if istable(specificData.activeReactions)
        coreRxns = cellstr(specificData.activeReactions.rxns);
    else
        coreRxns = cellstr(specificData.activeReactions);
    end
    keepOpen = {'EX_o2[e]', 'EX_co2[e]', 'EX_h2o[e]', 'EX_h[e]', 'EX_hco3[e]', 'EX_nh4[e]', 'EX_pi[e]', ...
        'EX_na1[e]', 'EX_k[e]', 'EX_cl[e]', 'EX_ca2[e]', 'EX_fe2[e]', 'EX_fe3[e]', 'EX_h2o2[e]', 'EX_i[e]'};
    given = {};
    for f = {'mediaData', 'rxnsHypothesis', 'rxns2constrain'}
        if isfield(specificData, f{1}) && istable(specificData.(f{1})) && ...
                ismember('rxns', specificData.(f{1}).Properties.VariableNames)
            given = [given; cellstr(specificData.(f{1}).rxns)]; %#ok<AGROW>
        end
    end
    secretionOnly = coreRxns(startsWith(coreRxns, 'EX_') & ~ismember(coreRxns, [keepOpen(:); given]));
    secretionOnly = unique(secretionOnly(ismember(secretionOnly, genericModel.rxns)));
    if ~isempty(secretionOnly)
        T = specificData.rxns2constrain;
        add = T(ones(numel(secretionOnly), 1), :);
        for v = T.Properties.VariableNames
            x = add.(v{1});
            if iscell(x), add.(v{1}) = repmat({''}, height(add), 1);
            elseif isstring(x), add.(v{1}) = repmat("", height(add), 1);
            elseif isnumeric(x), add.(v{1}) = NaN(height(add), 1);
            end
        end
        add.rxns = secretionOnly(:);
        add.lb = zeros(numel(secretionOnly), 1);
        add.ub = 1000 * ones(numel(secretionOnly), 1);
        if ismember('constraintDescription', T.Properties.VariableNames)
            add.constraintDescription = repmat({'Secretion only: not a medium component'}, height(add), 1);
        end
        specificData.rxns2constrain = [T; add];
        fprintf('%d core exchanges of non-medium metabolites made secretion-only\n', numel(secretionOnly));
    end

    % model extraction
    param.workingDirectory = outDir;
    param.diaryFilename = fullfile(outDir, 'modelGenerationReport.txt');
    [model, modelGenerationReport] = XomicsToModel(genericModel, specificData, param);
    model.modelID = name;
    model.description = sprintf('%s: context-specific model of a human %s neuron (FluxQuery, XomicsToModel, %s)', ...
        name, neuron, char(datetime('today', 'Format', 'yyyy-MM-dd')));

    % checks
    sol = optimizeCbModel(model);
    fprintf('%s: %d reactions, %d metabolites, %d genes; LP status %d\n', name, ...
        numel(model.rxns), numel(model.mets), numel(model.genes), sol.stat);

    save(fullfile(outDir, [name '.mat']), 'model', 'modelGenerationReport', 'specificData', 'param')
    fprintf('saved %s\n', fullfile(outDir, [name '.mat']));
end

disp('Models saved in results/matlab/models. Next: samplesGenFluxQuery.m')
