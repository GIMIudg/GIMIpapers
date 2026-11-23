%% FluxQuery - Step 2: generation of the dopaminergic, cholinergic and glutamatergic models
% Builds context-specific genome-scale models of human dopaminergic (iDopaNeuro),
% cholinergic (iCholiNeuro) and glutamatergic (iGlutaNeuro) neurons with
% XomicsToModel (COBRA Toolbox), following the iDopaNeuro pipeline. The three
% models are built by this script with the same parameters; the dopaminergic
% model uses the bibliomic, transcriptomic and exometabolomic data of the
% published iDopaNeuroC (data/models/iDopaNeuroC.mat), with the corrections of
% results/bibliomic/bibliomicCurationLog.xlsx.
%
% Inputs (project folder, chosen with a dialog; nothing is hard-coded):
%   results/bibliomic/bibliomicData_<neuron>Neuron.xlsx      bibliomic data (step 0, manual curation)
%   data/omics/<neuron>/transcriptomicData_<neuron>Neuron.txt  genes,expVal (FPKM)   (step 1)
%   data/omics/<neuron>/exometabolomicData_<neuron>Neuron.txt   rxns,...,mean,SD (umol/gDW/h), when measured
%   data/models/iDopaNeuroC.mat   source of the dopaminergic transcriptomic and exometabolomic
%                                 files, which are written to data/omics/dopaminergic/ when missing
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

neurons = {'dopaminergic', 'cholinergic', 'glutamatergic'};     % neurons to build
modelNames = containers.Map(neurons, {'iDopaNeuro', 'iCholiNeuro', 'iGlutaNeuro'});
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
param.TolMinBoundary = -1e4;                    % generic flux bounds of iDopaNeuro (limitBoundary 10000);
param.TolMaxBoundary = 1e4;                     % the XomicsToModel default (1000) truncates O2, CO2 and glucose
param.printLevel = 1;

%% Dopaminergic inputs from the published iDopaNeuroC (written once to data/omics/dopaminergic)
dopaDir = fullfile(dataDir, 'omics', 'dopaminergic');
dopaFiles = {fullfile(dopaDir, 'transcriptomicData_dopaminergicNeuron.txt'), ...
    fullfile(dopaDir, 'exometabolomicData_dopaminergicNeuron.txt')};
if ismember('dopaminergic', neurons) && ~all(cellfun(@isfile, dopaFiles))
    idn = load(fullfile(dataDir, 'models', 'iDopaNeuroC.mat'));
    sd = idn.iDopaNeuroC.XomicsToModelSpecificData;
    tr = sd.transcriptomicData;
    tr = table(string(tr.genes), tr.expVal, 'VariableNames', {'genes', 'expVal'});
    writetable(tr, dopaFiles{1});
    exo = sd.exoMet;                                    % iDopaNeuro names: rxnID, metID
    exo = renamevars(exo, intersect({'rxnID', 'metID'}, exo.Properties.VariableNames, 'stable'), ...
        replace(intersect({'rxnID', 'metID'}, exo.Properties.VariableNames, 'stable'), {'rxnID', 'metID'}, {'rxns', 'mets'}));
    exo = exo(:, intersect({'rxns', 'mets', 'metName', 'platform', 'mean', 'SD'}, exo.Properties.VariableNames, 'stable'));
    exo = exo(~cellfun(@isempty, cellstr(string(exo.rxns))), :);
    writetable(exo, dopaFiles{2});
    fprintf('dopaminergic transcriptomic (%d genes) and exometabolomic (%d rates) data written from iDopaNeuroC\n', ...
        height(tr), height(exo));
end

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
        for v = intersect({'rxns', 'mets', 'metName', 'platform'}, specificData.exoMet.Properties.VariableNames)
            specificData.exoMet.(v{1}) = cellstr(specificData.exoMet.(v{1}));   % text columns as cell arrays
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

    % XomicsToModel does not read the rxnsHypothesis sheet: its bounds (uptake hypotheses
    % with capped rates, required secretions) are applied as constrained reactions,
    % unless the reaction is already in rxns2constrain.
    if isfield(specificData, 'rxnsHypothesis') && istable(specificData.rxnsHypothesis) && ...
            ~isempty(specificData.rxnsHypothesis)
        H = specificData.rxnsHypothesis;
        H.rxns = strtrim(cellstr(H.rxns));
        % demands and sinks are not passed: thermoKernel does not accept them as constrained
        % reactions (they enter through sinkDemand)
        H = H(~ismember(H.rxns, cellstr(specificData.rxns2constrain.rxns)) & ismember(H.rxns, genericModel.rxns) & ...
            ~startsWith(H.rxns, {'DM_', 'sink_'}), :);
        T = specificData.rxns2constrain;
        add = T(ones(height(H), 1), :);
        for v = T.Properties.VariableNames
            x = add.(v{1});
            if iscell(x), add.(v{1}) = repmat({''}, height(add), 1);
            elseif isstring(x), add.(v{1}) = repmat("", height(add), 1);
            elseif isnumeric(x), add.(v{1}) = NaN(height(add), 1);
            end
        end
        add.rxns = H.rxns;
        add.lb = H.lb;
        add.ub = H.ub;
        if ismember('constraintDescription', T.Properties.VariableNames)
            add.constraintDescription = repmat({'Uptake hypothesis (rxnsHypothesis)'}, height(add), 1);
        end
        specificData.rxns2constrain = [T; add];
        fprintf('%d reactions of rxnsHypothesis applied as constrained reactions\n', height(add));
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
        'EX_na1[e]', 'EX_k[e]', 'EX_cl[e]', 'EX_ca2[e]', 'EX_fe2[e]', 'EX_fe3[e]', 'EX_i[e]'};   % H2O2: secretion only
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
        add.ub = param.TolMaxBoundary * ones(numel(secretionOnly), 1);
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

disp('Models saved in results/matlab/models. Next: python/modelValidationFluxQuery.py, then samplesGenFluxQuery.m')
