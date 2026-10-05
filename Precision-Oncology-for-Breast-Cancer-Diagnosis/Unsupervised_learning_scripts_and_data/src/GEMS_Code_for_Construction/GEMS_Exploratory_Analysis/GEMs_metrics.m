% ============================================================
% GEMs_metrics.m
% Purpose: Computes reaction-level metrics and pairwise similarity
%          across all 1226 patient-specific thermoKernel GEMs.
%
% Outputs (saved to All_models_created/ folder):
%   - Resumen_Modelos_ThermoKernel.csv     : Per-model summary (reactions + similarity)
%   - Intersection_Matrix_ThermoKernel.csv : Pairwise reaction intersection matrix (%)
%
% Requirements: COBRA Toolbox, Parallel Computing Toolbox (optional, for parfor)
% ============================================================


% Initialize CobraToolbox and solver
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');

% ── Resolve paths relative to this script's location ────────────────────────
% Hierarchy: GEMS_Exploratory_Analysis/ -> GEMS_Code_for_Construction/ -> src/ -> Unsupervised_learning_scripts_and_data/
scriptDir    = fileparts(mfilename('fullpath'));
dataRoot     = fullfile(scriptDir, '..', '..', '..', 'Clinical_data_and_models_ids');
modelsFolder = fullfile(dataRoot, 'All_models_created');
outputDir    = modelsFolder;   % Save results alongside the models

% ── Load all thermoKernel model files ────────────────────────────────────────
filesTK = dir(fullfile(modelsFolder, '*.mat'));
filesTK = filesTK(~[filesTK.isdir]);   % Exclude subdirectories

numFiles     = length(filesTK);
modelsThermo = cell(numFiles, 1);
rxnsThermo   = cell(numFiles, 1);
modelNames   = cell(numFiles, 1);

% Load models with error handling
fprintf('Loading %d models...\n', numFiles);
for i = 1:numFiles
    try
        modelTK        = readCbModel(fullfile(modelsFolder, filesTK(i).name));
        modelsThermo{i} = modelTK;
        rxnsThermo{i}  = modelTK.rxns;
        modelNames{i}  = erase(filesTK(i).name, '.mat');
    catch ME
        warning('Error loading model %s: %s', filesTK(i).name, ME.message);
        modelsThermo{i} = [];
        rxnsThermo{i}  = {};
        modelNames{i}  = erase(filesTK(i).name, '.mat');
    end
end

% ── Compute pairwise reaction intersection matrix ─────────────────────────────
interMatrix = zeros(numFiles);

fprintf('Computing intersection matrix...\n');
% Use parfor to speed up computation (requires Parallel Computing Toolbox)
parfor i = 1:numFiles
    tempRow = zeros(1, numFiles);
    for j = i:numFiles   % Upper triangle only (symmetric matrix)
        if isempty(rxnsThermo{i}) || isempty(rxnsThermo{j})
            tempRow(j) = 0;
        else
            common = intersect(rxnsThermo{i}, rxnsThermo{j});
            denom  = mean([length(rxnsThermo{i}), length(rxnsThermo{j})]);
            tempRow(j) = (denom > 0) * length(common) / max(denom, eps);
        end
    end
    interMatrix(i,:) = tempRow;
end

% Symmetrize: fill lower triangle from upper
for i = 2:numFiles
    for j = 1:i-1
        interMatrix(i,j) = interMatrix(j,i);
    end
end

% ── Figure 1: Intersection heatmap ───────────────────────────────────────────
figure;
imagesc(interMatrix * 100);
colormap(jet);
colorbar;
title('Reaction Intersection (%) Between ThermoKernel Models');
xlabel('Model');
ylabel('Model');
axis square;

% Manage tick labels: show all if ≤50 models, otherwise subsample
if numFiles <= 50
    xticks(1:numFiles); yticks(1:numFiles);
    xticklabels(modelNames); yticklabels(modelNames);
    xtickangle(45);
else
    % Show a subset of ticks to avoid label overlap
    stepTick = ceil(numFiles / 20);
    ticks    = 1:stepTick:numFiles;
    xticks(ticks); yticks(ticks);
    xticklabels(modelNames(ticks)); yticklabels(modelNames(ticks));
    xtickangle(45);
end

% ── Figure 2: Number of reactions per model (bar chart) ──────────────────────
numRxnsTK = cellfun(@length, rxnsThermo);

figure;
bar(numRxnsTK, 'FaceColor', [0.2 0.6 0.8]);
title('Number of Reactions per ThermoKernel Model');
xlabel('Model');
ylabel('Number of reactions');
if numFiles <= 50
    xticks(1:numFiles); xticklabels(modelNames);
else
    stepTick = ceil(numFiles / 20);
    ticks    = 1:stepTick:numFiles;
    xticks(ticks); xticklabels(modelNames(ticks));
end
xtickangle(45);
grid on;

% ── Figure 3: Reactions per model sorted descending ──────────────────────────
[~, sortIdx]  = sort(numRxnsTK, 'descend');
sortedNames   = modelNames(sortIdx);
sortedNumRxns = numRxnsTK(sortIdx);

figure;
bar(sortedNumRxns, 'FaceColor', [0.8 0.4 0.4]);
title('Number of Reactions per ThermoKernel Model (sorted)');
xlabel('Model (sorted)');
ylabel('Number of reactions');
if numFiles <= 50
    xticks(1:numFiles); xticklabels(sortedNames);
else
    stepTick = ceil(numFiles / 20);
    ticks    = 1:stepTick:numFiles;
    xticks(ticks); xticklabels(sortedNames(ticks));
end
xtickangle(45);
grid on;

% ── Figure 4: Histogram of pairwise intersection values ──────────────────────
upperTri = triu(interMatrix, 1);        % Upper triangle only (no diagonal)
vals     = upperTri(upperTri > 0);      % Exclude unfilled zeros

figure;
histogram(vals * 100, 20, 'FaceColor', [0.2 0.8 0.4], 'EdgeColor', 'k');
title('Distribution of Reaction Intersection (%) Between ThermoKernel Models');
xlabel('Intersection (%)');
ylabel('Frequency');
grid on;

% ── Figure 5: Average similarity per model ───────────────────────────────────
avgSimilarity = mean(interMatrix, 2, 'omitnan');

figure;
bar(avgSimilarity * 100, 'FaceColor', [0.4 0.4 0.8]);
title('Average Similarity of Each Model Relative to All Others');
xlabel('Model');
ylabel('Average similarity (%)');
if numFiles <= 50
    xticks(1:numFiles); xticklabels(modelNames);
else
    stepTick = ceil(numFiles / 20);
    ticks    = 1:stepTick:numFiles;
    xticks(ticks); xticklabels(modelNames(ticks));
end
xtickangle(45);
grid on;

% ── Export summary table ──────────────────────────────────────────────────────
T = table(modelNames(:), numRxnsTK(:), avgSimilarity(:), ...
    'VariableNames', {'Model', 'NumReactions', 'AverageSimilarity'});
writetable(T, fullfile(outputDir, 'Resumen_Modelos_ThermoKernel.csv'));

fprintf('✅ Analysis complete. Summary saved to Resumen_Modelos_ThermoKernel.csv\n');

% ── Export pairwise intersection matrix ──────────────────────────────────────
intersectionTable = array2table(interMatrix * 100, ...   % values in percent
    'VariableNames', modelNames, ...
    'RowNames',      modelNames);

writetable(intersectionTable, ...
    fullfile(outputDir, 'Intersection_Matrix_ThermoKernel.csv'), ...
    'WriteRowNames', true);

fprintf('✅ Intersection matrix exported to Intersection_Matrix_ThermoKernel.csv\n');
