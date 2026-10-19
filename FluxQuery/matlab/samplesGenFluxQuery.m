%% FluxQuery - Step 3: uniform sampling of the three neuron models
% Samples the feasible flux space of the dopaminergic (iDopaNeuroC), cholinergic
% and glutamatergic models with Coordinate Hit-and-Run with Rounding (CHRR,
% COBRA Toolbox). Three independent sets are generated per model; the first one
% computes the rounded polytope, which the other two reuse. A pilot set
% (smallSamples, randomly thinned from the pooled samples) is also saved for
% quick tests of step 4.
%
% Inputs (project folder, chosen with a dialog):
%   data/models/iDopaNeuroC.mat                         dopaminergic neuron (published)
%   results/matlab/models/<neuron>/<model>.mat          cholinergic and glutamatergic (step 2)
%
% Outputs, per neuron, in results/matlab/sampling/<neuron>/:
%   samples1.mat, samples2.mat, samples3.mat   (variables samples1..3: reactions x points)
%   smallSamples.mat                           (variable smallSamples: reactions x nSmall)
%   samplingReport.txt
%
% Run as a script, or open it as a Live Script (right click > Open as Live
% Script) and save as .mlx. Next step: samplingAnalysisFluxQuery.ipynb

%% Settings
clear
initCobraToolbox(false)

neurons = {'dopaminergic', 'cholinergic', 'glutamatergic'};
options.nFiles = 50;            % as in the original iDopaNeuro sampling
options.nPointsReturned = 8000;
nSets = 3;                      % independent sample sets per model
nSmall = 8000;                  % size of the pilot set
rng(123)                        % reproducibility

changeCobraSolver('gurobi', 'all', 0);

projectDir = uigetdir(pwd, 'Select the FluxQuery project folder (contains data/ and results/)');
if isequal(projectDir, 0), error('No project folder selected'); end

%% Sample each model
for i = 1:numel(neurons)
    neuron = neurons{i};
    fprintf('\n===== %s neuron =====\n', neuron);

    % model
    if strcmp(neuron, 'dopaminergic')
        modelFile = fullfile(projectDir, 'data', 'models', 'iDopaNeuroC.mat');
    else
        d = dir(fullfile(projectDir, 'results', 'matlab', 'models', neuron, 'i*Neuro.mat'));
        if isempty(d)
            warning('No model for the %s neuron: run modelGenFluxQuery.m first. Skipped.', neuron);
            continue
        end
        modelFile = fullfile(d(1).folder, d(1).name);
    end
    s = load(modelFile);
    if isfield(s, 'model'), model = s.model; else, fn = fieldnames(s); model = s.(fn{1}); end
    fprintf('model: %s (%d reactions, %d metabolites)\n', modelFile, numel(model.rxns), numel(model.mets));

    % feasibility check before sampling
    sol = optimizeCbModel(model);
    if sol.stat ~= 1
        warning('The %s model is not feasible (LP status %d). Skipped.', neuron, sol.stat);
        continue
    end

    outDir = fullfile(projectDir, 'results', 'matlab', 'sampling', neuron);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    workDir = fullfile(outDir, 'chrrTemp');
    if ~exist(workDir, 'dir'), mkdir(workDir); end
    oldDir = cd(workDir);                       % CHRR writes temporary files in the current folder

    % sampling: the first set computes the rounded polytope P, the others reuse it
    tic
    sets = cell(1, nSets);
    opt = options;
    opt.toRound = 1;
    [P, sets{1}] = sampleCbModel(model, [], 'CHRR', opt);
    opt.toRound = 0;
    for k = 2:nSets
        [~, sets{k}] = sampleCbModel(model, [], 'CHRR', opt, P);
    end
    elapsed = toc;
    cd(oldDir);
    rmdir(workDir, 's');

    % save (one variable per file, as in the original samplesGen.m)
    for k = 1:nSets
        varName = sprintf('samples%d', k);
        S.(varName) = sets{k};
        info = whos('-var', 'sets');
        if info.bytes / nSets > 1.9e9
            save(fullfile(outDir, [varName '.mat']), '-struct', 'S', varName, '-v7.3');
        else
            save(fullfile(outDir, [varName '.mat']), '-struct', 'S', varName);
        end
        S = rmfield(S, varName);
    end
    pooled = [sets{:}];
    smallSamples = pooled(:, randperm(size(pooled, 2), min(nSmall, size(pooled, 2))));
    save(fullfile(outDir, 'smallSamples.mat'), 'smallSamples');

    % report: steady state and bounds of the samples, agreement between sets
    res = max(abs(model.S * pooled), [], 1);
    viol = max(max(model.lb - pooled, 0) + max(pooled - model.ub, 0), [], 1);
    med = cellfun(@(x) median(x, 2), sets, 'UniformOutput', false);
    med = [med{:}];
    spread = max(med, [], 2) - min(med, [], 2);
    sorted = sort(pooled, 2);                    % interquartile range without toolboxes
    n = size(sorted, 2);
    iqrAll = sorted(:, max(1, round(0.75 * n))) - sorted(:, max(1, round(0.25 * n)));
    agree = mean(spread <= 0.1 * max(iqrAll, 1e-6));
    fid = fopen(fullfile(outDir, 'samplingReport.txt'), 'w');
    fprintf(fid, 'model: %s\nreactions: %d\nsets: %d x %d points (CHRR, nFiles = %d)\n', ...
        modelFile, numel(model.rxns), nSets, size(sets{1}, 2), options.nFiles);
    fprintf(fid, 'time: %.1f min\nmax |S v| over samples: %.3g\nmax bound violation: %.3g\n', ...
        elapsed / 60, max(res), max(viol));
    fprintf(fid, 'reactions whose median agrees between sets (within 10%% of the IQR): %.1f%%\n', 100 * agree);
    fclose(fid);
    fprintf('%s: %d x %d samples in %.1f min; median agreement between sets %.1f%%\n', ...
        neuron, nSets, size(sets{1}, 2), elapsed / 60, 100 * agree);
    clear sets pooled smallSamples P
end

disp('Samples saved in results/matlab/sampling. Next: samplingAnalysisFluxQuery.ipynb')
