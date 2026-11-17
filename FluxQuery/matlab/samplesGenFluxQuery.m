%% FluxQuery - Step 3: reproducible uniform sampling of the three neuron models
% Samples the feasible flux space of the dopaminergic (iDopaNeuroC), cholinergic
% (iCholiNeuro) and glutamatergic (iGlutaNeuro) models with Coordinate Hit-and-Run
% with Rounding (CHRR, COBRA Toolbox).
%
% REPRODUCIBILITY
% Every random number used by CHRR comes from MATLAB's global generator, which is
% reset with a fixed seed (Mersenne Twister) right before each sample set. Set k of
% neuron i always uses seed = masterSeed + 100*i + k, so any set can be regenerated
% on its own, in any order. The rounding of the polytope (set 1) is deterministic
% (flux variability + maximum-volume ellipsoid, no random numbers) and is reused by
% sets 2 and 3.
% Because the samples can be regenerated exactly, they are NOT part of the
% supplementary material. Instead, each neuron gets two small files that are:
%   samplingManifest.json   seeds, settings, software versions, SHA-256 of the
%                           model file and of every sample set
%   samplingSummary.csv     mean, SD and percentiles of every reaction
% When samplingManifest.json already exists (e.g. copied from the supplementary
% material), a new run compares its samples with it and reports whether they are
% bit-identical. Bit-identical results require the same model file, MATLAB release,
% COBRA Toolbox version and LP solver; on other platforms the samples are a new
% draw from the same distribution, and samplingSummary.csv shows how close they are.
%
% Inputs (project folder, chosen with a dialog):
%   data/models/iDopaNeuroC.mat                       dopaminergic neuron (published)
%   results/matlab/models/<neuron>/i*Neuro.mat        cholinergic and glutamatergic (step 2)
%
% Outputs, per neuron, in results/matlab/sampling/<neuron>/:
%   samples1.mat, samples2.mat, samples3.mat   variables samples1..3 (reactions x points)
%   smallSamples.mat                           variable smallSamples (pilot set)
%   samplingManifest.json, samplingSummary.csv, samplingReport.txt
%
% Run as a script, or open it as a Live Script (right click > Open as Live Script)
% and save as .mlx. Next step: code/python/samplingAnalysisFluxQuery.ipynb

%% Settings
clear
initCobraToolbox(false)
global CBTDIR

neurons = {'dopaminergic', 'cholinergic', 'glutamatergic'};   % neurons to sample
masterSeed = 20261113;          % fixed; do not change once results are reported
nSets = 3;                      % independent sample sets per model
options.nPointsReturned = 8000; % points per set
options.nStepsPerPoint = 2000;  % CHRR steps between stored points (200, the COBRA default, did not converge)
nSmall = 8000;                  % size of the pilot set (random subset of the pooled sets)
overwrite = false;              % false: skip a neuron whose samples match the current model, seeds and settings
solver = 'gurobi';

changeCobraSolver(solver, 'all', 0);

projectDir = uigetdir(pwd, 'Select the FluxQuery project folder (contains data/ and results/)');
if isequal(projectDir, 0), error('No project folder selected'); end

% software versions, recorded in every manifest
env.matlab = version;
env.platform = computer;
env.solver = solver;
try
    env.solverVersion = char(string(getCobraSolverVersion(solver, 0)));
catch
    env.solverVersion = 'unknown';
end
[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', CBTDIR));
if status == 0, env.cobraToolboxCommit = strtrim(commit); else, env.cobraToolboxCommit = 'unknown'; end

%% Sample each model
for i = 1:numel(neurons)
    neuron = neurons{i};
    fprintf('\n===== %s neuron =====\n', neuron);
    outDir = fullfile(projectDir, 'results', 'matlab', 'sampling', neuron);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    setFiles = arrayfun(@(k) fullfile(outDir, sprintf('samples%d.mat', k)), 1:nSets, 'UniformOutput', false);

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
    % skip only if these exact samples (same model file, seeds and settings) already exist
    seeds = masterSeed + 100 * i + (1:nSets);
    manifestFile = fullfile(outDir, 'samplingManifest.json');
    if ~overwrite && all(cellfun(@isfile, setFiles)) && isfile(manifestFile)
        old = jsondecode(fileread(manifestFile));
        if strcmp(old.modelSHA256, sha256File(modelFile)) && isequal(old.seeds(:)', seeds) ...
                && old.nPointsReturned == options.nPointsReturned && old.nStepsPerPoint == options.nStepsPerPoint
            fprintf('Samples of this model with these seeds and settings already exist. Skipped.\n');
            continue
        end
    end

    s = load(modelFile);
    if isfield(s, 'model'), model = s.model; else, fn = fieldnames(s); model = s.(fn{1}); end
    model = prepareForCHRR(model);
    fprintf('model: %s (%d reactions, %d metabolites, %d coupling constraints)\n', modelFile, ...
        numel(model.rxns), numel(model.mets), size(model.C, 1));

    sol = optimizeCbModel(model);
    if sol.stat ~= 1
        warning('The %s model is not feasible (LP status %d). Skipped.', neuron, sol.stat);
        continue
    end

    % sampling: set 1 computes the rounded polytope P, sets 2..nSets reuse it
    sets = cell(1, nSets);
    elapsed = zeros(1, nSets);
    for k = 1:nSets
        opt = options;
        rng(seeds(k), 'twister');
        tic
        if k == 1
            opt.toRound = 1;
            [P, sets{k}] = sampleCbModel(model, [], 'CHRR', opt);
        else
            opt.toRound = 0;
            [~, sets{k}] = sampleCbModel(model, [], 'CHRR', opt, P);
        end
        elapsed(k) = toc;
        fprintf('set %d (seed %d): %d points in %.1f min\n', k, seeds(k), size(sets{k}, 2), elapsed(k) / 60);
    end

    % pilot set
    pilotSeed = masterSeed + 100 * i + 99;
    rng(pilotSeed, 'twister');
    pooled = [sets{:}];
    smallSamples = pooled(:, randperm(size(pooled, 2), min(nSmall, size(pooled, 2))));

    % save (v7.3: no 2 GB limit; read in Python with h5py)
    for k = 1:nSets
        varName = sprintf('samples%d', k);
        S = struct(varName, sets{k});
        save(setFiles{k}, '-struct', 'S', varName, '-v7.3');
    end
    save(fullfile(outDir, 'smallSamples.mat'), 'smallSamples', '-v7.3');

    % manifest with fingerprints
    manifest = struct();
    manifest.neuron = neuron;
    manifest.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    manifest.modelFile = strrep(modelFile, [projectDir filesep], '');
    manifest.modelSHA256 = sha256File(modelFile);
    manifest.nReactions = numel(model.rxns);
    manifest.sampler = 'CHRR (sampleCbModel); set 1 toRound = 1, other sets reuse its rounded polytope';
    manifest.nPointsReturned = options.nPointsReturned;
    manifest.nStepsPerPoint = options.nStepsPerPoint;
    manifest.rng = 'twister';
    manifest.masterSeed = masterSeed;
    manifest.seeds = seeds;
    manifest.pilotSeed = pilotSeed;
    manifest.environment = env;
    manifest.minutesPerSet = round(elapsed / 60, 2);
    manifest.samplesSHA256 = cellfun(@sha256Matrix, sets, 'UniformOutput', false);
    manifest.smallSamplesSHA256 = sha256Matrix(smallSamples);

    % compare with a previous manifest, if any
    if isfile(manifestFile)
        old = jsondecode(fileread(manifestFile));
        sameInputs = strcmp(old.modelSHA256, manifest.modelSHA256) && isequal(old.seeds(:)', seeds) ...
            && old.nPointsReturned == options.nPointsReturned && old.nStepsPerPoint == options.nStepsPerPoint;
        same = isequal(cellstr(old.samplesSHA256)', manifest.samplesSHA256);
        if same
            fprintf('REPRODUCED: all sample sets are bit-identical to the previous manifest.\n');
        elseif sameInputs
            fprintf(['Same model, seeds and settings, but the samples differ from the previous manifest ' ...
                '(different MATLAB, COBRA Toolbox or solver version?). Compare samplingSummary.csv.\n']);
        else
            fprintf('The previous manifest used a different model or settings; it is replaced.\n');
        end
        copyfile(manifestFile, fullfile(outDir, 'samplingManifest_previous.json'));
    end
    fid = fopen(manifestFile, 'w');
    fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
    fclose(fid);

    % summary statistics per reaction (no toolboxes needed)
    n = size(pooled, 2);
    sorted = sort(pooled, 2);
    q = @(p) sorted(:, max(1, min(n, round(p * n))));

    % quality report: steady state, bounds, agreement between sets
    res = max(abs(model.S * pooled), [], 1);
    viol = max(max(model.lb - pooled, 0) + max(pooled - model.ub, 0), [], 1);
    med = cellfun(@(x) median(x, 2), sets, 'UniformOutput', false);
    med = [med{:}];
    spread = max(med, [], 2) - min(med, [], 2);
    iqrAll = q(0.75) - q(0.25);
    agree = mean(spread <= 0.1 * max(iqrAll, 1e-6));
    % Gelman-Rubin potential scale reduction factor, the three sets as chains
    m = numel(sets); nk = size(sets{1}, 2);
    chainMean = cell2mat(cellfun(@(x) mean(x, 2), sets, 'UniformOutput', false));
    chainVar = cell2mat(cellfun(@(x) var(x, 0, 2), sets, 'UniformOutput', false));
    W = mean(chainVar, 2);
    Bn = var(chainMean, 0, 2);                     % B/n
    Rhat = sqrt(((nk - 1) / nk * W + Bn) ./ max(W, 1e-12));
    Rhat(W < 1e-12 & Bn < 1e-12) = 1;              % fixed reactions
    fid = fopen(fullfile(outDir, 'samplingReport.txt'), 'w');
    fprintf(fid, 'model: %s\nmodel SHA-256: %s\nreactions: %d\n', manifest.modelFile, manifest.modelSHA256, numel(model.rxns));
    fprintf(fid, 'sets: %d x %d points (CHRR, %d steps per point), seeds %s\n', nSets, options.nPointsReturned, ...
        options.nStepsPerPoint, mat2str(seeds));
    fprintf(fid, 'time: %.1f min\nmax |S v| over samples: %.3g\nmax bound violation: %.3g\n', ...
        sum(elapsed) / 60, max(res), max(viol));
    fprintf(fid, 'reactions whose median agrees between sets (within 10%% of the IQR): %.1f%%\n', 100 * agree);
    fprintf(fid, 'Gelman-Rubin R-hat: median %.3f, max %.3f; reactions with R-hat < 1.1: %.1f%%, < 1.01: %.1f%%\n', ...
        median(Rhat), max(Rhat), 100 * mean(Rhat < 1.1), 100 * mean(Rhat < 1.01));
    fclose(fid);
    summary = table(model.rxns, mean(pooled, 2), std(pooled, 0, 2), q(0.05), q(0.25), q(0.5), q(0.75), q(0.95), Rhat, ...
        'VariableNames', {'rxns', 'mean', 'sd', 'p05', 'p25', 'p50', 'p75', 'p95', 'Rhat'});
    writetable(summary, fullfile(outDir, 'samplingSummary.csv'));
    manifest.medianAgreementBetweenSets = round(agree, 4);
    manifest.rhatFractionBelow1p1 = round(mean(Rhat < 1.1), 4);
    fid = fopen(manifestFile, 'w');                % rewrite the manifest with the convergence diagnostics
    fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
    fclose(fid);
    fprintf('%s done: median agreement between sets %.1f%%, R-hat < 1.1 in %.1f%% of reactions\n', ...
        neuron, 100 * agree, 100 * mean(Rhat < 1.1));
    clear sets pooled sorted smallSamples P S
end

disp('Samples saved in results/matlab/sampling. Next: samplingAnalysisFluxQuery.ipynb')

%% Local functions
function model = prepareForCHRR(model)
% CHRR (chrrParseModel) reads coupling constraints C*v (dsense) d only in 'G' or 'L'
% form, and fails on rows that are 'L'. Write every inequality as 'G' and move
% equalities to S, so that all coupling constraints are kept during sampling.
if ~isfield(model, 'b'), model.b = zeros(numel(model.mets), 1); end
if ~isfield(model, 'C') || isempty(model.C)
    model.C = sparse(0, numel(model.rxns)); model.d = zeros(0, 1); model.dsense = char(zeros(0, 1));
    return
end
dsense = char(model.dsense(:));
isL = dsense == 'L';
model.C(isL, :) = -model.C(isL, :);
model.d(isL) = -model.d(isL);
dsense(isL) = 'G';
isE = dsense == 'E';
if any(isE)
    nE = nnz(isE);
    model.S = [model.S; model.C(isE, :)];
    model.b = [model.b; model.d(isE)];
    model.mets = [model.mets; strcat('couplingEquality_', cellstr(num2str((1:nE)', '%d')))];
    if isfield(model, 'csense'), model.csense = [model.csense(:); repmat('E', nE, 1)]; end
    model.C = model.C(~isE, :); model.d = model.d(~isE); dsense = dsense(~isE);
end
model.dsense = dsense;
end

function h = sha256Matrix(X)
% SHA-256 of the IEEE-754 bytes of a double matrix (column order), in blocks
md = java.security.MessageDigest.getInstance('SHA-256');
X = double(X);
block = max(1, floor(2e6 / max(1, size(X, 1))));
for c = 1:block:size(X, 2)
    chunk = X(:, c:min(size(X, 2), c + block - 1));
    md.update(typecast(chunk(:), 'int8'));
end
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end

function h = sha256File(file)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file, 'r');
while true
    bytes = fread(fid, 2^24, '*int8');
    if isempty(bytes), break, end
    md.update(bytes);
end
fclose(fid);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
