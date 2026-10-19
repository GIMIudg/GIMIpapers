%% generateModelsAndSamples_flux2features
%
% Builds the two reconstructions of the flux2features study and samples every
% condition of its design with CHRR:
%
%   Original  RBC-GEM as distributed (data/RBC-GEM.mat)
%   Thermo    its stoichiometrically consistent subset, reduced to the
%             reactions that can carry flux free of thermodynamically
%             infeasible cycles under the same bounds
%
%   design    G6PDH2 upper bound at 0.50, 0.28 and 0.10 of the model's own
%             upper bound, crossed with EX_h2o2_e upper bound at 1.00, 0.75
%             and 0.50 of its own, plus the unconstrained reference
%             (10 conditions per model, 3 chains per condition)
%
% Writes, under <projectDir>/results/raw/:
%
%   models/   the constrained model of every condition
%   fluxes/   the FBA optimum of every condition
%   samples/  three CHRR chains per condition, 10,000 points each
%   sampling_log.txt   solver used for the rounding, seed and bound check
%                      of every chain
%
% The file and variable names are the ones the notebook
% (metabolic_flux_phenotyping.ipynb, stage 00) reads.
%
% Two rules this script enforces, because breaking either one silently
% produces chains that carry the wrong condition:
%
%   1. A rounded polytope P handed to sampleCbModel is what CHRR samples,
%      whatever model is passed beside it. P is therefore built afresh for
%      every condition from that condition's model and never reused across
%      conditions.
%   2. Every perturbed chain is checked against the G6PDH2 bound it is named
%      after before it is saved: it must stay at or under the bound and reach
%      at least 95% of it (the healthy distribution of G6PDH2 reaches ~900 in
%      both reconstructions, so every imposed bound binds). The peroxide
%      bound restricts secretion, which no sampled state approaches; it is
%      the design's negative control and is not checked here.
%
% Rounding the polytope (maximum-volume ellipsoid) solves linear programs to
% find an interior starting point. At G6PDH2 = 0.28 in Thermo that point fell
% on the boundary with Gurobi ("x0 not interior"); MOSEK finds an interior
% one. The rounding of every condition is therefore tried with Gurobi, then
% MOSEK, then GLPK, with the bound always imposed unchanged, and the solver
% that succeeded is logged.
%
% Requirements: COBRA Toolbox v3, Gurobi and/or MOSEK.
% Run from any folder; you are asked for the project folder (the one that
% holds data/ and results/), unless projectDir is set below.

clear
projectDir = '';                      % e.g. '/path/to/MCBCI_terapiasPersonalizadas'

initCobraToolbox(false)
changeCobraSolver('gurobi', 'all', 0);

if isempty(projectDir)
    projectDir = uigetdir(pwd, 'Select the project folder (holds data/ and results/)');
end
rawDir = fullfile(projectDir, 'results', 'raw');
for sub = {'fluxes', 'models', 'samples'}
    target = fullfile(rawDir, sub{1});
    if ~exist(target, 'dir'), mkdir(target); end
end
logFile = fullfile(rawDir, 'sampling_log.txt');
logId = fopen(logFile, 'a');
cleanup = onCleanup(@() fclose(logId));
fprintf(logId, '\n==== %s\n', datestr(now));

modelFile = fullfile(projectDir, 'data', 'RBC-GEM.mat');
if ~exist(modelFile, 'file')         % the repository keeps a copy of the model
    modelFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'RBC-GEM.mat');
end
loaded = load(modelFile);
modelOr = loaded.RBC_GEM;

opt.g6Levels = [.5, .28, .1];         % fractions of the model's own G6PDH2 upper bound
opt.h2Levels = [1, .75, .5];          % fractions of the model's own EX_h2o2_e upper bound
opt.nChains = 3;
opt.nPoints = 10000;
opt.roundingSolvers = {'gurobi', 'mosek', 'glpk'};
opt.baseSeed = 20260812;

%% Original: the unconstrained reference
v_or = optimizeCbModel(modelOr, 'max').v;
save(fullfile(rawDir, 'fluxes', 'v_or'), 'v_or')
chains = sampleCondition(modelOr, 'modelOr', [], opt, logId, 1);
for c = 1:opt.nChains
    name = sprintf('modelOr_samples%d', c);
    S.(name) = chains{c};
    save(fullfile(rawDir, 'samples', name), '-struct', 'S', name)
    S = struct();
end

%% Original: the perturbed conditions
sampleDesign(modelOr, '', '_model', '_sol', opt, rawDir, logId, 100);

%% Thermo: the consistent subset
feasTol = getCobraSolverParams('LP', 'feasTol');
[~, ~, ~, ~, ~, ~, ~, stoichConsistModel] = findStoichConsistentSubset(modelOr, ...
    0, 0, [], feasTol * 10);
[~, ~, ~, thermoConsistModel] = findThermoConsistentFluxSubset(stoichConsistModel);
save(fullfile(rawDir, 'models', 'thermoConsistModel'), 'thermoConsistModel')
fprintf(logId, 'Thermo: %d of %d reactions kept\n', numel(thermoConsistModel.rxns), ...
    numel(modelOr.rxns));

v_thermoOr = optimizeCbModel(thermoConsistModel, 'max').v;
save(fullfile(rawDir, 'fluxes', 'v_thermoOr'), 'v_thermoOr')
chains = sampleCondition(thermoConsistModel, 'thermoModelOr', [], opt, logId, 200);
for c = 1:opt.nChains
    name = sprintf('thermoModelOr_samples%d', c);
    S.(name) = chains{c};
    save(fullfile(rawDir, 'samples', name), '-struct', 'S', name)
    S = struct();
end

%% Thermo: the perturbed conditions
sampleDesign(thermoConsistModel, 'thermo_', '_thermoModel', '_thermoSol', opt, rawDir, ...
    logId, 300);

fprintf('done; log in %s\n', logFile);


%% ------------------------------------------------------------------------
function sampleDesign(base, prefix, modelSuffix, solSuffix, opt, rawDir, logId, seedBlock)
% Every condition of the 3 x 3 design: constrain, solve, sample, check, save.
g6 = findRxnIDs(base, 'G6PDH2');
h2 = findRxnIDs(base, 'EX_h2o2_e');
k = 0;
for i = 1:numel(opt.g6Levels)
    for j = 1:numel(opt.h2Levels)
        k = k + 1;
        model = changeRxnBounds(base, 'G6PDH2', base.ub(g6) * opt.g6Levels(i), 'u');
        model = changeRxnBounds(model, 'EX_h2o2_e', base.ub(h2) * opt.h2Levels(j), 'u');
        stem = sprintf('%sG6PD_%s_h2o2_%s', prefix, ...
            strrep(num2str(opt.g6Levels(i)), '.', '_'), strrep(num2str(opt.h2Levels(j)), '.', '_'));

        fbaSol = optimizeCbModel(model, 'max');
        save(fullfile(rawDir, 'models', [stem modelSuffix]), 'model', '-v7')
        save(fullfile(rawDir, 'fluxes', [stem solSuffix]), 'fbaSol', '-v7')

        bound = struct('index', g6, 'value', model.ub(g6));
        chains = sampleCondition(model, stem, bound, opt, logId, seedBlock + k);
        for c = 1:opt.nChains
            samples = chains{c};
            save(fullfile(rawDir, 'samples', sprintf('%s_samples%d', stem, c)), 'samples', '-v7')
        end
    end
end
end


function chains = sampleCondition(model, stem, bound, opt, logId, seedOffset)
% Three CHRR chains from this condition's own rounded polytope.
%
% The rounding is attempted with each solver in turn until it succeeds; the
% bound is never altered. The polytope P is local to this call, so it cannot
% leak into another condition. With a non-empty bound, every chain is checked
% against it before it is returned.
rng(opt.baseSeed + seedOffset, 'twister');
options.nPointsReturned = opt.nPoints;
options.toRound = 1;
P = [];
used = '';
for s = 1:numel(opt.roundingSolvers)
    solver = opt.roundingSolvers{s};
    if ~changeCobraSolver(solver, 'LP', 0)
        continue
    end
    try
        [P, first] = sampleCbModel(model, [], 'CHRR', options);
        used = solver;
    catch err
        fprintf(logId, '%s: rounding with %s failed (%s)\n', stem, solver, err.message);
        fprintf('%s: rounding with %s failed (%s)\n', stem, solver, err.message);
    end
    changeCobraSolver('gurobi', 'LP', 0);
    if ~isempty(used), break, end
end
if isempty(used)
    error('%s: the rounding failed with every solver; see the log', stem);
end

chains = cell(1, opt.nChains);
chains{1} = first;
options.toRound = 0;
for c = 2:opt.nChains
    [~, chains{c}] = sampleCbModel(model, [], 'CHRR', options, P);
end

for c = 1:opt.nChains
    if ~isempty(bound)
        top = checkBound(chains{c}, bound.index, bound.value, sprintf('%s chain %d', stem, c));
        fprintf(logId, '%s chain %d: rounding %s, seed %d, %d points, max G6PDH2 %.2f of %.2f\n', ...
            stem, c, used, opt.baseSeed + seedOffset, size(chains{c}, 2), top, bound.value);
    else
        fprintf(logId, '%s chain %d: rounding %s, seed %d, %d points\n', ...
            stem, c, used, opt.baseSeed + seedOffset, size(chains{c}, 2));
    end
end
end


function top = checkBound(samples, index, bound, label)
% The chain must stay at or under the imposed bound and must reach it.
top = max(samples(index, :));
assert(top <= bound * (1 + 1e-6) + 1e-9, ...
    '%s: G6PDH2 reaches %.2f, above its bound %.2f', label, top, bound);
assert(top >= 0.95 * bound, ...
    '%s: G6PDH2 tops out at %.2f, well below its bound %.2f (wrong polytope?)', label, top, bound);
end
