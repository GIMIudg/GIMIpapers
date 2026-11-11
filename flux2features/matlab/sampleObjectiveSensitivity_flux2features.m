%% sampleObjectiveSensitivity_flux2features
%
% Does the loss of oxidative reserve depend on holding NaKt at its maximum?
%
% Every chain of the main design (generateModelsAndSamples_flux2features.m)
% samples the face of the polytope on which the objective, the Na+/K+-ATPase
% (NaKt), runs at its maximum: sampleCbModel restricts the sampled box to the
% flux-variability ranges at 100% of the objective by default. This script
% samples two conditions of Original again with the objective held at no less
% than 95% and 90% of its maximum:
%
%   reference         RBC-GEM as distributed
%   G6PD_0_1_h2o2_1   G6PDH2 upper bound at 0.10 of its own, EX_h2o2_e at 1.00
%
% three CHRR chains of 10,000 points each, by the same procedure as the main
% design: same model file, same rounding, same solver fallback, same bound
% check. The fraction is imposed twice, so that it holds whatever the toolbox
% version does with it: as a lower bound on the objective reaction
% (alpha * max) and as options.optPercentage = 100 * alpha.
%
% Writes, under <projectDir>/results/raw/samples_objective/:
%
%   reference_alpha95_samples1.mat ... G6PD_0_1_h2o2_1_alpha90_samples3.mat
%       variable 'samples', reactions down the rows, one point per column
%   sampling_log.txt
%       solver of the rounding, seed, objective range and bound check of
%       every chain
%
% The notebook reads these files in stage 08c (run_08c_objective_sensitivity)
% and the manuscript reports them only once they exist.
%
% Run time: about the same per chain as the main design (12 chains in all).
% Requirements: COBRA Toolbox v3, Gurobi and/or MOSEK.

clear
projectDir = '';                      % e.g. '/path/to/MCBCI_terapiasPersonalizadas'

initCobraToolbox(false)
changeCobraSolver('gurobi', 'all', 0);

if isempty(projectDir)
    projectDir = uigetdir(pwd, 'Select the project folder (holds data/ and results/)');
end
outDir = fullfile(projectDir, 'results', 'raw', 'samples_objective');
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, 'sampling_log.txt');
logId = fopen(logFile, 'a');
cleanup = onCleanup(@() fclose(logId));
fprintf(logId, '\n==== %s\n', datestr(now));

modelFile = fullfile(projectDir, 'data', 'RBC-GEM.mat');
if ~exist(modelFile, 'file')         % the repository keeps a copy of the model
    modelFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'RBC-GEM.mat');
end
loaded = load(modelFile);
modelOr = loaded.RBC_GEM;

opt.alphas = [0.95, 0.90];            % objective held at no less than this fraction of its maximum
opt.g6Level = 0.1;                    % G6PDH2 upper bound, fraction of the model's own
opt.h2Level = 1;                      % EX_h2o2_e upper bound, fraction of the model's own
opt.nChains = 3;
opt.nPoints = 10000;
opt.roundingSolvers = {'gurobi', 'mosek', 'glpk'};
opt.baseSeed = 20260812;              % the main design's base seed; offsets 500+ are new

g6 = findRxnIDs(modelOr, 'G6PDH2');
h2 = findRxnIDs(modelOr, 'EX_h2o2_e');
restricted = changeRxnBounds(modelOr, 'G6PDH2', modelOr.ub(g6) * opt.g6Level, 'u');
restricted = changeRxnBounds(restricted, 'EX_h2o2_e', modelOr.ub(h2) * opt.h2Level, 'u');
conditions = { ...
    struct('stem', 'reference',       'model', modelOr,    'bound', []), ...
    struct('stem', 'G6PD_0_1_h2o2_1', 'model', restricted, ...
           'bound', struct('index', g6, 'value', restricted.ub(g6)))};

k = 0;
for a = 1:numel(opt.alphas)
    alpha = opt.alphas(a);
    for i = 1:numel(conditions)
        k = k + 1;
        condition = conditions{i};
        model = condition.model;
        objective = find(model.c);
        assert(numel(objective) == 1, 'expected a single objective reaction');
        best = optimizeCbModel(model, 'max');
        model.lb(objective) = alpha * best.f;
        stem = sprintf('%s_alpha%d', condition.stem, round(100 * alpha));
        fprintf(logId, '%s: objective %s, maximum %.4f, held at >= %.4f\n', stem, ...
            model.rxns{objective}, best.f, model.lb(objective));

        chains = sampleCondition(model, stem, condition.bound, alpha, opt, logId, 500 + k);
        for c = 1:opt.nChains
            samples = chains{c};
            values = samples(objective, :);
            fprintf(logId, '%s chain %d: objective from %.4f to %.4f (maximum %.4f)\n', ...
                stem, c, min(values), max(values), best.f);
            assert(min(values) >= alpha * best.f * (1 - 1e-6) - 1e-9, ...
                '%s chain %d: objective below the imposed fraction', stem, c);
            if min(values) > best.f * (1 - 1e-3)
                warning('%s chain %d: the objective never leaves its maximum; check optPercentage', ...
                    stem, c);
            end
            save(fullfile(outDir, sprintf('%s_samples%d', stem, c)), 'samples', '-v7')
        end
    end
end
fprintf('done; log in %s\n', logFile);


%% ------------------------------------------------------------------------
function chains = sampleCondition(model, stem, bound, alpha, opt, logId, seedOffset)
% Three CHRR chains from this condition's own rounded polytope, as in the main
% design, with the sampled box restricted at 100 * alpha % of the objective.
rng(opt.baseSeed + seedOffset, 'twister');
options.nPointsReturned = opt.nPoints;
options.toRound = 1;
options.optPercentage = 100 * alpha;
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
        top = max(chains{c}(bound.index, :));
        assert(top <= bound.value * (1 + 1e-6) + 1e-9, ...
            '%s chain %d: G6PDH2 reaches %.2f, above its bound %.2f', stem, c, top, bound.value);
        assert(top >= 0.95 * bound.value, ...
            '%s chain %d: G6PDH2 tops out at %.2f, well below its bound %.2f', stem, c, top, bound.value);
        fprintf(logId, '%s chain %d: rounding %s, seed %d, %d points, max G6PDH2 %.2f of %.2f\n', ...
            stem, c, used, opt.baseSeed + seedOffset, size(chains{c}, 2), top, bound.value);
    else
        fprintf(logId, '%s chain %d: rounding %s, seed %d, %d points\n', ...
            stem, c, used, opt.baseSeed + seedOffset, size(chains{c}, 2));
    end
end
end
