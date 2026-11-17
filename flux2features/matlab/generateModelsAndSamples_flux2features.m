%% generateModelsAndSamples_flux2features
%
% Every flux state of the flux2features study comes from this script. It
% builds the two reconstructions, imposes every condition and samples each
% one with three CHRR chains. Each chain is drawn from its own fixed seed,
% so running the script again with the same software returns the same
% chains, number for number.
%
%   Original   RBC-GEM as distributed (data/RBC-GEM.mat)
%   Pruned     its thermodynamically flux-consistent subset: the reactions
%              that can carry flux free of thermodynamically infeasible
%              cycles under the same bounds ("thermo" in file names)
%
%   design     G6PDH2 upper bound at 0.50, 0.28 and 0.10 of the model's own,
%              crossed with the EX_h2o2_e upper bound at 1.00, 0.75 and 0.50,
%              plus the unrestricted reference: 10 conditions per model,
%              with NaKt at its maximum (sampleCbModel's default)
%   objective  the reference and G6PDH2 0.10 (peroxide 1.00) of Original
%              again, with NaKt held at no less than 95% and 90% of its
%              maximum: 4 conditions
%
%   24 conditions x 3 chains x 10,000 points, about 7 minutes per chain on
%   a laptop (some 9 hours in all).
%
% Writes, under <projectDir>/results/raw/:
%
%   models/                    Pruned, and the constrained model of every condition
%   fluxes/                    the FBA optimum of every condition (v_or and
%                              v_thermoOr are the two references)
%   samples/                   the design chains, e.g. G6PD_0_1_h2o2_1_samples2.mat
%   samples_objective/         the objective chains, e.g. reference_alpha90_samples1.mat
%   sampling_log.txt           what happened, chain by chain
%   sampling_fingerprints.csv  seed, solver and SHA-256 of every chain
%
% The notebook (metabolic_flux_phenotyping.ipynb) reads the design chains in
% stage 01 and the objective chains in stage 08c.
%
% REPRODUCIBILITY
%
% Chain c of condition k is drawn from seed baseSeed + 100*k + c, and the
% rounding of condition k from seed baseSeed + 100*k (the table is printed
% when the script starts). The generator is reset to the chain's seed just
% before the chain is drawn, so a chain does not depend on anything that ran
% before it, and any one chain can be regenerated alone. CHRR in the COBRA
% Toolbox (chrrSampler, then genSamples of the Volume-and-Sampling package)
% draws its random numbers from MATLAB's global generator with rand, randn
% and randi and never reseeds it, so the seed fixes the chain.
%
% What a seed cannot fix is the arithmetic. Another version of MATLAB, of
% the COBRA Toolbox or of the solver can round the polytope or place the
% starting point slightly differently, and the chains then differ. The
% fingerprints tell you whether that happened: set mode = 'check' to sample
% a single chain (G6PDH2 0.10, peroxide 1.00, chain 1, about 7 minutes) and
% compare it with the fingerprint published with the code
% (matlab/sampling_fingerprints.csv). The fingerprint is the SHA-256 of the
% chain's numbers as saved (double precision, reactions down the rows).
%
% IF IT STOPS (MATLAB hangs, a solver error, a power cut), just run it again.
% Every chain whose file and fingerprint are both there is kept and skipped,
% conditions with all their chains are skipped without rounding, Pruned is
% loaded instead of rebuilt, and sampling continues with the first chain that
% is missing. A chain is saved under a temporary name and renamed only once
% it is complete, so a chain cut off while saving never passes for a good
% one; a leftover temporary file, or a chain without a fingerprint, is moved
% to results/raw/interrupted/ and drawn again from the same seed. Every line
% of the log carries the time, so the log shows where a stop happened.
% The rounded polytope of every condition is saved in models/ as soon as it is
% ready (<stem>_polytope.mat), so a run stopped while sampling starts again
% from the saved polytope instead of rounding for another 30-50 minutes.
% It is deleted once all the condition's chains are saved.
% If one condition keeps hanging while its polytope is rounded, stop MATLAB
% (Ctrl+C), put 'mosek' first in roundingSolvers and run again: the solver
% each condition used is written to the fingerprints.
%
% A sample file in a folder with no fingerprints at all comes from an older,
% unseeded run and stops the script, so that old and new chains never mix:
% move the old results/raw/ away first.
%
% TWO RULES THIS SCRIPT ENFORCES, because breaking either one silently
% produces chains that carry the wrong condition:
%
%   1. A rounded polytope is what CHRR samples, whatever model is passed
%      beside it. It is built afresh for every condition from that
%      condition's own model and never reused across conditions.
%   2. Every chain is checked before it is saved. A restricted chain must
%      stay at or under the G6PDH2 bound it is named after and reach at
%      least 95% of it. An objective chain must stay above its fraction of
%      the maximum and must actually leave the maximum. The peroxide bound is
%      never reached by any state (it is the design's inactive-bound
%      control) and is not checked.
%
% Rounding solves linear programs to find an interior starting point. At
% G6PDH2 0.28 in Pruned, Gurobi put that point on the boundary ("x0 not
% interior") and MOSEK found an interior one. Rounding is therefore tried
% with Gurobi, then MOSEK, then GLPK, the bounds always unchanged, and the
% solver that succeeded is written to the log and the fingerprints.
%
% RUNNING WITHOUT THE DESKTOP. If MATLAB itself stops responding (the desktop
% freezes after a waitbar fills), run the script from a terminal instead,
% where there is no desktop to freeze and less memory is used:
%
%   export FLUX2FEATURES_PROJECT=/path/to/MCBCI_terapiasPersonalizadas
%   /Applications/MATLAB_R2024b.app/bin/matlab -batch \
%       "run('/path/to/flux2features/matlab/generateModelsAndSamples_flux2features.m')"
%
% (use your MATLAB version in the path). The log and the fingerprints are the
% same either way, and the run resumes where the last one stopped.
%
% Requirements: COBRA Toolbox v3, Gurobi and/or MOSEK.
% Run from any folder; you are asked for the project folder (the one that
% holds data/ and results/), unless projectDir is set below.

clear
projectDir = '';       % e.g. '/path/to/MCBCI_terapiasPersonalizadas'; empty asks
mode = 'all';          % 'all' samples everything; 'check' samples one chain and compares it
baseSeed = 20260812;   % chain c of condition k uses baseSeed + 100*k + c

nChains = 3;
nPoints = 10000;
nStepsPerPoint = 200;                 % sampleCbModel's default, written out
roundingSolvers = {'gurobi', 'mosek', 'glpk'};
checkStem = 'G6PD_0_1_h2o2_1';         % the condition 'check' samples

% A chain that hangs. COBRA's sampler can get stuck for good when a numerical
% step leaves the polytope (a NaN in its rejection loop), and with a fixed
% seed it gets stuck at the same step every time. Give such a chain another
% seed here, by its file name without .mat, for example
%     seedOffset.G6PD_0_1_h2o2_1_alpha90_samples2 = 1000;
% The seed it then uses (the usual one plus the offset) is written to the log
% and the fingerprints, so the chain stays reproducible.
seedOffset = struct();

initCobraToolbox(false)
changeCobraSolver('gurobi', 'all', 0);

%% Folders, log and fingerprints

if isempty(projectDir)                 % from a terminal, the environment can name it
    projectDir = getenv('FLUX2FEATURES_PROJECT');
end
if isempty(projectDir)
    projectDir = uigetdir(pwd, 'Select the project folder (holds data/ and results/)');
end
codeDir = fileparts(mfilename('fullpath'));
rawDir = fullfile(projectDir, 'results', 'raw');
if strcmp(mode, 'check')              % a check writes nothing into the project
    rawDir = fullfile(tempdir, 'flux2features_check');
    if exist(rawDir, 'dir'), rmdir(rawDir, 's'); end
end
for sub = {'models', 'fluxes', 'samples', 'samples_objective', 'interrupted'}
    if ~exist(fullfile(rawDir, sub{1}), 'dir'), mkdir(fullfile(rawDir, sub{1})); end
end

logId = fopen(fullfile(rawDir, 'sampling_log.txt'), 'a');
cleanup = onCleanup(@() fclose(logId));
software = softwareVersions();
note(logId, '\n==== %s, mode %s', datestr(now), mode);
note(logId, 'MATLAB %s; COBRA Toolbox %s; Gurobi %s; MOSEK %s', software.matlab, ...
    software.cobra, software.gurobi, software.mosek);

fingerprintFile = fullfile(rawDir, 'sampling_fingerprints.csv');
written = readFingerprints(fingerprintFile);                         % this folder, to resume
published = readFingerprints(fullfile(codeDir, 'sampling_fingerprints.csv'));  % to compare

%% The two reconstructions

modelFile = fullfile(projectDir, 'data', 'RBC-GEM.mat');
if ~exist(modelFile, 'file')         % the repository keeps a copy of the model
    modelFile = fullfile(fileparts(codeDir), 'data', 'RBC-GEM.mat');
end
loaded = load(modelFile);
models.Original = loaded.RBC_GEM;
note(logId, 'Original: %s, %d reactions', modelFile, numel(models.Original.rxns));

prunedFile = fullfile(rawDir, 'models', 'thermoConsistModel.mat');
if strcmp(mode, 'all') && exist(prunedFile, 'file') && ~isempty(written.file)
    % resuming a seeded run: Pruned was built by this script on the first run
    loaded = load(prunedFile);
    models.Pruned = loaded.thermoConsistModel;
    note(logId, 'Pruned: loaded from %s (%d reactions)', prunedFile, numel(models.Pruned.rxns));
elseif strcmp(mode, 'all')
    feasTol = getCobraSolverParams('LP', 'feasTol');
    [~, ~, ~, ~, ~, ~, ~, stoichConsistModel] = findStoichConsistentSubset(models.Original, ...
        0, 0, [], feasTol * 10);
    [~, ~, ~, thermoConsistModel] = findThermoConsistentFluxSubset(stoichConsistModel);
    save(prunedFile, 'thermoConsistModel')
    models.Pruned = thermoConsistModel;
    note(logId, 'Pruned: %d of %d reactions kept', numel(models.Pruned.rxns), ...
        numel(models.Original.rxns));
end

if strcmp(mode, 'all')
    v_or = optimizeCbModel(models.Original, 'max').v;
    save(fullfile(rawDir, 'fluxes', 'v_or'), 'v_or')
    v_thermoOr = optimizeCbModel(models.Pruned, 'max').v;
    save(fullfile(rawDir, 'fluxes', 'v_thermoOr'), 'v_thermoOr')
end

%% The conditions, in a fixed order (the order fixes the seeds)

level = @(x) strrep(num2str(x), '.', '_');
conditions = struct('model', {}, 'g6', {}, 'h2', {}, 'alpha', {}, 'stem', {}, 'folder', {});
for name = {'Original', 'Pruned'}
    if strcmp(name{1}, 'Original'), prefix = ''; reference = 'modelOr';
    else, prefix = 'thermo_'; reference = 'thermoModelOr'; end
    conditions(end + 1) = condition(name{1}, 1, 1, 1, reference, 'samples');
    for g6 = [0.50, 0.28, 0.10]
        for h2 = [1, 0.75, 0.50]
            stem = sprintf('%sG6PD_%s_h2o2_%s', prefix, level(g6), level(h2));
            conditions(end + 1) = condition(name{1}, g6, h2, 1, stem, 'samples');
        end
    end
end
for alpha = [0.95, 0.90]
    tag = sprintf('_alpha%d', round(100 * alpha));
    conditions(end + 1) = condition('Original', 1, 1, alpha, ['reference' tag], 'samples_objective');
    conditions(end + 1) = condition('Original', 0.10, 1, alpha, ['G6PD_0_1_h2o2_1' tag], 'samples_objective');
end

note(logId, '\n  k  model     G6PDH2  H2O2  NaKt   seeds of chains 1-%d   file stem', nChains);
for k = 1:numel(conditions)
    c = conditions(k);
    note(logId, '%3d  %-8s  %5.2f  %5.2f  %4.0f%%  %d-%d     %s', k, c.model, c.g6, c.h2, ...
        100 * c.alpha, baseSeed + 100 * k + 1, baseSeed + 100 * k + nChains, c.stem);
end

%% Sampling

for k = 1:numel(conditions)
    c = conditions(k);
    chains = 1:nChains;
    if strcmp(mode, 'check')
        if ~strcmp(c.stem, checkStem), continue, end
        chains = 1;
    end

    % which chains of this condition are still missing
    missing = [];
    for chain = chains
        file = fullfile(c.folder, sprintf('%s_samples%d.mat', c.stem, chain));
        seed = chainSeed(baseSeed, k, chain, c.stem, seedOffset);
        setAside(fullfile(rawDir, [file(1:end - 4) '_partial.mat']), rawDir, logId);
        if exist(fullfile(rawDir, file), 'file')
            if any(strcmp(written.file, file) & written.seed == seed)
                continue                       % done in an earlier run
            end
            if isempty(written.file)
                error(['%s exists but this folder has no fingerprints: it comes from an ' ...
                    'older, unseeded run. Move the old results/raw/ away before sampling.'], file);
            end
            setAside(fullfile(rawDir, file), rawDir, logId);   % saved, but never fingerprinted
        end
        missing(end + 1) = chain; %#ok<AGROW>
    end
    if isempty(missing)
        note(logId, '%s: all %d chains already sampled, skipped', c.stem, numel(chains));
        continue
    end

    % the condition's own model, and its FBA optimum
    model = constrain(models.(c.model), c);
    fba = optimizeCbModel(model, 'max');
    save(fullfile(rawDir, 'models', [c.stem '_model']), 'model', '-v7')
    save(fullfile(rawDir, 'fluxes', [c.stem '_fba']), 'fba', '-v7')
    objective = find(model.c);
    note(logId, '\n%s (%s): NaKt maximum %.4f, held at >= %.4f', c.stem, c.model, ...
        fba.f, model.lb(objective));

    % the condition's own rounded polytope, from its own seed
    % The rounding is the slow step (30-50 minutes for some conditions) and does
    % not depend on the chain seeds, so it is saved as soon as it is done and
    % loaded on the next run instead of being repeated. The same rounding seed
    % gives the same polytope, so a loaded one changes no chain.
    options = struct('nStepsPerPoint', nStepsPerPoint, 'optPercentage', 100 * c.alpha);
    polytopeFile = fullfile(rawDir, 'models', [c.stem '_polytope.mat']);
    roundingSeed = baseSeed + 100 * k;
    [P, solver] = loadPolytope(polytopeFile, roundingSeed, logId, c.stem);
    if isempty(P)
        note(logId, '%s: rounding the polytope (chains to sample: %s); this can take 30-50 minutes', ...
            c.stem, mat2str(missing));
        started = tic;
        [P, solver] = roundPolytope(model, options, roundingSeed, roundingSolvers, logId, c.stem);
        save(polytopeFile, 'P', 'solver', 'roundingSeed', '-v7')
        note(logId, '%s: rounded with %s in %.0f minutes, saved to %s', c.stem, solver, ...
            toc(started) / 60, polytopeFile);
    end

    for chain = missing
        seed = chainSeed(baseSeed, k, chain, c.stem, seedOffset);
        file = fullfile(c.folder, sprintf('%s_samples%d.mat', c.stem, chain));
        label = sprintf('%s chain %d', c.stem, chain);

        note(logId, ['%s: sampling %d points with seed %d (about 6 minutes; if the waitbar ' ...
            'stops moving for 20 minutes the chain is stuck: see seedOffset)'], label, nPoints, seed);
        rng(seed, 'twister');
        options.nPointsReturned = nPoints;
        options.toRound = 0;
        [~, samples] = sampleCbModel(model, [], 'CHRR', options, P);

        checkChain(samples, model, c, fba.f, label, logId);
        % saved under a temporary name, renamed only once complete
        partial = fullfile(rawDir, [file(1:end - 4) '_partial.mat']);
        save(partial, 'samples', '-v7')
        movefile(partial, fullfile(rawDir, file));
        sha = sha256(samples);
        appendFingerprint(fingerprintFile, file, c, chain, seed, solver, ...
            baseSeed + 100 * k, size(samples, 2), sha, software);
        note(logId, '%s: seed %d, rounding %s, %d points, SHA-256 %s', label, seed, ...
            solver, size(samples, 2), sha);
        compareWithPublished(published, file, seed, sha, logId);
    end

    % every chain of the condition is saved: the polytope is no longer needed
    if exist(polytopeFile, 'file')
        delete(polytopeFile);
        note(logId, '%s: all chains saved, rounded polytope deleted', c.stem);
    end
end

if strcmp(mode, 'check')
    note(logId, '\nCheck done; nothing was written into the project (files in %s).', rawDir);
else
    note(logId, ['\nDone. Fingerprints in %s; copy that file to matlab/ to publish ' ...
        'them with the code.'], fingerprintFile);
end


%% ------------------------------------------------------------------------
function c = condition(model, g6, h2, alpha, stem, folder)
% One row of the design.
c = struct('model', model, 'g6', g6, 'h2', h2, 'alpha', alpha, 'stem', stem, 'folder', folder);
end


function model = constrain(base, c)
% The bounds of condition c, as fractions of the model's own upper bounds;
% below the maximal objective, NaKt is held at no less than alpha * maximum.
model = base;
g6 = findRxnIDs(base, 'G6PDH2');
h2 = findRxnIDs(base, 'EX_h2o2_e');
model = changeRxnBounds(model, 'G6PDH2', base.ub(g6) * c.g6, 'u');
model = changeRxnBounds(model, 'EX_h2o2_e', base.ub(h2) * c.h2, 'u');
if c.alpha < 1
    objective = find(model.c);
    assert(numel(objective) == 1, 'expected a single objective reaction');
    best = optimizeCbModel(model, 'max');
    model.lb(objective) = c.alpha * best.f;
end
end


function [P, solver] = loadPolytope(file, roundingSeed, logId, stem)
% The polytope rounded on an earlier run, if it was rounded from the same seed.
P = [];
solver = '';
if ~exist(file, 'file'), return, end
saved = load(file);
if saved.roundingSeed ~= roundingSeed
    note(logId, '%s: saved polytope has rounding seed %d, not %d; rounding again', ...
        stem, saved.roundingSeed, roundingSeed);
    return
end
P = saved.P;
solver = saved.solver;
note(logId, '%s: rounded polytope loaded from %s (rounded with %s)', stem, file, solver);
end


function [P, used] = roundPolytope(model, options, seed, solvers, logId, stem)
% Round the condition's polytope (maximum-volume ellipsoid), trying each
% solver in turn. This is exactly what chrrSampler does before it samples
% (chrrParseModel, then preprocess with the same options), stopped before
% any point is drawn: the polytope can then be saved at once, and no random
% walk runs until a chain's own seed is set. Rounding draws no random
% numbers; the generator is reset anyway, so nothing depends on which
% attempts failed.
P = [];
used = '';
for s = 1:numel(solvers)
    if ~changeCobraSolver(solvers{s}, 'LP', 0)
        continue
    end
    try
        rng(seed, 'twister');
        parsed = chrrParseModel(model);
        if (~isfield(parsed, 'A_eq') || isempty(parsed.A_eq)) && ...
                (~isfield(parsed, 'b_eq') || isempty(parsed.b_eq))
            parsed.A_eq = [];
            parsed.b_eq = [];
        end
        settings = struct('toRound', 1, 'bioModel', 1, 'fullDim', 0, 'model', model, ...
            'useFastFVA', false, 'optPercentage', options.optPercentage);
        P = preprocess(parsed, settings);
        used = solvers{s};
    catch err
        note(logId, '%s: rounding with %s failed (%s)', stem, solvers{s}, err.message);
    end
    changeCobraSolver('gurobi', 'LP', 0);
    if ~isempty(used), break, end
end
if isempty(used)
    error('%s: the rounding failed with every solver; see the log', stem);
end
end


function checkChain(samples, model, c, best, label, logId)
% Rule 2: the chain carries the condition it is named after.
if c.g6 < 1
    g6 = findRxnIDs(model, 'G6PDH2');
    bound = model.ub(g6);
    top = max(samples(g6, :));
    assert(top <= bound * (1 + 1e-6) + 1e-9, ...
        '%s: G6PDH2 reaches %.2f, above its bound %.2f', label, top, bound);
    assert(top >= 0.95 * bound, ...
        '%s: G6PDH2 tops out at %.2f, well below its bound %.2f (wrong polytope?)', ...
        label, top, bound);
    note(logId, '%s: max G6PDH2 %.2f of %.2f', label, top, bound);
end
if c.alpha < 1
    values = samples(find(model.c), :); %#ok<FNDSB>
    assert(min(values) >= c.alpha * best * (1 - 1e-6) - 1e-9, ...
        '%s: NaKt falls below %.0f%% of its maximum', label, 100 * c.alpha);
    if min(values) > best * (1 - 1e-3)
        warning('%s: NaKt never leaves its maximum; check optPercentage', label);
    end
    note(logId, '%s: NaKt from %.4f to %.4f (maximum %.4f)', label, min(values), ...
        max(values), best);
end
end


function hex = sha256(samples)
% SHA-256 of the chain as saved: double precision, column by column. Fed to
% Java in pieces so that a 260 MB chain never sits in the Java heap at once.
digest = java.security.MessageDigest.getInstance('SHA-256');
values = samples(:);
step = 1e6;
for first = 1:step:numel(values)
    piece = values(first:min(first + step - 1, numel(values)));
    digest.update(typecast(piece, 'int8'));
end
hex = sprintf('%02x', typecast(digest.digest(), 'uint8'));
end


function fingerprints = readFingerprints(file)
% The fingerprints already written (empty if the file does not exist).
fingerprints = struct('file', {{}}, 'seed', [], 'sha256', {{}});
if ~exist(file, 'file'), return, end
rows = readtable(file, 'TextType', 'char', 'Delimiter', ',');
if isempty(rows), return, end
fingerprints.file = rows.file;
fingerprints.seed = rows.seed;
fingerprints.sha256 = rows.sha256;
end


function appendFingerprint(file, name, c, chain, seed, solver, roundingSeed, points, sha, software)
% One line per chain, written as soon as the chain is saved.
isNew = ~exist(file, 'file');
id = fopen(file, 'a');
if isNew
    fprintf(id, ['file,model,g6pd,h2o2,objective_fraction,chain,seed,rounding_seed,' ...
        'rounding_solver,points,sha256,matlab,cobra,gurobi,mosek,written\n']);
end
fprintf(id, '%s,%s,%.2f,%.2f,%.2f,%d,%d,%d,%s,%d,%s,%s,%s,%s,%s,%s\n', name, c.model, ...
    c.g6, c.h2, c.alpha, chain, seed, roundingSeed, solver, points, sha, software.matlab, ...
    software.cobra, software.gurobi, software.mosek, datestr(now, 'yyyy-mm-dd HH:MM'));
fclose(id);
end


function compareWithPublished(published, file, seed, sha, logId)
% Is this the chain published with the code?
match = strcmp(published.file, file) & published.seed == seed;
if ~any(match)
    return
end
if strcmp(published.sha256{find(match, 1)}, sha)
    note(logId, '  identical to the published chain');
else
    note(logId, ['  DIFFERS from the published chain: compare your MATLAB, COBRA ' ...
        'Toolbox and solver versions with those in matlab/sampling_fingerprints.csv']);
end
end


function software = softwareVersions()
% What produced the chains, for the log and the fingerprints.
software.matlab = strrep(version, ',', ';');
cobraRoot = fileparts(which('initCobraToolbox'));
[status, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', cobraRoot));
if status == 0
    software.cobra = strtrim(commit);
else
    software.cobra = 'unknown';
end
software.gurobi = solverPath('gurobi');
software.mosek = solverPath('mosekopt');
end


function where = solverPath(name)
% The solver's installation path, which names its version.
where = which(name);
if isempty(where), where = 'not installed'; end
where = strrep(where, ',', ';');
end


function seed = chainSeed(baseSeed, k, chain, stem, seedOffset)
% Chain c of condition k: baseSeed + 100*k + c, plus any offset set for it above.
seed = baseSeed + 100 * k + chain;
name = sprintf('%s_samples%d', stem, chain);
if isfield(seedOffset, name)
    seed = seed + seedOffset.(name);
end
end


function setAside(file, rawDir, logId)
% Move an incomplete chain out of the way, keeping it for inspection.
if ~exist(file, 'file'), return, end
[~, name, ext] = fileparts(file);
target = fullfile(rawDir, 'interrupted', [name '_' datestr(now, 'yyyymmdd_HHMMSS') ext]);
movefile(file, target);
note(logId, '%s: incomplete (interrupted run), moved to %s; it is drawn again', ...
    [name ext], target);
end


function note(logId, varargin)
% Print a line, with the time, and write it to the log.
text = sprintf(varargin{:});
stamp = datestr(now, 'HH:MM');
if startsWith(text, newline)
    text = [newline stamp '  ' text(2:end)];
else
    text = [stamp '  ' text];
end
fprintf('%s\n', text);
fprintf(logId, '%s\n', text);
end
