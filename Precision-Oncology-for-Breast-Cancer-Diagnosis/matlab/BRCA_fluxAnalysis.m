function BRCA_fluxAnalysis()
% BRCA_fluxAnalysis  Fluxes, metabolic metrics and Pareto surfaces for every model built by modelGen.m.
%
% Replaces Fluxes_calculation.m, Norms_1_2_fluxes-nurvo.mlx, prepareModel.m and countCarbons.m
% (same solutions, same metrics, same column names) and writes everything that
% BRCA_Metabolic_Flux_Pipeline_02.ipynb reads.
%
% HOW TO RUN: edit the CONFIGURATION section below (or leave dataFolder empty to pick it in a
% dialog, as in modelGen.m) and press Run. COBRA Toolbox v3 and Gurobi (LP and QP) are required.
% The run can be interrupted: models already processed are kept in matlabResults/perModel and
% skipped when the script is run again.
%
% INPUT   <dataFolder>/matlabResults/models/<sample>_Xomics/Model.mat   (output of modelGen.m)
%         (.mat files directly in the models folder are also accepted)
% OUTPUT  <dataFolder>/matlabResults/
%   FeatureMatrix_TumorPhenotype_All.csv  one row per model: Model, Flux_<rxn>_<sol>, metrics <root>_<sol>,
%                                         SA_<subsystem>_<sol>, Oncomet_<met>_<sol>, nRxns, nMets, nGenes,
%                                         MaxBiomass_FBA, QC_pass
%   ModelQC.csv                           PASS / WARN / FAIL per model and why
%   ParetoSurface_100sol.csv              Model, Solution, Biomass, ATP, CU, EA (10 x 10 epsilon-constraint grid)
%   FeatureDictionary.csv                 column -> type, solution, Recon3D reaction / subsystem
%   MaxBiomass_perModel.csv               as in Fluxes_calculation.m
%   perModel/  run_log.txt  run_config.mat
%
% FLUX SOLUTIONS (as in Fluxes_calculation.m; biomass >= 0.99 * FBA optimum in solutions 2-5)
%   FBA   max biomass_ob
%   pFBA  min sum |v|
%   L1w   min sum w.*|v|
%   L2    min 1/2 sum v.^2
%   L2w   min 1/2 sum w.*v.^2
%
% METRICS: definitions of Fluxes_calculation.m, translated to Recon3D IDs (EX_glc__D_e -> EX_glc_D[e],
% atp_c -> atp[c], ATPS4m -> ATPS4minew, ICDHym -> ICDHyrm, ACACT7m -> r0653). Changes, all bug fixes:
%   - countCarbons no longer counts Cl, Co, Ca, Cu... as carbon.
%   - ratios return NaN instead of dividing by eps (AnabolismScore, GlnDependence) when the denominator is 0.
%   - a failed solution gives NaN (the old scripts wrote zeros, which looked like real "no flux" models).
%   - prepareModel.m is not needed: the new models already have biomass_ob without product and ATPM.
%     The old reversible ATP_HYDROLYSIS (-1000..1000) is not added because it creates ATP from nothing.
%   - expression weights: see cfg.weightScheme.
%
% MODEL QC (never stops the run; a FAIL model gets NaN features and QC_pass = 0)
%   FAIL  no Model.mat, biomass_ob missing, FBA infeasible or biomass = 0, XomicsToModel relaxed the
%         constraints (relaxationUsed = 1).
%   WARN  metric reactions absent from the model (metric computed with the reactions present, as before),
%         biomass at its upper bound, ATP produced with every uptake closed (energy-generating cycle;
%         FAIL if cfg.failOnEGC = true), no usable expression weights.
%
% PARETO: replaces Pareto_multiobjective_1226_models.m with a fixed 10 x 10 grid (biomass level x ATP level,
% minimum carbon uptake at each point) instead of the Latin-hypercube CU-EA front.

%% ======================================= CONFIGURATION =======================================
cfg = struct();
cfg.dataFolder   = '';            % folder that contains matlabResults/ (the one modelGen.m asks for); '' = dialog
cfg.useParallel  = false;         % true: parfor over models (Parallel Computing Toolbox)
cfg.overwrite    = false;         % true: recompute models already in perModel/
cfg.doPareto     = true;          % ParetoSurface_100sol.csv
cfg.lpSolver     = 'gurobi';
cfg.qpSolver     = 'gurobi';
cfg.fracOpt      = 0.99;          % biomass kept in pFBA / L1w / L2 / L2w (as in Fluxes_calculation.m)
cfg.failOnRelaxation = true;      % relaxed models were part of the old "Divergent" artefact
cfg.failOnEGC    = false;         % ATP made with every uptake closed: FAIL (true) or WARN (false)
% Expression weights for L1w / L2w (model.expressionRxns, linear FPKM):
%   'lowExpressionCost'  w = 1 - 0.99*e, e = log(1+FPKM) scaled to [0,1]: weakly expressed reactions cost
%                        more, so flux is routed through expressed ones (default).
%   'original'           w = expr + |min(expr)| + 1e-6 as in Fluxes_calculation.m: the cost GROWS with
%                        expression, i.e. it penalises the most expressed reactions. Only to reproduce
%                        the submitted manuscript.
cfg.weightScheme = 'lowExpressionCost';
%% ==============================================================================================

cfg = defaultConfig(cfg);
if isempty(cfg.dataFolder)
    cfg.dataFolder = uigetdir(pwd, 'Select the folder that contains matlabResults');
    if isequal(cfg.dataFolder, 0), error('No folder selected.'); end
end
cfg.modelDir = fullfile(cfg.dataFolder, 'matlabResults', 'models');
cfg.outDir   = fullfile(cfg.dataFolder, 'matlabResults');
assert(isfolder(cfg.modelDir), 'Model folder not found: %s', cfg.modelDir);
pmDir = fullfile(cfg.outDir, 'perModel');
if ~isfolder(pmDir), mkdir(pmDir); end
logFile = fullfile(cfg.outDir, 'run_log.txt');
logmsg(logFile, '=== BRCA_fluxAnalysis started ===');
save(fullfile(cfg.outDir, 'run_config.mat'), 'cfg');

if exist('optimizeCbModel', 'file') ~= 2, initCobraToolbox(false); end
setSolvers(cfg);

% ------------------------------------------------------------------ list models
jobs = findModels(cfg.modelDir);
assert(~isempty(jobs), 'No models found in %s', cfg.modelDir);
logmsg(logFile, sprintf('%d model folders/files found (%d without Model.mat)', numel(jobs), sum(cellfun(@isempty, {jobs.path}))));

todo = false(numel(jobs), 1);
for i = 1:numel(jobs)
    todo(i) = cfg.overwrite || ~isfile(fullfile(pmDir, [jobs(i).id '.mat']));
end
idx = find(todo);
logmsg(logFile, sprintf('%d to process, %d already done', numel(idx), sum(~todo)));

% ------------------------------------------------------------------ per model
if cfg.useParallel && ~isempty(idx)
    if isempty(gcp('nocreate')), parpool; end
    parfevalOnAll(@setSolvers, 0, cfg);
    sub = jobs(idx);
    parfor k = 1:numel(sub)
        processOne(sub(k), cfg, pmDir);
    end
else
    for k = 1:numel(idx)
        processOne(jobs(idx(k)), cfg, pmDir);
        if mod(k, 10) == 0 || k == numel(idx)
            logmsg(logFile, sprintf('  %d / %d', k, numel(idx)));
        end
    end
end

% ------------------------------------------------------------------ assemble
S = assembleOutputs(jobs, pmDir, cfg, logFile);
logmsg(logFile, sprintf('Done: %d models | PASS %d | WARN %d | FAIL %d', S.nModels, S.nPass, S.nWarn, S.nFail));
end

% =============================================================================================
function jobs = findModels(modelDir)
% one job per sample: <modelDir>/<name>/Model.mat (modelGen.m) or <modelDir>/<name>.mat
jobs = struct('id', {}, 'path', {});
d = dir(modelDir);
for i = 1:numel(d)
    nm = d(i).name;
    if startsWith(nm, '.') || strcmpi(nm, 'Icon') || startsWith(nm, 'Icon'), continue, end
    if d(i).isdir
        f = fullfile(modelDir, nm, 'Model.mat');
        if ~isfile(f)
            m = dir(fullfile(modelDir, nm, '*.mat'));
            if ~isempty(m), f = fullfile(m(1).folder, m(1).name); else, f = ''; end
        end
        jobs(end+1) = struct('id', sampleId(nm), 'path', f); %#ok<AGROW>
    elseif endsWith(nm, '.mat')
        jobs(end+1) = struct('id', sampleId(nm), 'path', fullfile(modelDir, nm)); %#ok<AGROW>
    end
end
[~, u] = unique({jobs.id}, 'stable');
jobs = jobs(u);
end

function id = sampleId(name)
% 'TCGA_3C_AAAU_01A_Xomics' or 'TCGA-3C-AAAU-01A_Xomics_specificModel.mat' -> 'TCGA-3C-AAAU-01A'
t = regexp(name, 'TCGA[-_.][A-Za-z0-9]{2}[-_.][A-Za-z0-9]{4}[-_.]\d{2}[A-Z]?', 'match', 'once');
if isempty(t), [~, id] = fileparts(name); else, id = regexprep(t, '[_.]', '-'); end
end

% =============================================================================================
function processOne(job, cfg, pmDir)
R = struct('id', job.id, 'file', job.path, 'messages', {{}});
t0 = tic;
try
    if isempty(job.path), error('no Model.mat (XomicsToModel did not finish)'); end
    model = loadModel(job.path);
    [model, V] = validateModel(model, cfg);
    R.V = V;
    R.rxns = model.rxns;
    R.size = [numel(model.rxns), numel(model.mets), numel(model.genes)];
    R.sol = emptySol(cfg, numel(model.rxns));
    R.metrics = []; R.SA = []; R.pareto = [];
    if ~strcmp(V.status, 'FAIL')
        R.sol = computeSolutions(model, V, cfg);
        [R.metrics, R.SA] = computeMetrics(model, R.sol, V, cfg);
        if cfg.doPareto
            try
                R.pareto = paretoSurface(model, V, cfg);
            catch ME
                R.messages{end+1} = ['Pareto failed: ' ME.message];
            end
        end
    end
catch ME
    R.V = struct('status', 'FAIL', 'reasons', {{ME.message}}, 'warnings', {{}});
    R.messages{end+1} = getReport(ME, 'basic');
end
R.seconds = toc(t0);
saveResult(fullfile(pmDir, [job.id '.mat']), R);
end

function sol = emptySol(cfg, n)
sol = struct('name', cfg.solNames, 'v', nan(n, 1), 'stat', NaN, 'msg', '');
end

function model = loadModel(path)
S = load(path);
fn = fieldnames(S);
model = [];
for k = {'Model', 'specificModel', 'model'}          % names used by modelGen.m and older scripts
    if isfield(S, k{1}) && isstruct(S.(k{1})), model = S.(k{1}); break, end
end
for i = 1:numel(fn)
    if isempty(model) && isstruct(S.(fn{i})) && isfield(S.(fn{i}), 'S') && isfield(S.(fn{i}), 'rxns')
        model = S.(fn{i});
    end
end
if isempty(model), error('no COBRA model in %s', path); end
n = numel(model.rxns); m = numel(model.mets);
model.rxns = model.rxns(:); model.mets = model.mets(:);
model.lb = full(model.lb(:)); model.ub = full(model.ub(:));
if ~isfield(model, 'b') || numel(model.b) ~= m, model.b = zeros(m, 1); end
if ~isfield(model, 'csense') || numel(model.csense) ~= m, model.csense = repmat('E', m, 1); end
model.b = full(model.b(:)); model.csense = model.csense(:);
model.c = zeros(n, 1); model.osenseStr = 'max';
if ~isfield(model, 'genes'), model.genes = {}; end
if ~isfield(model, 'metFormulas'), model.metFormulas = repmat({''}, m, 1); end
if ~isfield(model, 'subSystems'), model.subSystems = repmat({''}, n, 1); end
model.subSystems = flatSubsystems(model.subSystems);
end

function ss = flatSubsystems(ss)
% getUniqueSubsystems of Fluxes_calculation.m: nested cells -> first entry, empty -> ''
ss = ss(:);
for i = 1:numel(ss)
    if iscell(ss{i}), if isempty(ss{i}), ss{i} = ''; else, ss{i} = ss{i}{1}; end, end
    if isempty(ss{i}) || ~ischar(ss{i}), ss{i} = ''; end
end
end

% =============================================================================================
function [model, V] = validateModel(model, cfg)
L = cfg.lists;
V = struct('status', 'PASS', 'reasons', {{}}, 'warnings', {{}});
V.bioIdx = find(strcmp(model.rxns, cfg.biomassRxn));
V.atpIdx = find(strcmp(model.rxns, cfg.atpDemandRxn));
if isempty(V.atpIdx) && ismember('ATP_HYDROLYSIS', model.rxns)      % legacy models (prepareModel.m)
    V.atpIdx = find(strcmp(model.rxns, 'ATP_HYDROLYSIS'));
    model.lb(V.atpIdx) = max(model.lb(V.atpIdx), 0);                % irreversible: no ATP synthesis
    V.warnings{end+1} = 'ATPM absent: legacy ATP_HYDROLYSIS used (made irreversible)';
end
if isempty(V.bioIdx), V.reasons{end+1} = [cfg.biomassRxn ' missing']; end

allR = unique([L.warburgGly, L.warburgOx, L.cu, L.ea, L.redoxNAD, L.redoxNADH, L.mfi, L.anab, ...
    L.catab, L.oncomet, L.nadph, L.lipSat, L.lipUnsat, L.lipPL]);
V.missingMetricRxns = setdiff(allR, model.rxns);
if ~isempty(V.missingMetricRxns)
    V.warnings{end+1} = ['metric reactions absent: ' strjoin(V.missingMetricRxns, ';')];
end

V.relaxationUsed = NaN;
if isfield(model, 'relaxationUsed') && (isnumeric(model.relaxationUsed) || islogical(model.relaxationUsed))
    V.relaxationUsed = double(any(model.relaxationUsed(:)));
    if V.relaxationUsed == 1
        msg = 'XomicsToModel relaxed the constraints (relaxationUsed = 1)';
        if cfg.failOnRelaxation, V.reasons{end+1} = msg; else, V.warnings{end+1} = msg; end
    end
end

% FBA
V.maxBiomass = NaN; V.fbaStat = NaN; V.vFBA = []; V.biomassLB = NaN;
if ~isempty(V.bioIdx)
    V.biomassLB = model.lb(V.bioIdx);
    m1 = model; m1.c(V.bioIdx) = 1;
    s = optimizeCbModel(m1, 'max');
    V.fbaStat = s.stat;
    if s.stat == 1 && s.f > cfg.tol
        V.maxBiomass = s.f; V.vFBA = s.x;
        if s.f >= 0.999 * model.ub(V.bioIdx)
            V.warnings{end+1} = 'FBA biomass at its upper bound';
        end
    elseif s.stat == 1
        V.reasons{end+1} = 'FBA biomass = 0';
    else
        V.reasons{end+1} = sprintf('FBA not optimal (stat = %d)', s.stat);
    end
end

% energy-generating cycles: close every uptake, release forced fluxes, maximise ATP demand
V.egcATP = NaN;
if ~isempty(V.atpIdx)
    m2 = model;
    one = find(full(sum(m2.S ~= 0, 1)) == 1);                      % exchange / demand / sink reactions
    for j = one
        if full(m2.S(m2.S(:, j) ~= 0, j)) < 0, m2.lb(j) = max(m2.lb(j), 0); else, m2.ub(j) = min(m2.ub(j), 0); end
    end
    m2.lb(m2.lb > 0) = 0; m2.ub(m2.ub < 0) = 0;
    m2.c(V.atpIdx) = 1;
    s2 = optimizeCbModel(m2, 'max');
    if s2.stat == 1
        V.egcATP = s2.f;
        if s2.f > cfg.tol
            msg = sprintf('ATP produced with all uptakes closed (%.3g): energy-generating cycle', s2.f);
            if cfg.failOnEGC, V.reasons{end+1} = msg; else, V.warnings{end+1} = msg; end
        end
    end
else
    V.warnings{end+1} = 'no ATP demand reaction: EGC test and Pareto skipped';
end

% weights for L1w / L2w
V.w = expressionWeights(model, cfg);
if isempty(V.w)
    V.w = ones(numel(model.rxns), 1);
    V.warnings{end+1} = 'no usable expressionRxns: L1w/L2w with uniform weights';
end

if ~isempty(V.reasons), V.status = 'FAIL'; elseif ~isempty(V.warnings), V.status = 'WARN'; end
end

function w = expressionWeights(model, cfg)
w = [];
if ~isfield(model, 'expressionRxns') || numel(model.expressionRxns) ~= numel(model.rxns), return, end
e = double(model.expressionRxns(:));
switch cfg.weightScheme
    case 'original'                                   % Fluxes_calculation.m
        e(isnan(e)) = 0;
        w = e + abs(min(e)) + 1e-6;
    case 'lowExpressionCost'
        ok = isfinite(e);
        if sum(ok) < 2, return, end
        e = log1p(max(e, 0));
        if max(e(ok)) <= min(e(ok)), return, end
        en = (e - min(e(ok))) ./ (max(e(ok)) - min(e(ok)));
        w = 1 - 0.99 * en;                            % in [0.01, 1]
        w(~ok) = median(w(ok));                       % reactions without GPR: median cost
    otherwise
        error('unknown cfg.weightScheme %s', cfg.weightScheme);
end
end

% =============================================================================================
function sol = computeSolutions(model, V, cfg)
n = numel(model.rxns);
sol = emptySol(cfg, n);
bmin = cfg.fracOpt * V.maxBiomass;
for s = 1:numel(cfg.solNames)
    try
        switch cfg.solNames{s}
            case 'FBA',  v = V.vFBA; st = V.fbaStat;
            case 'pFBA', [v, st] = weightedL1(model, V.bioIdx, bmin, ones(n, 1));
            case 'L1w',  [v, st] = weightedL1(model, V.bioIdx, bmin, V.w);
            case 'L2',   [v, st] = weightedL2(model, V.bioIdx, bmin, ones(n, 1));
            case 'L2w',  [v, st] = weightedL2(model, V.bioIdx, bmin, V.w);
        end
        if st == 1 && numel(v) == n
            sol(s).v = v(:); sol(s).stat = 1;
        else
            sol(s).stat = st; sol(s).msg = sprintf('stat = %d', st);
        end
    catch ME
        sol(s).stat = -99; sol(s).msg = ME.message;
    end
end
end

function [v, stat] = weightedL1(model, bioIdx, bioMin, w)
% min sum w.*|v|  s.t.  S v = b, lb <= v <= ub, v_bio >= bioMin   (LP with t >= |v|)
[m, n] = size(model.S);
lb = model.lb; ub = model.ub; lb(bioIdx) = max(lb(bioIdx), bioMin);
I = speye(n);
LP.A = [model.S, sparse(m, n); I, -I; -I, -I];
LP.b = [model.b; zeros(2 * n, 1)];
LP.csense = [model.csense; repmat('L', 2 * n, 1)];
LP.c = [zeros(n, 1); w(:)];
LP.lb = [lb; zeros(n, 1)];
LP.ub = [ub; max(abs(lb), abs(ub))];
LP.osense = 1;
s = solveCobraLP(LP);
stat = s.stat; v = [];
if stat == 1, v = s.full(1:n); end
end

function [v, stat] = weightedL2(model, bioIdx, bioMin, w)
% min 1/2 sum w.*v.^2  s.t.  S v = b, lb <= v <= ub, v_bio >= bioMin   (QP)
n = numel(model.rxns);
lb = model.lb; ub = model.ub; lb(bioIdx) = max(lb(bioIdx), bioMin);
QP.A = model.S; QP.b = model.b; QP.csense = model.csense;
QP.F = spdiags(w(:), 0, n, n); QP.c = zeros(n, 1);
QP.lb = lb; QP.ub = ub; QP.osense = 1;
s = solveCobraQP(QP);
stat = s.stat; v = [];
if stat == 1, v = s.full(1:n); end
end

% =============================================================================================
function [M, SA] = computeMetrics(model, sol, V, cfg)
% metrics of Fluxes_calculation.m, one struct per solution
L = cfg.lists; tol = cfg.tol;
ssNames = unique(model.subSystems(~cellfun(@isempty, model.subSystems)));
SA = struct('names', {ssNames}, 'values', nan(numel(ssNames), numel(sol)));
atpRow = find(strcmp(model.mets, cfg.atpMet));
tcaBool = strcmp(model.subSystems, 'Citric acid cycle');
bioFBA = V.maxBiomass;
glnIdx = find(strcmp(model.rxns, cfg.glnRxn)); glcIdx = find(strcmp(model.rxns, cfg.glcRxn));
M = cell(1, numel(sol));
for s = 1:numel(sol)
    v = sol(s).v; r = struct();
    if all(isnan(v))
        for k = 1:numel(cfg.metricRoots), r.(cfg.metricRoots{k}) = NaN; end
        r.Oncomet = nan(1, numel(L.oncomet)); M{s} = r; continue
    end
    r.CU = calcCU(v, model, L.cu);
    r.EA = calcEA(v, model, L.ea, cfg.EAcoeffs);
    gly = calcRobustSum(v, model, L.warburgGly); oxp = calcRobustSum(v, model, L.warburgOx);
    r.WarburgIndex = gly.value - oxp.value;
    if ~isempty(atpRow)
        bal = full(model.S(atpRow, :))' .* v;
        r.ATPProduction = sum(bal(bal > 1e-9)); r.ATPConsumption = abs(sum(bal(bal < -1e-9)));
    else
        r.ATPProduction = NaN; r.ATPConsumption = NaN;
    end
    a = calcRobustSum(v, model, L.redoxNAD); b = calcRobustSum(v, model, L.redoxNADH);
    r.RedoxIndex = ternary((a.coverage + b.coverage) / 2 >= 0.5, a.value - b.value, NaN);
    cnt = 0; nF = 0;
    for k = 1:numel(L.mfi)
        j = find(strcmp(model.rxns, L.mfi{k}));
        if ~isempty(j), nF = nF + 1; cnt = cnt + (v(j) < -cfg.mfiThreshold); end
    end
    r.MFI = ternary(nF >= 2, cnt, NaN);
    a = calcRobustSum(v, model, L.anab); b = calcRobustSum(v, model, L.catab);
    r.AnabolismScore = ternary(a.coverage >= 0.5 && b.coverage >= 0.5 && b.value > tol, a.value / b.value, NaN);
    onc = nan(1, numel(L.oncomet));
    for k = 1:numel(L.oncomet)
        j = find(strcmp(model.rxns, L.oncomet{k}));
        if ~isempty(j), onc(k) = max(0, v(j)); end
    end
    r.Oncomet = onc;
    a = calcRobustSum(v, model, L.nadph);
    r.NADPHdemand = ternary(a.coverage >= 0.5, a.value / bioFBA, NaN);
    r.TCA_completeness = ternary(any(tcaBool), sum(abs(v(tcaBool)) > tol) / sum(tcaBool), NaN);
    a = calcRobustSum(v, model, L.lipSat);   r.LipidSat   = ternary(a.coverage >= 0.5, a.value, NaN);
    a = calcRobustSum(v, model, L.lipUnsat); r.LipidUnsat = ternary(a.coverage >= 0.5, a.value, NaN);
    a = calcRobustSum(v, model, L.lipPL);    r.LipidPL    = ternary(a.coverage >= 0.5, a.value, NaN);
    if ~isempty(glnIdx) && ~isempty(glcIdx)
        fq = abs(min(0, v(glnIdx))); fc = abs(min(0, v(glcIdx)));
        r.GlnDependence = ternary(fq + fc > tol, fq / (fq + fc), NaN);
    else
        r.GlnDependence = NaN;
    end
    for k = 1:numel(ssNames)
        sel = strcmp(model.subSystems, ssNames{k});
        SA.values(k, s) = sum(abs(v(sel))) / sum(sel);
    end
    M{s} = r;
end
end

function result = calcRobustSum(v, model, ids)
[tf, loc] = ismember(ids, model.rxns);
result.value = sum(abs(v(loc(tf))));
result.coverage = mean(tf);
end

function CU = calcCU(v, model, ids)
% carbon uptake: for each listed exchange with uptake (v < 0), carbons of its substrates x coef x |v|
CU = 0;
for j = 1:numel(ids)
    r = find(strcmp(model.rxns, ids{j}));
    if isempty(r) || v(r) >= 0, continue, end
    mi = find(model.S(:, r) < 0);
    for q = mi(:)'
        CU = CU + countCarbons(model.metFormulas{q}) * (-full(model.S(q, r))) * (-v(r));
    end
end
end

function EA = calcEA(v, model, ids, coeffs)
% enzyme allocation: sum coeff_k * v_k over amino-acid exchanges (signed flux, as in the original)
EA = 0;
for k = 1:numel(ids)
    r = find(strcmp(model.rxns, ids{k}));
    if ~isempty(r), EA = EA + coeffs(k) * v(r); end
end
end

function n = countCarbons(formula)
% number of C atoms; 'C' followed by a lowercase letter (Cl, Co, Ca, Cu...) is another element
n = 0;
if isempty(formula) || ~ischar(formula), return, end
tok = regexp(formula, 'C(?![a-z])(\d*)', 'tokens');
for k = 1:numel(tok)
    if isempty(tok{k}{1}), n = n + 1; else, n = n + str2double(tok{k}{1}); end
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

% =============================================================================================
function P = paretoSurface(model, V, cfg)
% Biomass x ATP grid (10 x 10): biomass >= fb*max, ATP demand >= fa*max(ATP | biomass), minimise carbon uptake
if isempty(V.atpIdx), P = []; return, end
b = V.bioIdx; a = V.atpIdx;
[m, n] = size(model.S);
cuIdx = find(ismember(model.rxns, cfg.lists.cu));
cC = zeros(n, 1);                                     % carbons per unit uptake
for j = cuIdx(:)'
    mi = find(model.S(:, j) < 0);
    for q = mi(:)', cC(j) = cC(j) + countCarbons(model.metFormulas{q}) * (-full(model.S(q, j))); end
end
fb = linspace(cfg.paretoMinFrac, 1, cfg.paretoN(1));
fa = linspace(cfg.paretoMinFrac, 1, cfg.paretoN(2));
P = nan(numel(fb) * numel(fa), 5); r = 0;
for i = 1:numel(fb)
    m1 = model; m1.lb(b) = max(m1.lb(b), fb(i) * V.maxBiomass); m1.c(a) = 1;
    s = optimizeCbModel(m1, 'max');
    for k = 1:numel(fa)
        r = r + 1; P(r, 1) = r;
        if s.stat ~= 1, continue, end
        lb = m1.lb; ub = m1.ub; lb(a) = max(lb(a), fa(k) * s.f);
        % variables [v; u], u >= -v (uptake part of exchanges written met[e] <=>), minimise cC'*u
        LP.A = [m1.S, sparse(m, n); -speye(n), -speye(n)];
        LP.b = [m1.b; zeros(n, 1)];
        LP.csense = [m1.csense; repmat('L', n, 1)];
        LP.c = [zeros(n, 1); cC];
        LP.lb = [lb; zeros(n, 1)]; LP.ub = [ub; max(abs(lb), abs(ub))];
        LP.osense = 1;
        q = solveCobraLP(LP);
        if q.stat ~= 1, continue, end
        v = q.full(1:n);
        P(r, 2:5) = [v(b), v(a), calcCU(v, model, cfg.lists.cu), calcEA(v, model, cfg.lists.ea, cfg.EAcoeffs)];
    end
end
P = P(~all(isnan(P(:, 2:5)), 2), :);
end

% =============================================================================================
function S = assembleOutputs(jobs, pmDir, cfg, logFile)
nS = numel(cfg.solNames);
Rs = {};
for i = 1:numel(jobs)
    f = fullfile(pmDir, [jobs(i).id '.mat']);
    if isfile(f), t = load(f, 'R'); Rs{end+1} = t.R; end %#ok<AGROW>
end
nM = numel(Rs);
allR = {}; allSS = {};
for i = 1:nM
    if isfield(Rs{i}, 'rxns'), allR = union(allR, Rs{i}.rxns); end
    if isfield(Rs{i}, 'SA') && ~isempty(Rs{i}.SA), allSS = union(allSS, Rs{i}.SA.names); end
end
allR = allR(:); allSS = allSS(:); nR = numel(allR); nSS = numel(allSS);

% ---- columns (same order as Fluxes_calculation.m)
H = {}; D = cell(0, 4);
for s = 1:nS
    for k = 1:nR
        h = sprintf('Flux_%s_%s', sanitize(allR{k}), cfg.solNames{s});
        H{end+1} = h; D(end+1, :) = {h, 'flux', cfg.solNames{s}, allR{k}}; %#ok<AGROW>
    end
end
scalA = {'CU', 'EA', 'WarburgIndex', 'ATPConsumption', 'ATPProduction'};
scalB = {'RedoxIndex', 'MFI', 'AnabolismScore'};
scalC = {'NADPHdemand', 'TCA_completeness', 'LipidSat', 'LipidUnsat', 'LipidPL', 'GlnDependence'};
[H, D] = addRoots(H, D, scalA, scalA, 'metric', cfg.solNames);
[H, D] = addRoots(H, D, cellfun(@(x) ['SA_' sanitize(x)], allSS, 'UniformOutput', false), allSS, 'subsystem', cfg.solNames);
[H, D] = addRoots(H, D, scalB, scalB, 'metric', cfg.solNames);
[H, D] = addRoots(H, D, cellfun(@(x) ['Oncomet_' x], cfg.oncometNames, 'UniformOutput', false), cfg.lists.oncomet, 'oncometabolite', cfg.solNames);
[H, D] = addRoots(H, D, scalC, scalC, 'metric', cfg.solNames);
extra = {'nRxns', 'nMets', 'nGenes', 'MaxBiomass_FBA', 'QC_pass'};
for k = 1:numel(extra), H{end+1} = extra{k}; D(end+1, :) = {extra{k}, 'qc', '', ''}; end %#ok<AGROW>

% ---- values
X = nan(nM, numel(H)); ids = cell(nM, 1);
qc = cell(nM, 10); par = []; parId = {};
for i = 1:nM
    R = Rs{i}; V = R.V; ids{i} = R.id; c = 0;
    ok = isfield(R, 'sol') && ~strcmp(V.status, 'FAIL');
    for s = 1:nS
        col = nan(nR, 1);
        if ok && ~all(isnan(R.sol(s).v))
            col(:) = 0;                                   % reaction absent from this model: no flux
            [~, loc] = ismember(R.rxns, allR); col(loc) = R.sol(s).v;
        end
        X(i, c + (1:nR)) = col; c = c + nR;
    end
    for k = 1:numel(scalA), for s = 1:nS, c = c + 1; X(i, c) = metric(R, scalA{k}, s); end, end
    for k = 1:nSS, for s = 1:nS
        c = c + 1;
        if ok && ~isempty(R.SA) && ~all(isnan(R.sol(s).v))
            j = find(strcmp(R.SA.names, allSS{k}));
            if isempty(j), X(i, c) = 0; else, X(i, c) = R.SA.values(j, s); end
        end
    end, end
    for k = 1:numel(scalB), for s = 1:nS, c = c + 1; X(i, c) = metric(R, scalB{k}, s); end, end
    for k = 1:numel(cfg.oncometNames), for s = 1:nS
        c = c + 1;
        if ok && ~isempty(R.metrics), X(i, c) = R.metrics{s}.Oncomet(k); end
    end, end
    for k = 1:numel(scalC), for s = 1:nS, c = c + 1; X(i, c) = metric(R, scalC{k}, s); end, end
    sz = [NaN NaN NaN]; if isfield(R, 'size'), sz = R.size; end
    mb = getf(V, 'maxBiomass', NaN);
    X(i, c + (1:5)) = [sz, mb, double(~strcmp(V.status, 'FAIL'))];
    solStat = '';
    if isfield(R, 'sol'), solStat = strjoin(arrayfun(@(q) sprintf('%s=%g', q.name, q.stat), R.sol, 'UniformOutput', false), ';'); end
    qc(i, :) = {R.id, V.status, strjoin(getf(V, 'reasons', {}), ' | '), strjoin(getf(V, 'warnings', {}), ' | '), ...
        getf(V, 'relaxationUsed', NaN), getf(V, 'biomassLB', NaN), mb, getf(V, 'egcATP', NaN), solStat, getf(R, 'seconds', NaN)};
    if isfield(R, 'pareto') && ~isempty(R.pareto)
        par = [par; R.pareto]; parId = [parId; repmat({R.id}, size(R.pareto, 1), 1)]; %#ok<AGROW>
    end
end

% ---- write
writeCSV(fullfile(cfg.outDir, 'FeatureMatrix_TumorPhenotype_All.csv'), [{'Model'}, H], ids, X);
writetable(cell2table(qc, 'VariableNames', {'Model', 'Status', 'FailReasons', 'Warnings', 'RelaxationUsed', ...
    'BiomassLB', 'MaxBiomass_FBA', 'EGC_ATP_noUptake', 'SolutionStatus', 'Seconds'}), fullfile(cfg.outDir, 'ModelQC.csv'));
writetable(cell2table(D, 'VariableNames', {'Column', 'Type', 'Solution', 'Recon3D_ID_or_Name'}), fullfile(cfg.outDir, 'FeatureDictionary.csv'));
writetable(cell2table([ids, num2cell(X(:, end - 1))], 'VariableNames', {'Model', 'MaxBiomass_FBA'}), fullfile(cfg.outDir, 'MaxBiomass_perModel.csv'));
if cfg.doPareto
    writeCSV(fullfile(cfg.outDir, 'ParetoSurface_100sol.csv'), {'Model', 'Solution', 'Biomass', 'ATP', 'CU', 'EA'}, parId, par);
end
S.nModels = nM;
S.nPass = sum(strcmp(qc(:, 2), 'PASS')); S.nWarn = sum(strcmp(qc(:, 2), 'WARN')); S.nFail = sum(strcmp(qc(:, 2), 'FAIL'));
logmsg(logFile, sprintf('FeatureMatrix: %d models x %d columns', nM, numel(H) + 1));
end

function [H, D] = addRoots(H, D, roots, keys, type, solNames)
% one column per (root, solution), root-major: CU_FBA, CU_pFBA, ..., EA_FBA, ...
for k = 1:numel(roots)
    for s = 1:numel(solNames)
        h = sprintf('%s_%s', roots{k}, solNames{s});
        H{end+1} = h; D(end+1, :) = {h, type, solNames{s}, keys{k}}; %#ok<AGROW>
    end
end
end

function x = metric(R, root, s)
x = NaN;
if isfield(R, 'metrics') && ~isempty(R.metrics) && isfield(R.metrics{s}, root), x = R.metrics{s}.(root); end
end

function writeCSV(fname, header, idcol, X)
% plain CSV (writetable would rename columns such as Flux_EX_glc_D_e_FBA)
fid = fopen(fname, 'w'); assert(fid > 0, 'cannot write %s', fname);
fprintf(fid, '%s\n', strjoin(header, ','));
if ~isempty(X)
    fmt = [repmat('%.10g,', 1, size(X, 2) - 1), '%.10g\n'];
    for i = 1:size(X, 1), fprintf(fid, '%s,', idcol{i}); fprintf(fid, fmt, X(i, :)); end
end
fclose(fid);
end

% =============================================================================================
function cfg = defaultConfig(cfg)
d.solNames = {'FBA', 'pFBA', 'L1w', 'L2', 'L2w'};
d.tol = 1e-6;
d.mfiThreshold = 1e-6;
d.paretoN = [10 10];
d.paretoMinFrac = 0.1;
d.biomassRxn = 'biomass_ob';
d.atpDemandRxn = 'ATPM';
d.atpMet = 'atp[c]';
d.glnRxn = 'EX_gln_L[e]';
d.glcRxn = 'EX_glc_D[e]';
aa = {'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val'};
ex = @(a) ternary(strcmp(a, 'gly'), 'EX_gly[e]', sprintf('EX_%s_L[e]', a));
% reaction lists of Fluxes_calculation.m in Recon3D IDs
L.warburgGly = {'EX_glc_D[e]', 'PFK', 'LDH_L', 'EX_lac_L[e]'};
L.warburgOx  = {'CSm', 'SUCD1m', 'ICDHxm', 'ATPS4minew'};
L.cu = [{'EX_glc_D[e]'}, cellfun(ex, {'ile', 'leu', 'lys', 'met', 'phe', 'thr', 'trp', 'val', 'ala', 'asn', ...
    'gln', 'tyr', 'arg', 'gly', 'ser', 'glu'}, 'UniformOutput', false), {'EX_pyr[e]'}];
L.ea = cellfun(ex, aa, 'UniformOutput', false);
L.redoxNAD  = {'LDH_L', 'LDH_Lm', 'MDH', 'MDHm'};
L.redoxNADH = {'PDHm', 'GAPD', 'ICDHxm', 'ICDHyrm'};
L.mfi = {'EX_gln_L[e]', 'EX_lac_L[e]', 'EX_pyr[e]', 'EX_glc_D[e]'};
L.anab  = {'PGI', 'TALA', 'TKT1', 'RPI'};
L.catab = {'PDHm', 'PCm', 'CSm', 'PYK'};
L.oncomet = {'EX_lac_L[e]', 'EX_succ[e]', 'EX_akg[e]'};
L.nadph = {'MTHFD', 'MTHFDm', 'ICDHxm', 'ICDHyrm'};
L.lipSat   = {'FAOXC140', 'r0653'};
L.lipUnsat = {'C161CPT22', 'C30CPT1'};
L.lipPL    = {'PSSA1_hs', 'CEPTC'};
d.lists = L;
d.oncometNames = {'Lactate', 'Succinate', 'AlphaKG'};
kcatKM = [1.4, 140, 3.2, 4.5, 0.7, 1.9, 2.3, 0.35, 1.1, 2.6, 2.6, 0.8, 3.0, 4.2, 1.5, 0.95, 1.2, 0.6, 1.3, 2.6];
d.EAcoeffs = 10 ./ kcatKM;                            % W ./ kcatKM, W = 10 (original order = aa)
d.metricRoots = {'CU', 'EA', 'WarburgIndex', 'ATPConsumption', 'ATPProduction', 'RedoxIndex', 'MFI', ...
    'AnabolismScore', 'NADPHdemand', 'TCA_completeness', 'LipidSat', 'LipidUnsat', 'LipidPL', 'GlnDependence'};
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(cfg, fn{i}), cfg.(fn{i}) = d.(fn{i}); end
end
end

function setSolvers(cfg)
changeCobraSolver(cfg.lpSolver, 'LP', 0);
changeCobraSolver(cfg.qpSolver, 'QP', 0);
end

function s = sanitize(str)
s = regexprep(str, '\[(\w+)\]$', '_$1');
s = regexprep(s, '[^A-Za-z0-9]+', '_');
s = regexprep(s, '^_+|_+$', '');
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function saveResult(f, R) %#ok<INUSD>
save(f, 'R', '-v7');
end

function logmsg(file, msg)
line = sprintf('[%s] %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg);
fprintf('%s\n', line);
fid = fopen(file, 'a'); if fid > 0, fprintf(fid, '%s\n', line); fclose(fid); end
end
