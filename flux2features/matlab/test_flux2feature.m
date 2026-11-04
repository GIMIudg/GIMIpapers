function test_flux2feature(modelFile, sampleFile)
% TEST_FLUX2FEATURE  Check flux2feature, by hand and against the reference.
%
%   TEST_FLUX2FEATURE
%       Runs the hand-checkable part: a three-reaction model whose every
%       descriptor can be worked out on paper, the absence rule, and a toy
%       network with a mitochondrion for the compartment descriptors. If the
%       repository's data/RBC-GEM.mat and matlab/reference_state.csv are
%       present, it also runs the reference check below on them.
%
%   TEST_FLUX2FEATURE(MODELFILE, SAMPLEFILE)
%       Also compares against flux2feature_reference.csv, which holds the
%       values the validated Python implementation produces for the first
%       retained state of the healthy reference chain. Pass the RBC-GEM model
%       and the matching sample archive, for example
%
%           test_flux2feature('data/input/RBC-GEM.mat', ...
%                             'results/matlab/samples/modelOr_samples1.mat')
%
%       The reference state is column 1 of that archive. The two
%       implementations were written from the same registry but not from the
%       same code, so agreeing descriptor by descriptor on a real state of a
%       3,275-reaction model is the evidence that they are one instrument.
%
%   Prints what passed and raises on the first disagreement.

fprintf('flux2feature test\n');
toyModel();
compartmentToy();
here = fileparts(mfilename('fullpath'));
if nargin >= 2
    referenceCheck(modelFile, sampleFile);
elseif exist(fullfile(here, '..', 'data', 'RBC-GEM.mat'), 'file') == 2 ...
        && exist(fullfile(here, 'reference_state.csv'), 'file') == 2
    % the repository carries the model and the reference state
    referenceCheck(fullfile(here, '..', 'data', 'RBC-GEM.mat'), ...
                   fullfile(here, 'reference_state.csv'));
else
    fprintf('  (pass a model and a sample archive to run the reference check)\n');
end
fprintf('all checks passed\n');
end


% ==========================================================================
function toyModel()
% A three-reaction model:  EX_a -> a -> b -> EX_b
S = [ 1 -1  0; ...
      0  1 -1];
model.S = sparse(S);
model.rxns = {'EX_a'; 'R1'; 'EX_b'};
model.mets = {'a_c'; 'b_c'};
model.lb = [-10; 0; 0];
model.ub = [0; 10; 10];
model.subSystems = {''; 'Glycolysis / Gluconeogenesis'; ''};

v = [-4; 4; 4];
f = flux2feature(model, v);

check(abs(f.int_L1 - 4) < 1e-12, 'int_L1');
check(abs(f.ext_L1 - 8) < 1e-12, 'ext_L1');
check(abs(f.ex_total_uptake - 4) < 1e-12, 'ex_total_uptake');
check(abs(f.ex_total_secretion - 4) < 1e-12, 'ex_total_secretion');
check(abs(f.cpl_boundary_share - 8/12) < 1e-12, 'cpl_boundary_share');
check(abs(f.path_glycolysis_activity - 4) < 1e-12, 'path_glycolysis_activity');
check(abs(f.path_glycolysis_share - 1) < 1e-12, 'path_glycolysis_share');
% R1 runs at 4 of its 0..10 span and sits at neither bound
check(abs(f.cap_int_mean_saturation - 0.4) < 1e-12, 'cap_int_mean_saturation');
check(abs(f.cap_int_at_bound - 0) < 1e-12, 'cap_int_at_bound');
fprintf('  toy model: %d descriptors, values as computed by hand\n', ...
        numel(fieldnames(f)));

% the absence rule
info = flux2feature(model);
idx = find(strcmp({info.name}, 'rxn_g6pd'), 1);
check(~isempty(idx) && ~info(idx).computable, 'rxn_g6pd declared and not computable');
check(~isempty(strfind(info(idx).reason, 'G6PDH2')), 'reason names the missing reaction');
every = flux2feature(model, v, 'Which', 'all');
check(isnan(every.rxn_g6pd), 'an absent reaction gives NaN, not 0');
check(isnan(every.eff_ppp_flux_split), 'a ratio of absent parts gives NaN');
fprintf('  absence rule: missing requirements give NaN and carry a reason\n');

% many states at once, and the reference comparison
V = [v, 2*v];
many = flux2feature(model, V);
check(numel(many.int_L1) == 2, 'one value per state');
check(abs(many.int_L1(1) - f.int_L1) < 1e-12, 'first state matches the single call');
check(abs(many.int_L1(2) - 8) < 1e-12, 'second state scales');
[~, ~, cmp] = flux2feature(model, V, v);
k = find(strcmp({cmp.name}, 'int_L1'), 1);
check(abs(cmp(k).log2ratio(2) - 1) < 1e-12, 'log2 ratio against the reference');
fprintf('  matrices and the reference comparison behave\n');
end


% ==========================================================================
function compartmentToy()
% Glycolysis in the cytosol; PDH and oxidative phosphorylation in the
% mitochondrion; an ATP/ADP translocase, a pyruvate and a CO2 carrier between
% them. The same network as python/flux2feature.py and the notebook's check.
mets = {'glc_e','glc_c','pyr_c','pyr_m','atp_c','adp_c','atp_m','adp_m', ...
        'nad_c','nadh_c','nad_m','nadh_m','nadp_c','nadph_c','nadp_m','nadph_m', ...
        'co2_c','co2_m','co2_e'};
rxns = {'EX_glc_e','GLCt','GLY','PYRt','PDH','OXPHOS','ANT','ATPM','G6PD_c', ...
        'IDH_m','CO2t_m','CO2t_e','EX_co2_e'};
entries = { ...
    'EX_glc_e', {'glc_e', -1}; ...
    'GLCt',     {'glc_e', -1; 'glc_c', 1}; ...
    'GLY',      {'glc_c', -1; 'adp_c', -2; 'nad_c', -2; 'pyr_c', 2; 'atp_c', 2; 'nadh_c', 2}; ...
    'PYRt',     {'pyr_c', -1; 'pyr_m', 1}; ...
    'PDH',      {'pyr_m', -1; 'nad_m', -1; 'nadh_m', 1; 'co2_m', 1}; ...
    'OXPHOS',   {'nadh_m', -1; 'adp_m', -3; 'nad_m', 1; 'atp_m', 3}; ...
    'ANT',      {'atp_m', -1; 'adp_c', -1; 'atp_c', 1; 'adp_m', 1}; ...
    'ATPM',     {'atp_c', -1; 'adp_c', 1}; ...
    'G6PD_c',   {'glc_c', -1; 'nadp_c', -1; 'nadph_c', 1; 'co2_c', 1}; ...
    'IDH_m',    {'nadp_m', -1; 'nadph_m', 1; 'co2_m', 1}; ...
    'CO2t_m',   {'co2_m', -1; 'co2_c', 1}; ...
    'CO2t_e',   {'co2_c', -1; 'co2_e', 1}; ...
    'EX_co2_e', {'co2_e', -1}};
S = zeros(numel(mets), numel(rxns));
for j = 1:size(entries, 1)
    pairs = entries{j, 2};
    for k = 1:size(pairs, 1)
        S(strcmp(mets, pairs{k, 1}), j) = pairs{k, 2};
    end
end
model.S = sparse(S);
model.rxns = rxns(:);
model.mets = mets(:);
v = [-1.5; 1.5; 1.0; 2.0; 2.0; 2.0; 5.0; 7.0; 0.5; 1.5; 3.5; 4.0; 4.0];
f = flux2feature(model, v);
expected = { ...
    'cmp_atp_made_in_mitochondrion', 6/8; 'cmp_atp_made_in_cytosol', 2/8; ...
    'cmp_nadh_made_in_mitochondrion', 2/4; 'cmp_nadh_made_in_cytosol', 2/4; ...
    'cmp_nadph_made_in_mitochondrion', 1.5/2; 'cmp_nadph_made_in_cytosol', 0.5/2; ...
    'cmp_co2_made_in_mitochondrion', 3.5/4; 'cmp_co2_made_in_cytosol', 0.5/4; ...
    'cmp_flux_in_mitochondrion', 5.5/30; 'cmp_flux_in_cytosol', 14/30; ...
    'cmp_flux_between_compartments', 10.5/30};
for k = 1:size(expected, 1)
    check(isfield(f, expected{k, 1}), [expected{k, 1} ' computed']);
    check(abs(f.(expected{k, 1}) - expected{k, 2}) < 1e-12, expected{k, 1});
end
% a model without organelles declares them and says why it cannot compute them
one.S = sparse([1 -1 0; 0 1 -1]);
one.rxns = {'EX_a'; 'R1'; 'EX_b'};
one.mets = {'a_c'; 'b_c'};
info = flux2feature(one);
idx = find(strcmp({info.name}, 'cmp_flux_between_compartments'), 1);
check(~isempty(idx) && ~info(idx).computable, 'compartment flux unavailable without organelles');
fprintf('  compartments: %d descriptors, values as computed by hand\n', size(expected, 1));
end


% ==========================================================================
function referenceCheck(modelFile, sampleFile)
here = fileparts(mfilename('fullpath'));
referenceFile = fullfile(here, 'flux2feature_reference.csv');
if exist(referenceFile, 'file') ~= 2
    error('flux2feature:reference', 'no reference at %s', referenceFile);
end

model = loadModel(modelFile);
v = loadFirstState(sampleFile, numel(model.rxns));
mine = flux2feature(model, v, 'Which', 'all');
[names, expected] = readReference(referenceFile);

worst = 0; worstName = ''; compared = 0; missing = {};
for k = 1:numel(names)
    if ~isfield(mine, names{k})
        missing{end+1} = names{k};                                      %#ok<AGROW>
        continue
    end
    a = mine.(names{k});
    b = expected(k);
    if isnan(a) && isnan(b), compared = compared + 1; continue, end
    relative = abs(a - b) / max(1, abs(b));
    compared = compared + 1;
    if relative > worst, worst = relative; worstName = names{k}; end
end
if ~isempty(missing)
    error('flux2feature:reference', 'not produced by the MATLAB version: %s', ...
        strjoin(missing, ', '));
end
check(worst <= 1e-5, sprintf(['every descriptor within 1e-5 of the reference ' ...
    '(worst: %s at %.2e)'], worstName, worst));
fprintf(['  reference: %d of %d descriptors agree with the Python ' ...
         'implementation to 1e-5 (worst %.2e, %s)\n'], compared, numel(names), ...
        worst, worstName);
end


function model = loadModel(modelFile)
raw = load(modelFile);
names = fieldnames(raw);
model = [];
for k = 1:numel(names)
    candidate = raw.(names{k});
    if isstruct(candidate) && isfield(candidate, 'S') && isfield(candidate, 'rxns')
        model = candidate;
        break
    end
end
if isempty(model)
    error('flux2feature:model', 'no COBRA model struct in %s', modelFile);
end
end


function v = loadFirstState(sampleFile, nRxns)
if endsWith(sampleFile, '.csv')
    % reaction,flux: the reference state the repository carries
    handle = fopen(sampleFile, 'r');
    cleanup = onCleanup(@() fclose(handle));
    fgetl(handle);
    v = [];
    while true
        line = fgetl(handle);
        if ~ischar(line), break, end
        parts = strsplit(strtrim(line), ',');
        v(end+1, 1) = str2double(parts{end});                           %#ok<AGROW>
    end
    if numel(v) ~= nRxns
        error('flux2feature:samples', '%d fluxes in %s for %d reactions', ...
            numel(v), sampleFile, nRxns);
    end
    fprintf('  reference state: %s\n', sampleFile);
    return
end
raw = load(sampleFile);
names = fieldnames(raw);
v = [];
for k = 1:numel(names)
    candidate = raw.(names{k});
    if isnumeric(candidate) && size(candidate, 1) == nRxns && size(candidate, 2) >= 1
        v = full(double(candidate(:, 1)));
        fprintf('  reference state: %s(:,1) from %s\n', names{k}, sampleFile);
        break
    end
    if isstruct(candidate) && isfield(candidate, 'points') ...
            && size(candidate.points, 1) == nRxns
        v = full(double(candidate.points(:, 1)));
        fprintf('  reference state: points(:,1) from %s\n', sampleFile);
        break
    end
end
if isempty(v)
    error('flux2feature:samples', ...
        'no %d-row flux matrix in %s', nRxns, sampleFile);
end
end


function [names, values] = readReference(file)
handle = fopen(file, 'r');
if handle < 0, error('flux2feature:reference', 'cannot open %s', file); end
cleanup = onCleanup(@() fclose(handle));
fgetl(handle);                                   % header
names = {}; values = [];
while true
    line = fgetl(handle);
    if ~ischar(line), break, end
    parts = strsplit(strtrim(line), ',');
    if numel(parts) < 2, continue, end
    names{end+1} = parts{1};                                            %#ok<AGROW>
    values(end+1) = str2double(parts{2});                               %#ok<AGROW>
end
end


function check(condition, what)
if ~condition
    error('flux2feature:test', 'FAILED: %s', what);
end
end
