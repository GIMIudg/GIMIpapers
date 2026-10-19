function test_flux2feature(modelFile, sampleFile)
% TEST_FLUX2FEATURE  Check flux2feature, by hand and against the reference.
%
%   TEST_FLUX2FEATURE
%       Runs the hand-checkable part only: a three-reaction model whose every
%       descriptor can be worked out on paper, plus the absence rule.
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
if nargin >= 2
    referenceCheck(modelFile, sampleFile);
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
