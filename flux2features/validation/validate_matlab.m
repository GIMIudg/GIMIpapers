function validate_matlab()
% VALIDATE_MATLAB  Independent validation of the MATLAB flux2feature.
%
%   Run from MATLAB, in any folder:
%       run('<path to flux2features>/validation/validate_matlab.m')
%   or, with validation/ as the current folder:
%       validate_matlab
%
%   The MATLAB twin of validate_python.py. It reports PASS or FAIL for:
%     1. matlab/test_flux2feature.m (values worked by hand, the absence rule,
%        the compartment network and the reference state of RBC-GEM);
%     2. every case in validation/cases/ and two RBC-GEM states (the reference
%        state and the zero vector), written to validation/results/matlab/ for
%        the cross-check with Python (validation/cross_check.py);
%     3. properties on the RBC-GEM state: a vector, a row and a column give the
%        same values; a matrix gives, state by state, what each state gives
%        alone; every declared descriptor is returned with 'Which','all', the
%        unavailable ones as NaN with a reason; duplicate identifiers, non-finite
%        fluxes and a vector of the wrong length are refused;
%     4. the compartment descriptors are the same with BiGG (atp_m) and VMH
%        (atp[m]) identifiers.
%   The report, with the MATLAB version, goes to validation/results/matlab_report.md.
%   It changes nothing in the tool: a FAIL is a finding to document.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'matlab'));
out = fullfile(here, 'results', 'matlab');
if ~exist(out, 'dir'), mkdir(out); end
results = cell(0, 3);
environment = sprintf('MATLAB %s, %s', version, computer);
fprintf('%s\n', environment);

% -- 1. the unit tests ------------------------------------------------------
fprintf('\n1. unit tests\n');
try
    test_flux2feature();
    results = record(results, 'test_flux2feature', true, '');
catch err
    results = record(results, 'test_flux2feature', false, err.message);
end

% -- 2. cases for the cross-check -------------------------------------------
fprintf('\n2. cases written for the cross-check\n');
listing = dir(fullfile(here, 'cases'));
for k = 1:numel(listing)
    if ~listing(k).isdir || listing(k).name(1) == '.', continue, end
    [model, v] = loadCase(fullfile(here, 'cases', listing(k).name));
    writeValues(out, listing(k).name, model, v);
    results = record(results, ['case ' listing(k).name ' mapped'], true, '');
end
model = loadModel(fullfile(root, 'data', 'RBC-GEM.mat'));
state = readState(fullfile(root, 'matlab', 'reference_state.csv'));
writeValues(out, 'rbc_reference_state', model, state);
writeValues(out, 'rbc_zero_vector', model, zeros(size(state)));
results = record(results, 'RBC-GEM reference and zero vector mapped', true, '');

% -- 3. properties ----------------------------------------------------------
fprintf('\n3. properties, on the RBC-GEM state\n');
alone = flux2feature(model, state, 'Which', 'all');
asRow = flux2feature(model, state(:)', 'Which', 'all');
asColumn = flux2feature(model, state(:), 'Which', 'all');
results = record(results, 'vector, row and column give the same values', ...
    sameValues(alone, asRow, 1) && sameValues(alone, asColumn, 1), '');
batch = flux2feature(model, [state(:), 0.5 * state(:), state(:)], 'Which', 'all');
half = flux2feature(model, 0.5 * state, 'Which', 'all');
results = record(results, 'a matrix gives what each state gives alone', ...
    sameValues(alone, batch, 1) && sameValues(half, batch, 2), '');

info = flux2feature(model);
names = {info.name};
results = record(results, sprintf('all %d declared descriptors returned with Which all', ...
    numel(names)), isequal(sort(fieldnames(alone))', sort(names)), '');
unavailable = names(~[info.computable]);
honest = true;
for k = 1:numel(unavailable)
    entry = info(strcmp(names, unavailable{k}));
    honest = honest && isnan(alone.(unavailable{k})) && ~isempty(entry.reason);
end
results = record(results, sprintf('%d unavailable descriptors are NaN with a reason, never 0', ...
    numel(unavailable)), honest, '');

duplicated = model;
duplicated.rxns{end} = duplicated.rxns{1};
results = record(results, 'duplicate reaction identifiers are refused', ...
    refuses(@() flux2feature(duplicated, state)), '');
broken = state;
broken(1) = NaN;
results = record(results, 'non-finite fluxes are refused', ...
    refuses(@() flux2feature(model, broken)), '');
results = record(results, 'a vector of the wrong length is refused', ...
    refuses(@() flux2feature(model, state(1:end-1))), '');

% -- 4. identifier conventions ----------------------------------------------
fprintf('\n4. identifier conventions\n');
[bigg, vb] = loadCase(fullfile(here, 'cases', 'compartment_bigg'));
[vmh, vv] = loadCase(fullfile(here, 'cases', 'compartment_vmh'));
a = flux2feature(bigg, vb, 'Which', 'all');
b = flux2feature(vmh, vv, 'Which', 'all');
fields = fieldnames(a);
compartment = fields(strncmp(fields, 'cmp_', 4));
agree = true;
for k = 1:numel(compartment)
    agree = agree && abs(a.(compartment{k}) - b.(compartment{k})) ...
        <= 1e-12 * max(1, abs(a.(compartment{k})));
end
results = record(results, sprintf(['compartment descriptors (%d): atp_m (BiGG) and ' ...
    'atp[m] (VMH) give the same values'], numel(compartment)), agree, '');

% -- report -----------------------------------------------------------------
passed = sum([results{:, 2}]);
handle = fopen(fullfile(here, 'results', 'matlab_report.md'), 'w');
fprintf(handle, '# MATLAB validation of flux2feature\n\n%s\n\n', environment);
fprintf(handle, '%d of %d checks passed.\n\n| check | result | detail |\n|---|---|---|\n', ...
    passed, size(results, 1));
for k = 1:size(results, 1)
    verdict = 'FAIL';
    if results{k, 2}, verdict = 'PASS'; end
    fprintf(handle, '| %s | %s | %s |\n', results{k, 1}, verdict, results{k, 3});
end
fclose(handle);
fprintf('\n%d of %d checks passed; report in validation/results/matlab_report.md\n', ...
    passed, size(results, 1));
end


% ==========================================================================
function results = record(results, name, ok, detail)
verdict = 'FAIL';
if ok, verdict = 'PASS'; end
if isempty(detail)
    fprintf('  %s  %s\n', verdict, name);
else
    fprintf('  %s  %s  (%s)\n', verdict, name, detail);
end
results(end+1, :) = {name, logical(ok), detail};
end


function ok = refuses(call)
ok = false;
try
    call();
catch
    ok = true;
end
end


function ok = sameValues(reference, values, column)
% Every field of REFERENCE equals column COLUMN of the same field of VALUES.
ok = true;
names = fieldnames(reference);
for k = 1:numel(names)
    x = reference.(names{k});
    y = values.(names{k})(column);
    if isnan(x) && isnan(y), continue, end
    if isnan(x) ~= isnan(y) || abs(x - y) > 1e-9 * max(1, abs(x))
        ok = false;
        fprintf('    %s: %.17g vs %.17g\n', names{k}, x, y);
        return
    end
end
end


function writeValues(out, name, model, v)
values = flux2feature(model, v, 'Which', 'all');
info = flux2feature(model);
handle = fopen(fullfile(out, [name '.csv']), 'w');
fprintf(handle, 'descriptor,value,computable\n');
for k = 1:numel(info)
    x = values.(info(k).name);
    if isnan(x)
        fprintf(handle, '%s,NaN,%d\n', info(k).name, info(k).computable);
    else
        fprintf(handle, '%s,%.17g,%d\n', info(k).name, x, info(k).computable);
    end
end
fclose(handle);
end


function [model, v] = loadCase(folder)
% reactions.csv (reaction,lb,ub,subsystem), stoichiometry.csv
% (metabolite,reaction,coefficient) and flux.csv (reaction,flux).
rows = readRows(fullfile(folder, 'reactions.csv'));
model.rxns = rows(:, 1);
model.lb = str2double(rows(:, 2));
model.ub = str2double(rows(:, 3));
model.subSystems = rows(:, 4);
triplets = readRows(fullfile(folder, 'stoichiometry.csv'));
[~, first] = unique(triplets(:, 1), 'stable');
model.mets = triplets(first, 1);
S = zeros(numel(model.mets), numel(model.rxns));
for k = 1:size(triplets, 1)
    S(strcmp(model.mets, triplets{k, 1}), strcmp(model.rxns, triplets{k, 2})) = ...
        str2double(triplets{k, 3});
end
model.S = sparse(S);
fluxes = readRows(fullfile(folder, 'flux.csv'));
v = zeros(numel(model.rxns), 1);
for k = 1:size(fluxes, 1)
    v(strcmp(model.rxns, fluxes{k, 1})) = str2double(fluxes{k, 2});
end
end


function rows = readRows(file)
% A small CSV reader: a header, then comma-separated fields without quotes.
handle = fopen(file, 'r');
if handle < 0, error('validate:read', 'cannot open %s', file); end
cleanup = onCleanup(@() fclose(handle));
header = strsplit(fgetl(handle), ',');
rows = cell(0, numel(header));
while true
    line = fgetl(handle);
    if ~ischar(line), break, end
    if isempty(strtrim(line)), continue, end
    parts = strsplit(line, ',', 'CollapseDelimiters', false);
    parts(end+1:numel(header)) = {''};
    rows(end+1, :) = parts(1:numel(header));                            %#ok<AGROW>
end
end


function v = readState(file)
rows = readRows(file);
v = str2double(rows(:, 2));
end


function model = loadModel(file)
raw = load(file);
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
    error('validate:model', 'no COBRA model struct in %s', file);
end
model.rxns = cellstr(model.rxns(:));
end
