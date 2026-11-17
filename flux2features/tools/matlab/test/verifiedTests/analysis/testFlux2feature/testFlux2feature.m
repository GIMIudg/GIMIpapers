% The COBRAToolbox: testFlux2feature.m
%
% Purpose:
%     - tests flux2feature: every descriptor of a three-reaction model, of a
%       two-compartment network and of a toy exchanging glucose, glutamine,
%       lactate and CO2 (with peroxide made and removed two ways) against
%       values worked by hand; the absence rule (a missing requirement gives
%       NaN with its reason, never 0); matrices of states and the comparison
%       with a reference state; the tiers of the release; and a sampled state
%       of RBC-GEM against the reference values of the Python implementation
%       (all 179 descriptors to 1e-5 relative), which the analysis notebook
%       of the study checks against its own engine. No solver is needed.
%
% Authors:
%     - Germán Preciat, Ramyet Sotelo Rodriguez, 2026
%

global CBTDIR

% save the current path and switch to the test path
currentDir = cd(fileparts(which(mfilename)));

% get the path of the test folder
testPath = pwd;

% set the tolerances: hand-worked values, and agreement with Python
tol = 1e-12;
referenceTol = 1e-5;

fprintf(' -- Running testFlux2feature ... \n');

% a three-reaction model:  EX_a -> a -> b -> EX_b
S = [ 1 -1  0; ...
      0  1 -1];
model.S = sparse(S);
model.rxns = {'EX_a'; 'R1'; 'EX_b'};
model.mets = {'a_c'; 'b_c'};
model.lb = [-10; 0; 0];
model.ub = [0; 10; 10];
model.subSystems = {''; 'Glycolysis / Gluconeogenesis'; ''};

v = [-4; 4; 4];
f = flux2feature(model, v, 'Which', 'computable');

assert(abs(f.int_L1 - 4) < tol, 'int_L1');
assert(abs(f.ext_L1 - 8) < tol, 'ext_L1');
assert(abs(f.ex_total_uptake - 4) < tol, 'ex_total_uptake');
assert(abs(f.ex_total_secretion - 4) < tol, 'ex_total_secretion');
assert(abs(f.cpl_boundary_share - 8/12) < tol, 'cpl_boundary_share');
assert(abs(f.path_glycolysis_activity - 4) < tol, 'path_glycolysis_activity');
assert(abs(f.path_glycolysis_share - 1) < tol, 'path_glycolysis_share');
% R1 runs at 4 of its 0..10 span and sits at neither bound
assert(abs(f.cap_int_mean_saturation - 0.4) < tol, 'cap_int_mean_saturation');
assert(abs(f.cap_int_at_bound - 0) < tol, 'cap_int_at_bound');
fprintf('  toy model: %d descriptors, values as computed by hand\n', ...
        numel(fieldnames(f)));

% the absence rule
info = flux2feature(model);
idx = find(strcmp({info.name}, 'rxn_g6pd'), 1);
assert(~isempty(idx) && ~info(idx).computable, 'rxn_g6pd declared and not computable');
assert(~isempty(strfind(info(idx).reason, 'G6PDH2')), 'reason names the missing reaction');
every = flux2feature(model, v, 'Which', 'all');
assert(isnan(every.rxn_g6pd), 'an absent reaction gives NaN, not 0');
assert(isnan(every.eff_ppp_flux_split), 'a ratio of absent parts gives NaN');
fprintf('  absence rule: missing requirements give NaN and carry a reason\n');

% many states at once, and the reference comparison
V = [v, 2*v];
many = flux2feature(model, V, 'Which', 'computable');
assert(numel(many.int_L1) == 2, 'one value per state');
assert(abs(many.int_L1(1) - f.int_L1) < tol, 'first state matches the single call');
assert(abs(many.int_L1(2) - 8) < tol, 'second state scales');
[~, ~, cmp] = flux2feature(model, V, v, 'Which', 'computable');
k = find(strcmp({cmp.name}, 'int_L1'), 1);
assert(abs(cmp(k).log2ratio(2) - 1) < tol, 'log2 ratio against the reference');
fprintf('  matrices and the reference comparison behave\n');

% the tiers of the release: the default is the core, an alias only on request
if ~all(strcmp({info.tier}, 'unclassified'))
    core = flux2feature(model, v);
    coreNames = {info(strcmp({info.tier}, 'core')).name}';
    assert(isequal(sort(fieldnames(core)), sort(coreNames)), 'the default output is the core');
    assert(~isfield(core, 'int_L1') && isfield(every, 'int_L1'), ...
          'a generic statistic is left out of the core');
    aliases = info(strcmp({info.tier}, 'alias'));
    assert(all(~cellfun(@isempty, {aliases.canonical})), 'every alias names its canonical form');
    both = flux2feature(model, v, 'Which', {'core', 'module'});
    assert(numel(fieldnames(both)) > numel(fieldnames(core)), 'the modules add to the core');
    fprintf('  tiers: %d core, %d module, %d generic, %d alias\n', ...
            nnz(strcmp({info.tier}, 'core')), nnz(strcmp({info.tier}, 'module')), ...
            nnz(strcmp({info.tier}, 'generic')), numel(aliases));
end

% glycolysis in the cytosol; PDH and oxidative phosphorylation in the
% mitochondrion; an ATP/ADP translocase, a pyruvate and a CO2 carrier between
% them (the same network as the Python tests and the notebook)
clear model
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
f = flux2feature(model, v, 'Which', 'all');
expected = { ...
    'cmp_atp_made_in_mitochondrion', 6/8; 'cmp_atp_made_in_cytosol', 2/8; ...
    'cmp_nadh_made_in_mitochondrion', 2/4; 'cmp_nadh_made_in_cytosol', 2/4; ...
    'cmp_nadph_made_in_mitochondrion', 1.5/2; 'cmp_nadph_made_in_cytosol', 0.5/2; ...
    'cmp_co2_made_in_mitochondrion', 3.5/4; 'cmp_co2_made_in_cytosol', 0.5/4; ...
    'cmp_flux_in_mitochondrion', 5.5/30; 'cmp_flux_in_cytosol', 14/30; ...
    'cmp_flux_between_compartments', 10.5/30};
for k = 1:size(expected, 1)
    assert(isfield(f, expected{k, 1}), [expected{k, 1} ' computed']);
    assert(abs(f.(expected{k, 1}) - expected{k, 2}) < tol, expected{k, 1});
end
% a model without organelles declares them and says why it cannot compute them
one.S = sparse([1 -1 0; 0 1 -1]);
one.rxns = {'EX_a'; 'R1'; 'EX_b'};
one.mets = {'a_c'; 'b_c'};
info = flux2feature(one);
idx = find(strcmp({info.name}, 'cmp_flux_between_compartments'), 1);
assert(~isempty(idx) && ~info(idx).computable, 'compartment flux unavailable without organelles');
fprintf('  compartments: %d descriptors, values as computed by hand\n', size(expected, 1));

% glucose and glutamine in, lactate and CO2 out; peroxide made by an oxidase
% (3) and removed by catalase (one turn, 2) and glutathione peroxidase (two
% turns, 2); glutathione reductase uses 2 NADPH, G6PD and GND make 5 and 1,
% an oxidase uses the other 2 (the same toy as the Python tests; the
% stoichiometry is not balanced, the descriptors read only the rows they name)
clear model
mets = {'glc__D_e','gln__L_e','lac__L_e','co2_e','h2o2_c','o2_c','h2o_c', ...
        'gthrd_c','gthox_c','nadph_c','nadp_c'};
formulas = {'C6H12O6','C5H10N2O3','C3H5O3','CO2','H2O2','O2','H2O', ...
            'C10H16N3O6S','C20H30N6O12S2','C21H26N7O17P3','C21H25N7O17P3'};
entries = { ...
    'EX_glc__D_e', {'glc__D_e', -1}, -2; ...
    'EX_gln__L_e', {'gln__L_e', -1}, -1; ...
    'EX_lac__L_e', {'lac__L_e', -1}, 3; ...
    'EX_co2_e',    {'co2_e', -1}, 8; ...
    'OXIDASE',     {'o2_c', -1; 'h2o2_c', 1}, 3; ...
    'CAT',         {'h2o2_c', -2; 'o2_c', 1; 'h2o_c', 2}, 1; ...
    'GTHP',        {'h2o2_c', -1; 'gthrd_c', -2; 'gthox_c', 1; 'h2o_c', 2}, 2; ...
    'GTHOy',       {'gthox_c', -1; 'nadph_c', -1; 'gthrd_c', 2; 'nadp_c', 1}, 2; ...
    'G6PDH2',      {'nadp_c', -1; 'nadph_c', 1}, 5; ...
    'GND',         {'nadp_c', -1; 'nadph_c', 1}, 1; ...
    'NADPHOX',     {'nadph_c', -1; 'nadp_c', 1}, 2};
S = zeros(numel(mets), size(entries, 1));
for j = 1:size(entries, 1)
    pairs = entries{j, 2};
    for k = 1:size(pairs, 1)
        S(strcmp(mets, pairs{k, 1}), j) = pairs{k, 2};
    end
end
model.S = sparse(S);
model.rxns = entries(:, 1);
model.mets = mets(:);
model.metFormulas = formulas(:);
v = cell2mat(entries(:, 3));
f = flux2feature(model, v, 'Which', 'all');
expected = { ...
    'elem_n_uptake', 2; 'elem_n_from_amino_acids', 1; ...
    'elem_c_from_glucose', 12/17; 'elem_c_from_amino_acids', 5/17; ...
    'elem_c_from_lactate_pyruvate', 0; 'elem_c_to_lactate', 9/17; 'elem_c_to_co2', 8/17; ...
    'div_c_effective_sources', 289/169; 'div_c_dominant_source', 12/17; ...
    'aa_glutamine_net_uptake', 1; 'aa_total_uptake', 1; 'aa_n_taken_up', 1; ...
    'ros_h2o2_production', 3; 'ros_h2o2_removal', 4; 'ros_h2o2_removal_per_production', 4/3; ...
    'ros_h2o2_removed_by_catalase', 0.5; 'ros_h2o2_removed_by_glutathione_peroxidase', 0.5; ...
    'ros_glutathione_regeneration_cover', 6/2; 'ros_glutathione_cycle_closure', 1; ...
    'conv_lactate_per_glucose_carbon', 9/12; 'conv_nadph_per_carbon', 6/17; ...
    'bud_nadph_to_glutathione_reduction', 2/4; 'bud_oxppp_cover_of_glutathione_reduction', 6/2};
for k = 1:size(expected, 1)
    assert(isfield(f, expected{k, 1}), [expected{k, 1} ' declared']);
    assert(abs(f.(expected{k, 1}) - expected{k, 2}) <= 1e-9 * max(1, abs(expected{k, 2})), ...
           expected{k, 1});
end
info = flux2feature(model);
idx = find(strcmp({info.name}, 'ros_h2o2_removed_by_peroxiredoxin'), 1);
assert(~isempty(idx) && ~info(idx).computable, 'no peroxiredoxin, so not computable');
assert(isnan(f.ros_h2o2_removed_by_peroxiredoxin), 'absent is NaN, never 0');
fprintf(['  elements, amino acids, oxidative defence, budgets and conversion: ' ...
         '%d descriptors, values as computed by hand\n'], size(expected, 1));

% a sampled state of RBC-GEM against the reference values of the Python
% implementation. RBC-GEM_flux2feature.mat is RBC-GEM (Haiman et al., PLoS
% Comput Biol 2025) reduced to the fields flux2feature reads; the state and
% the values are written by the study's notebook after checking the Python
% tool against its engine.
modelFile = [testPath filesep 'RBC-GEM_flux2feature.mat'];
sampleFile = [testPath filesep 'reference_state.csv'];
if exist(modelFile, 'file') == 2 && exist(sampleFile, 'file') == 2
    referenceFile = [testPath filesep 'flux2feature_reference.csv'];
    if exist(referenceFile, 'file') ~= 2
        error('flux2feature:reference', 'no reference at %s', referenceFile);
    end

    model = loadModel(modelFile);
    v = loadFirstState(sampleFile, model.rxns);
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
    assert(isempty(missing), sprintf('declared in Python, not in MATLAB: %s', ...
        strjoin(missing, ', ')));
    assert(numel(fieldnames(mine)) == numel(names), 'the same descriptors in both languages');
    assert(worst <= referenceTol, sprintf(['every descriptor within 1e-5 of the reference ' ...
        '(worst: %s at %.2e)'], worstName, worst));
    fprintf(['  RBC-GEM: %d of %d descriptors agree with the Python ' ...
             'implementation to 1e-5 (worst %.2e, %s)\n'], compared, numel(names), ...
            worst, worstName);
else
    fprintf('  (RBC-GEM_flux2feature.mat or reference_state.csv not found: reference check skipped)\n');
end

% output a success message
fprintf('Done.\n');

% change the directory
cd(currentDir)


% ==========================================================================
% helpers
% ==========================================================================
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


function v = loadFirstState(sampleFile, rxns)
% reaction,flux: the reference state, in the order of the model's reactions
handle = fopen(sampleFile, 'r');
if handle < 0, error('flux2feature:samples', 'cannot open %s', sampleFile); end
cleanup = onCleanup(@() fclose(handle));
fgetl(handle);                                   % header
names = {}; v = [];
while true
    line = fgetl(handle);
    if ~ischar(line), break, end
    parts = strsplit(strtrim(line), ',');
    if numel(parts) < 2, continue, end
    names{end+1, 1} = parts{1};                                         %#ok<AGROW>
    v(end+1, 1) = str2double(parts{end});                               %#ok<AGROW>
end
assert(isequal(names, rxns(:)), 'the state and the model list the same reactions');
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
