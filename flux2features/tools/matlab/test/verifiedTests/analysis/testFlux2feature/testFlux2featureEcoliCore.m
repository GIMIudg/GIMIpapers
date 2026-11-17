% The COBRAToolbox: testFlux2featureEcoliCore.m
%
% Purpose:
%     - tests flux2feature on a reconstruction it was not written for,
%       E. coli core: identifiers mapped through the declarations, descriptor
%       values against the same quantities computed directly from the model
%       and a parsimonious flux balance solution, and the absence rule (no
%       peroxide and no mitochondrion give NaN with a reason, never 0)
%
% Authors:
%     - Germán Preciat, Ramyet Sotelo Rodriguez, 2026
%

global CBTDIR

% save the current path and switch to the test path
currentDir = cd(fileparts(which(mfilename)));

% set the tolerance
tol = 1e-9;

% this test needs an LP solver
solverPkgs = prepareTest('needsLP', true);

fprintf(' -- Running testFlux2featureEcoliCore ... \n');

model = getDistributedModel('ecoli_core_model.mat');
% flux2feature names metabolites as BiGG does (atp_c, glc__D_e)
model.mets = strrep(regexprep(model.mets, '\[(\w+)\]$', '_$1'), '-', '__');

% subsystems as one string per reaction, whichever way COBRA stored them
subs = model.subSystems;
for k = 1:numel(subs)
    if iscell(subs{k})
        if isempty(subs{k}), subs{k} = ''; else, subs{k} = subs{k}{1}; end
    end
end
glycolysis = unique(subs(~cellfun(@isempty, regexpi(subs, 'glycolysis', 'once'))));
pentose = unique(subs(~cellfun(@isempty, regexpi(subs, 'pentose', 'once'))));
assert(~isempty(glycolysis) && ~isempty(pentose), 'E. coli core annotates both pathways');

% the same labels, E. coli core's identifiers
pathways = { ...
    'glycolysis',        glycolysis(:)'; ...
    'pentose_phosphate', pentose(:)'; ...
    'glutathione',       {'Glutathione metabolism'}; ...
    'redox_detox',       {'Reactive species formation and detoxification'}; ...
    'purine',            {'Purine metabolism'}; ...
    'pyrimidine',        {'Pyrimidine metabolism'}; ...
    'transport',         {'Transport, Extracellular', 'Transport, extracellular'}};
named = { ...
    'g6pd','G6PDH2r'; 'pgl','PGL'; 'gnd','GND'; 'catalase','CAT'; 'gpx','GTHP'; ...
    'sod','SPODM'; 'hexokinase','HEX1'; 'pfk','PFK'; 'pyruvate_kinase','PYK'; ...
    'ldh','LDH_D'; 'adenylate_kinase','ADK1'};
declaration = {'Pathways', pathways, 'NamedReactions', named};

% the audit needs no solver
info = flux2feature(model, [], [], declaration{:});
assert(numel(info) == 179, 'every declared descriptor is audited');
reasons = containers.Map({info.name}, {info.reason});
assert(~info(strcmp({info.name}, 'turnover_h2o2_c')).computable, 'no peroxide in E. coli core');
assert(~isempty(reasons('cmp_flux_between_compartments')), 'no organelle, and the reason why');

exchange = strncmp(model.rxns, 'EX_', 3);
internal = ~exchange;
absS = abs(model.S);
inPPP = ismember(subs, pentose) & internal;

for p = 1:numel(solverPkgs.LP)
    solverOK = changeCobraSolver(solverPkgs.LP{p}, 'LP', 0);
    if ~solverOK, continue, end
    fprintf('   Testing flux2feature on E. coli core using %s ... \n', solverPkgs.LP{p});

    % the parsimonious optimum: the smallest total flux at maximal growth
    solution = optimizeCbModel(model, 'max', 'one');
    assert(solution.stat == 1, 'E. coli core grows');
    v = solution.x;

    f = flux2feature(model, v, 'Which', 'all', declaration{:});
    assert(numel(fieldnames(f)) == 179, 'every declared descriptor is returned');

    % values that need nothing named
    intL1 = sum(abs(v(internal)));
    extL1 = sum(abs(v(exchange)));
    assert(abs(f.int_L1 - intL1) <= tol * max(1, intL1), 'int_L1');
    assert(abs(f.ext_L1 - extL1) <= tol * max(1, extL1), 'ext_L1');
    assert(abs(f.cpl_boundary_share - extL1 / (intL1 + extL1)) <= tol, 'cpl_boundary_share');

    % values that need the mapped identifiers
    assert(abs(f.rxn_g6pd - v(strcmp(model.rxns, 'G6PDH2r'))) <= tol, 'rxn_g6pd, mapped');
    assert(abs(f.rxn_ldh - v(strcmp(model.rxns, 'LDH_D'))) <= tol, 'rxn_ldh, mapped');
    ppp = sum(abs(v(inPPP)));
    assert(abs(f.path_pentose_phosphate_activity - ppp) <= tol * max(1, ppp), ...
           'path_pentose_phosphate_activity');
    assert(abs(f.path_pentose_phosphate_share - ppp / intL1) <= tol, ...
           'path_pentose_phosphate_share');
    nadph = abs(v(internal))' * full(absS(strcmp(model.mets, 'nadph_c'), internal))';
    nadp = abs(v(internal))' * full(absS(strcmp(model.mets, 'nadp_c'), internal))';
    assert(abs(f.turnover_nadph_c - nadph) <= tol * max(1, nadph), 'turnover_nadph_c');
    assert(abs(f.cof_nadph_nadp_share - nadph / (nadph + nadp)) <= tol, 'cof_nadph_nadp_share');

    % the absence rule
    assert(isnan(f.rxn_catalase) && isnan(f.turnover_h2o2_c), 'absent is NaN, never 0');
    assert(isnan(f.cmp_atp_made_in_mitochondrion), 'no mitochondrion, so NaN');
    unmapped = flux2feature(model, v, 'Which', {'rxn_g6pd'});
    assert(isnan(unmapped.rxn_g6pd), 'G6PDH2 is not an E. coli core identifier');

    % the default output is the core
    core = flux2feature(model, v, declaration{:});
    assert(numel(fieldnames(core)) == nnz(strcmp({info.tier}, 'core')), 'the core by default');
end

% output a success message
fprintf('Done.\n');

% change the directory
cd(currentDir)
