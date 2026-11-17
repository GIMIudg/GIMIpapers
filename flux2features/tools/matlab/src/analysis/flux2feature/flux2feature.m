function [features, info, comparison] = flux2feature(model, fluxData, flux2Compare, varargin)
% Maps flux vectors of a constraint-based model to declared functional
% descriptors, the representation Phi_M(v), and audits which descriptors the
% model can support. One file, base MATLAB, no toolbox.
%
% USAGE:
%
%    info = flux2feature(model)
%    features = flux2feature(model, fluxData)
%    [features, info, comparison] = flux2feature(model, fluxData, flux2Compare, varargin)
%
% INPUTS:
%    model:          COBRA model structure with at least the fields `S`,
%                    `rxns` and `mets`. The fields `lb`, `ub`, `subSystems`,
%                    `metFormulas` and `metNames` are used when present.
%
% OPTIONAL INPUTS:
%    fluxData:       Flux states, nRxns x nStates (one state per column, as
%                    `sampleCbModel` returns them) or one flux vector. If
%                    empty or omitted, only the audit is returned.
%    flux2Compare:   Reference state (nRxns x 1), for example a healthy
%                    reference or a wild-type solution. A matrix is reduced
%                    to its median state.
%    varargin:       Name-value pairs:
%
%                      * Which - which descriptors to return: a tier of the
%                        release, 'core' (default), 'module', 'generic' or
%                        'alias', or a cell array of tiers; 'all' for every
%                        declared descriptor; 'computable' for those this
%                        model supports; or a cell array of descriptor names.
%                        A requested descriptor the model cannot support is
%                        NaN, never 0.
%                      * Tolerance - below this a flux counts as zero
%                        (default 1e-6).
%                      * Exchange - logical vector marking the boundary
%                        reactions. By default a reaction is a boundary
%                        reaction if its identifier starts with EX_, and
%                        failing that if it touches exactly one metabolite.
%                      * Pathways, PathwayReactions, NamedReactions,
%                        TrackedExchanges, CofactorPairs - labelled
%                        declarations, as cell arrays with one row per
%                        label ({'g6pd', 'G6PDH2r'; ...}). They are merged
%                        with the defaults: a label given replaces its
%                        default identifiers, a new label adds a descriptor,
%                        so mapping identifiers never removes a coordinate.
%                      * TurnoverMetabolites, EfficiencyRatios, Currency,
%                        CompartmentNames, CompartmentMetabolites - replace
%                        the defaults. The defaults describe human
%                        erythrocyte metabolism (RBC-GEM identifiers).
%
% OUTPUTS:
%    features:       Structure whose fields are descriptor names, each a
%                    1 x nStates row of values (the audit if `fluxData` is
%                    empty).
%    info:           Structure array, one element per declared descriptor:
%                    `name`, `family`, `block`, `units`, `description`,
%                    `requires`, `computable`, `reason` (why not, when it is
%                    not computable), `scope`, `tier` and `canonical` (for an
%                    alias, the descriptor it repeats).
%    comparison:     Structure array, one element per returned descriptor,
%                    with the reference value, the difference and the log2
%                    ratio of every state against `flux2Compare`.
%
% EXAMPLE:
%
%    info = flux2feature(model);                  % what this model supports, and why not
%    nnz([info.computable])
%    f = flux2feature(model, samples);            % the core of Phi, one value per state
%    f.path_pentose_phosphate_activity
%    f = flux2feature(model, samples, 'Which', {'core', 'module'});
%    [f, ~, cmp] = flux2feature(model, samples, vReference, 'Which', 'all');
%    % another reconstruction: map the labels to its identifiers
%    f = flux2feature(model, samples, 'NamedReactions', {'g6pd', 'G6PDH2r'; 'ldh', 'LDH_D'});
%
% NOTE:
%
%    Three rules are enforced. Every descriptor is declared before it is
%    computed (name, units, block, requirements, reading); a descriptor whose
%    requirement the model lacks is NaN, never 0, because zero is a
%    measurement and absence is not; nothing is renamed, reordered or dropped
%    silently. The tiers of the release come from the audit of the study
%    (Supplementary Table S8 of the paper): `core`, functional coordinates
%    with a biochemical reading; `module`, descriptors whose use depends on
%    the model, the medium or the kind of state; `generic`, statistics of the
%    flux vector; `alias`, descriptors that repeat a canonical one. The core
%    is the release's set of interpretable coordinates for the study, not a
%    claim that these are the most important descriptors of any system.
%    Ranges over sampled states are empirical ranges of those states, not
%    confidence intervals or physiological limits.
%
%    Version 1.0.0. The same 179 descriptors, tiers and values as the Python
%    package flux2feature, tested against the same reference values of a
%    sampled state of RBC-GEM, to 1e-5 relative (testFlux2feature.m). Metabolite
%    identifiers are read in BiGG style (atp_c); convert atp[c] first, e.g.
%    model.mets = strrep(regexprep(model.mets, '\[(\w+)\]$', '_$1'), '-', '__');
%
% .. Author: - Germán Preciat, Ramyet Sotelo Rodriguez, Carolina Elizabeth Gómez Márquez, 2026

% ---------------------------------------------------------------- arguments
if nargin < 2, fluxData = []; end
if nargin < 3, flux2Compare = []; end
if nargin >= 3 && ischar(flux2Compare)
    % called as flux2feature(model, V, 'Which', 'all', ...)
    varargin = [{flux2Compare}, varargin];
    flux2Compare = [];
end

opt = parseOptions(varargin);
model = normaliseModel(model, opt);
declared = attachStatus(declare(model, opt));

info = rmfield(declared, 'compute');
if isempty(fluxData)
    features = info;               % flux2feature(model) -> the audit
    comparison = struct([]);
    return
end

V = orientStates(fluxData, numel(model.rxns));          % nStates x nRxns
keep = selectDescriptors(declared, opt.which);
everything = computeAll(declared, V);
features = keepFields(everything, {declared(keep).name});

comparison = struct([]);
if ~isempty(flux2Compare)
    R = orientStates(flux2Compare, numel(model.rxns));
    if size(R, 1) ~= 1
        R = median(R, 1);                               % a reference is one state
    end
    reference = keepFields(computeAll(declared, R), {declared(keep).name});
    names = fieldnames(features);
    comparison = repmat(struct('name', '', 'family', '', 'block', '', ...
        'value', [], 'reference', NaN, 'difference', [], 'log2ratio', []), ...
        numel(names), 1);
    picked = declared(keep);
    for k = 1:numel(names)
        value = features.(names{k});
        base  = reference.(names{k});
        comparison(k).name      = names{k};
        comparison(k).family    = picked(k).family;
        comparison(k).block     = picked(k).block;
        comparison(k).value     = value;
        comparison(k).reference = base;
        comparison(k).difference = value - base;
        ratio = safeDivide(value, repmat(base, size(value)));
        ratio(ratio <= 0) = NaN;                        % a sign change is not a ratio
        comparison(k).log2ratio = log2(ratio);
    end
end
end


% ==========================================================================
% options and model
% ==========================================================================
function opt = parseOptions(args)
opt.which      = 'core';
opt.tolerance  = 1e-6;
opt.exchange   = [];
opt.currency   = {'h_c','h_e','h2o_c','h2o_e','pi_c','pi_e','ppi_c', ...
                  'co2_c','co2_e','na1_c','na1_e','k_c','k_e','cl_c','cl_e'};
% pathway groups: {label, {subsystem strings}}
opt.pathways = { ...
    'glycolysis',        {'Glycolysis / Gluconeogenesis'}; ...
    'pentose_phosphate', {'Pentose phosphate pathway'}; ...
    'glutathione',       {'Glutathione metabolism'}; ...
    'redox_detox',       {'Reactive species formation and detoxification'}; ...
    'purine',            {'Purine metabolism'}; ...
    'pyrimidine',        {'Pyrimidine metabolism'}; ...
    'transport',         {'Transport, extracellular'}};
% groups defined by reaction identity. A reaction claimed here is removed
% from any subsystem-defined group: the Rapoport-Luebering shunt is annotated
% as glycolysis in RBC-GEM and is not glycolysis.
opt.pathwayReactions = {'rapoport_luebering', {'DPGM','DPGase'}};
opt.turnoverMetabolites = {'atp_c','adp_c','amp_c','nadh_c','nad_c','nadph_c', ...
    'nadp_c','gthrd_c','gthox_c','glc__D_c','g6p_c','f6p_c','pyr_c','lac__L_c', ...
    'h2o2_c','o2_c','13dpg_c','r5p_c','ru5p__D_c','6pgc_c'};
opt.namedReactions = { ...
    'g6pd','G6PDH2'; 'pgl','PGL'; 'gnd','GND'; 'catalase','CAT'; 'gpx','GTHP'; ...
    'sod','SPODM'; 'hexokinase','HEX1'; 'pfk','PFK'; 'pyruvate_kinase','PYK'; ...
    'ldh','LDH_L'; 'adenylate_kinase','ADK1'};
opt.trackedExchanges = { ...
    'glc__D_e','glucose'; 'lac__L_e','lactate'; 'o2_e','oxygen'; ...
    'co2_e','carbon_dioxide'; 'h2o2_e','hydrogen_peroxide'; ...
    'gthrd_e','reduced_glutathione'; 'pyr_e','pyruvate'; 'cys__L_e','cysteine'; ...
    'pi_e','phosphate'};
% compartments the compartment shares name, and the metabolites they divide
opt.compartmentNames = {'m', 'mitochondrion'; 'c', 'cytosol'};
opt.compartmentMetabolites = {'atp', 'nadh', 'nadph', 'co2'};
opt.extracellular = {'e', 'b'};
opt.cofactorPairs = { ...
    'nadph_nadp', {'nadph_c','nadp_c'}; ...
    'nadh_nad',   {'nadh_c','nad_c'}; ...
    'glutathione',{'gthrd_c','gthox_c'}};
% ratios, as {label, numerator descriptor, denominator descriptor}
opt.efficiencyRatios = { ...
    'nadph_per_glucose', 'turnover_nadph_c', 'turnover_glc__D_c'; ...
    'nadph_per_atp',     'turnover_nadph_c', 'turnover_atp_c'; ...
    'atp_per_glucose',   'turnover_atp_c',   'turnover_glc__D_c'; ...
    'nadph_per_h2o2',    'turnover_nadph_c', 'turnover_h2o2_c'; ...
    'redox_buffer_ratio','turnover_gthrd_c', 'turnover_gthox_c'; ...
    'ppp_per_glycolysis','path_pentose_phosphate_activity', 'path_glycolysis_activity'; ...
    'ppp_flux_split',    'rxn_g6pd',         'rxn_hexokinase'};

for k = 1:2:numel(args)
    name = lower(args{k});
    switch name
        case 'which',               opt.which = args{k+1};
        case 'tolerance',           opt.tolerance = args{k+1};
        case 'exchange',            opt.exchange = logical(args{k+1}(:))';
        case 'currency',            opt.currency = args{k+1};
        case 'pathways',            opt.pathways = mergeByLabel(opt.pathways, args{k+1}, 1);
        case 'pathwayreactions',    opt.pathwayReactions = mergeByLabel(opt.pathwayReactions, args{k+1}, 1);
        case 'turnovermetabolites', opt.turnoverMetabolites = args{k+1};
        case 'namedreactions',      opt.namedReactions = mergeByLabel(opt.namedReactions, args{k+1}, 1);
        case 'trackedexchanges',    opt.trackedExchanges = mergeByLabel(opt.trackedExchanges, args{k+1}, 2);
        case 'cofactorpairs',       opt.cofactorPairs = mergeByLabel(opt.cofactorPairs, args{k+1}, 1);
        case 'efficiencyratios',    opt.efficiencyRatios = args{k+1};
        case 'compartmentnames',    opt.compartmentNames = args{k+1};
        case 'compartmentmetabolites', opt.compartmentMetabolites = args{k+1};
        otherwise, error('flux2feature:option', 'unknown option "%s"', args{k});
    end
end
end


function table = mergeByLabel(defaults, given, labelColumn)
% A labelled declaration merged with its defaults: a label given replaces the
% default row of that label, a new label is added, the others are kept, so
% mapping a model's identifiers never removes a coordinate.
table = defaults;
if isempty(given), return, end
if size(given, 2) ~= size(defaults, 2)
    error('flux2feature:option', 'expected a cell array with %d columns', size(defaults, 2));
end
for k = 1:size(given, 1)
    hit = find(strcmp(table(:, labelColumn), given{k, labelColumn}), 1);
    if isempty(hit)
        table(end+1, :) = given(k, :);                                  %#ok<AGROW>
    else
        table(hit, :) = given(k, :);
    end
end
end


function model = normaliseModel(model, opt)
% The parts of a reconstruction these descriptors use, in one shape.
required = {'S','rxns','mets'};
for k = 1:numel(required)
    if ~isfield(model, required{k})
        error('flux2feature:model', 'model needs a .%s field', required{k});
    end
end
model.rxns = cellstr(model.rxns(:));
model.mets = cellstr(model.mets(:));
if numel(unique(model.rxns)) ~= numel(model.rxns)
    error('flux2feature:model', 'duplicate reaction identifiers');
end
if numel(unique(model.mets)) ~= numel(model.mets)
    error('flux2feature:model', 'duplicate metabolite identifiers');
end
nRxns = numel(model.rxns);
nMets = numel(model.mets);
if size(model.S, 1) ~= nMets || size(model.S, 2) ~= nRxns
    error('flux2feature:model', 'S is %dx%d, expected %dx%d', ...
        size(model.S, 1), size(model.S, 2), nMets, nRxns);
end
if ~isfield(model, 'lb') || isempty(model.lb), model.lb = -inf(nRxns, 1); end
if ~isfield(model, 'ub') || isempty(model.ub), model.ub =  inf(nRxns, 1); end
model.lb = double(model.lb(:));
model.ub = double(model.ub(:));
model.subSystems = flattenSubSystems(model, nRxns);
model.metFormulas = textField(model, 'metFormulas', nMets);
model.metNames = textField(model, 'metNames', nMets);

if ~isempty(opt.exchange)
    model.exchangeMask = opt.exchange(:)';
else
    prefixed = false(1, nRxns);
    for k = 1:nRxns
        prefixed(k) = strncmp(model.rxns{k}, 'EX_', 3);
    end
    if any(prefixed)
        model.exchangeMask = prefixed;
    else
        model.exchangeMask = (sum(abs(model.S) > 0, 1) == 1);
    end
end
model.exchangeMask = logical(model.exchangeMask(:))';
end


function values = textField(model, field, n)
% A cell array of n character vectors, '' where the model has nothing.
values = repmat({''}, n, 1);
if ~isfield(model, field) || isempty(model.(field)), return, end
raw = model.(field);
if ischar(raw), raw = cellstr(raw); end
if isstring(raw), raw = cellstr(raw); end
for k = 1:min(n, numel(raw))
    if ischar(raw{k}), values{k} = raw{k}; end
end
end


function subs = flattenSubSystems(model, nRxns)
% COBRA has stored subSystems as a cellstr and as a cell of cellstr.
subs = repmat({''}, nRxns, 1);
if ~isfield(model, 'subSystems') || isempty(model.subSystems), return, end
raw = model.subSystems;
if ischar(raw), raw = cellstr(raw); end
for k = 1:min(nRxns, numel(raw))
    entry = raw{k};
    if iscell(entry)
        if isempty(entry), entry = ''; else, entry = entry{1}; end
    end
    if ~ischar(entry), entry = ''; end
    subs{k} = entry;
end
end


function V = orientStates(fluxData, nRxns)
% Return nStates x nRxns, whichever way round the caller had it.
fluxData = double(full(fluxData));
if any(~isfinite(fluxData(:)))
    error('flux2feature:values', 'fluxes contain NaN or infinite values');
end
if isvector(fluxData)
    if numel(fluxData) ~= nRxns
        error('flux2feature:shape', '%d fluxes for %d reactions', ...
            numel(fluxData), nRxns);
    end
    V = fluxData(:)';
    return
end
if size(fluxData, 1) == nRxns
    V = fluxData.';                      % reactions down the rows, as sampleCbModel
elseif size(fluxData, 2) == nRxns
    V = fluxData;                        % one state per row
else
    error('flux2feature:shape', 'fluxData is %dx%d; neither side is %d reactions', ...
        size(fluxData, 1), size(fluxData, 2), nRxns);
end
end


% ==========================================================================
% the registry: declared first, computed later
% ==========================================================================
function d = declare(model, opt)
nRxns = numel(model.rxns);
external = model.exchangeMask;
internal = ~external;
absS = abs(model.S);
keepMet = ~ismember(model.mets, opt.currency);
tol = opt.tolerance;

d = emptyDescriptor();

% -- shape of the state, needing nothing named -------------------------
blocks = {'internal', internal, 'int'; 'external', external, 'ext'};
for b = 1:size(blocks, 1)
    block = blocks{b, 1}; mask = blocks{b, 2}; tag = blocks{b, 3};
    if ~any(mask), continue, end
    d(end+1) = mk([tag '_L1'], 'magnitude', block, 'flux', ...
        sprintf('total %s flux magnitude', block), ...
        @(V, vals) sum(abs(V(:, mask)), 2));                                %#ok<*AGROW>
    d(end+1) = mk([tag '_L2'], 'magnitude', block, 'flux', ...
        sprintf('Euclidean norm of the %s fluxes', block), ...
        @(V, vals) sqrt(sum(V(:, mask).^2, 2)));
    d(end+1) = mk([tag '_Linf'], 'magnitude', block, 'flux', ...
        sprintf('largest single %s flux', block), ...
        @(V, vals) max(abs(V(:, mask)), [], 2));
    d(end+1) = mk([tag '_active_fraction'], 'sparsity', block, 'ratio', ...
        sprintf('fraction of %s reactions carrying flux', block), ...
        @(V, vals) mean(abs(V(:, mask)) > tol, 2));
    d(end+1) = mk([tag '_gini'], 'sparsity', block, 'dimensionless', ...
        sprintf('Gini concentration of the %s flux magnitudes', block), ...
        @(V, vals) giniOf(abs(V(:, mask))));
    d(end+1) = mk([tag '_shannon_entropy'], 'entropy', block, 'nats', ...
        sprintf('entropy of the %s flux shares', block), ...
        @(V, vals) entropyOf(abs(V(:, mask))));
    d(end+1) = mk([tag '_directionality'], 'architecture', block, 'dimensionless', ...
        sprintf('net directional bias of the %s fluxes', block), ...
        @(V, vals) safeDivide( ...
            sum(max(V(:, mask), 0), 2) - sum(max(-V(:, mask), 0), 2), ...
            sum(abs(V(:, mask)), 2)));
    d(end+1) = mk([tag '_cv_abs'], 'distribution', block, 'dimensionless', ...
        sprintf('coefficient of variation of the %s flux magnitudes', block), ...
        @(V, vals) cvOf(abs(V(:, mask))));
    d(end+1) = mk([tag '_top10pct_share'], 'distribution', block, 'ratio', ...
        sprintf('share of %s flux carried by its busiest tenth', block), ...
        @(V, vals) topShare(abs(V(:, mask)), 0.10));
end

% -- network, needing only the stoichiometry ---------------------------
if any(internal)
    M = absS(keepMet, internal) > 0;
    T = absS(keepMet, internal);
    d(end+1) = mk('net_active_metabolites', 'network', 'internal', 'count', ...
        'metabolites touched by a carrying reaction, currency excluded', ...
        @(V, vals) sum((double(abs(V(:, internal)) > tol) * double(M)') > 0, 2), ...
        'stoichiometric matrix');
    d(end+1) = mk('net_top_decile_turnover_share', 'network', 'internal', 'ratio', ...
        ['share of gross metabolite turnover carried by the busiest tenth of ' ...
         'metabolites, currency excluded'], ...
        @(V, vals) topShare(abs(V(:, internal)) * T', 0.10), ...
        'stoichiometric matrix');
end

% -- pathways ----------------------------------------------------------
claimed = false(1, nRxns);
identityGroups = cell(size(opt.pathwayReactions, 1), 2);
for g = 1:size(opt.pathwayReactions, 1)
    members = opt.pathwayReactions{g, 2};
    columns = [];
    for k = 1:numel(members)
        hit = find(strcmp(model.rxns, members{k}), 1);
        if ~isempty(hit), columns(end+1) = hit; end
    end
    claimed(columns) = true;
    identityGroups{g, 1} = opt.pathwayReactions{g, 1};
    identityGroups{g, 2} = columns;
end

groups = [opt.pathways(:, 1), cell(size(opt.pathways, 1), 1)];
for g = 1:size(opt.pathways, 1)
    groups{g, 2} = find(ismember(model.subSystems, opt.pathways{g, 2})' ...
                        & internal & ~claimed);
end
groups = [groups; identityGroups];

for g = 1:size(groups, 1)
    label = groups{g, 1};
    columns = groups{g, 2};
    present = ~isempty(columns);
    reason = '';
    if ~present, reason = sprintf('no reaction belongs to %s in this model', label); end
    d(end+1) = mk(['path_' label '_activity'], 'pathway', 'internal', 'flux', ...
        sprintf('total flux magnitude carried by %s (%d reactions)', label, numel(columns)), ...
        @(V, vals) sum(abs(V(:, columns)), 2), ...
        sprintf('a reaction in %s', label), present, reason, 'conditional');
    d(end+1) = mk(['path_' label '_share'], 'pathway', 'internal', 'ratio', ...
        sprintf('share of internal flux carried by %s', label), ...
        @(V, vals) safeDivide(sum(abs(V(:, columns)), 2), sum(abs(V(:, internal)), 2)), ...
        sprintf('a reaction in %s', label), present, reason, 'conditional');
    d(end+1) = mk(['path_' label '_engagement'], 'pathway', 'internal', 'ratio', ...
        sprintf('fraction of %s reactions carrying flux', label), ...
        @(V, vals) mean(abs(V(:, columns)) > tol, 2), ...
        sprintf('a reaction in %s', label), present, reason, 'conditional');
    two = present && numel(columns) >= 2;
    reasonTwo = reason;
    if present && ~two, reasonTwo = 'only one reaction in the group'; end
    d(end+1) = mk(['path_' label '_directionality'], 'pathway', 'internal', 'dimensionless', ...
        sprintf('net directional bias of %s', label), ...
        @(V, vals) safeDivide( ...
            sum(max(V(:, columns), 0), 2) - sum(max(-V(:, columns), 0), 2), ...
            sum(abs(V(:, columns)), 2)), ...
        sprintf('at least two reactions in %s', label), two, reasonTwo, 'conditional');
end

% -- metabolite turnover ----------------------------------------------
for k = 1:numel(opt.turnoverMetabolites)
    met = opt.turnoverMetabolites{k};
    row = find(strcmp(model.mets, met), 1);
    present = ~isempty(row);
    reason = '';
    if ~present, reason = sprintf('%s is not in this model', met); end
    weights = [];
    if present, weights = full(absS(row, internal))'; end
    d(end+1) = mk(['turnover_' met], 'turnover', 'internal', 'flux', ...
        sprintf(['gross stoichiometric throughput of %s; not a measurable ' ...
                 'turnover'], met), ...
        @(V, vals) abs(V(:, internal)) * weights, ...
        sprintf('the metabolite %s', met), present, reason, 'conditional');
end

% -- named reactions ---------------------------------------------------
for k = 1:size(opt.namedReactions, 1)
    label = opt.namedReactions{k, 1};
    identifier = opt.namedReactions{k, 2};
    column = find(strcmp(model.rxns, identifier), 1);
    present = ~isempty(column);
    reason = '';
    if ~present, reason = sprintf('%s is not in this model', identifier); end
    d(end+1) = mk(['rxn_' label], 'biochemical', 'internal', 'flux', ...
        sprintf('flux through %s, in the model''s own sign convention', identifier), ...
        @(V, vals) V(:, column), ...
        sprintf('the reaction %s', identifier), present, reason, 'conditional');
end

% -- efficiency ratios -------------------------------------------------
for k = 1:size(opt.efficiencyRatios, 1)
    label = opt.efficiencyRatios{k, 1};
    numerator = opt.efficiencyRatios{k, 2};
    denominator = opt.efficiencyRatios{k, 3};
    haveNumerator = isComputable(d, numerator);
    haveDenominator = isComputable(d, denominator);
    present = haveNumerator && haveDenominator;
    reason = '';
    if ~present
        missing = {};
        if ~haveNumerator, missing{end+1} = numerator; end
        if ~haveDenominator, missing{end+1} = denominator; end
        reason = ['needs ' strjoin(missing, ', ')];
    end
    d(end+1) = mk(['eff_' label], 'efficiency', 'internal', 'ratio', ...
        sprintf('%s relative to %s', numerator, denominator), ...
        @(V, vals) safeDivide(abs(vals.(numerator)), abs(vals.(denominator))), ...
        sprintf('%s and %s', numerator, denominator), present, reason, 'conditional');
end

% -- cofactor allocation ----------------------------------------------
for k = 1:size(opt.cofactorPairs, 1)
    label = opt.cofactorPairs{k, 1};
    pair = opt.cofactorPairs{k, 2};
    rowA = find(strcmp(model.mets, pair{1}), 1);
    rowB = find(strcmp(model.mets, pair{2}), 1);
    present = ~isempty(rowA) && ~isempty(rowB);
    reason = '';
    if ~present, reason = 'one member of the pair is absent'; end
    wA = []; wB = [];
    if present
        wA = full(absS(rowA, internal))';
        wB = full(absS(rowB, internal))';
    end
    d(end+1) = mk(['cof_' label '_share'], 'cofactor', 'internal', 'ratio', ...
        sprintf('share of the %s/%s pool turnover carried by %s', ...
                pair{1}, pair{2}, pair{1}), ...
        @(V, vals) safeDivide(abs(V(:, internal)) * wA, ...
                              abs(V(:, internal)) * wA + abs(V(:, internal)) * wB), ...
        sprintf('both %s and %s', pair{1}, pair{2}), present, reason, 'conditional');
end

% -- compartments: where in the cell a metabolite is made --------------
% Production of x in compartment K, booked by the sign of each flux, over the
% reactions that make or use x chemically: a reaction in which x appears in
% two compartments only moves it. Structurally unavailable, with the reason,
% unless x lives in two intracellular compartments (shares) or the model has
% two (flux shares). Identifiers carry the compartment as atp_m or atp[m].
[metBase, metComp] = splitMetabolites(model.mets);
outside = (ismember(metComp, opt.extracellular) | cellfun(@isempty, metComp))';
inside = reshape(unique(metComp(~outside)), 1, []);
internalColumns = find(internal);
placed = absS(:, internalColumns) > 0;
home = repmat({''}, 1, numel(internalColumns));
between = false(1, numel(internalColumns));
for j = 1:numel(internalColumns)
    places = unique(metComp(placed(:, j) & ~outside));
    if numel(places) == 1, home{j} = places{1}; end
    between(j) = numel(places) >= 2;
end
for k = 1:numel(opt.compartmentMetabolites)
    met = opt.compartmentMetabolites{k};
    rows = find(strcmp(metBase(:), met) & ~outside);
    places = reshape(unique(metComp(rows)), 1, []);
    for c = 1:size(opt.compartmentNames, 1)
        code = opt.compartmentNames{c, 1};
        place = opt.compartmentNames{c, 2};
        present = any(strcmp(places, code)) && numel(places) >= 2;
        reason = '';
        if ~present
            reason = sprintf(['%s is in %d intracellular compartment(s) (%s); ' ...
                              'needs %s and one other'], met, numel(places), ...
                              strjoin(places, ', '), code);
        end
        d(end+1) = mk(sprintf('cmp_%s_made_in_%s', met, place), 'compartment', ...
            'internal', 'ratio', ...
            sprintf('share of the cell''s gross %s production made in the %s (stoichiometric, not net)', upper(met), place), ...
            @(V, vals) compartmentShare(V, model.S, rows, metComp, internalColumns, code), ...
            sprintf('%s in the %s and in another intracellular compartment', met, place), ...
            present, reason, 'conditional');
    end
end
two = numel(inside) >= 2;
why = '';
if ~two
    why = sprintf('%d intracellular compartment(s) (%s); needs two', numel(inside), ...
                  strjoin(inside, ', '));
end
for c = 1:size(opt.compartmentNames, 1)
    code = opt.compartmentNames{c, 1};
    place = opt.compartmentNames{c, 2};
    present = two && any(strcmp(inside, code));
    reason = why;
    if two && ~present, reason = sprintf('no %s in this model', place); end
    columns = internalColumns(strcmp(home, code));
    d(end+1) = mk(['cmp_flux_in_' place], 'compartment', 'internal', 'ratio', ...
        sprintf('share of internal flux carried by reactions inside the %s', place), ...
        @(V, vals) safeDivide(sum(abs(V(:, columns)), 2), sum(abs(V(:, internal)), 2)), ...
        'at least two intracellular compartments', present, reason, 'conditional');
end
columns = internalColumns(between);
d(end+1) = mk('cmp_flux_between_compartments', 'compartment', 'internal', 'ratio', ...
    'share of internal flux carried by reactions joining two intracellular compartments', ...
    @(V, vals) safeDivide(sum(abs(V(:, columns)), 2), sum(abs(V(:, internal)), 2)), ...
    'at least two intracellular compartments', two, why, 'conditional');

% -- elements, amino acids, oxidative defence, budgets and conversion ---
% The export of an exchanged metabolite is read from the sign of its
% coefficient, so it is positive out of the cell whichever way the exchange is
% written. Shares and ratios are NaN where their denominator is zero.
exCols = [];
exRows = [];
for c = find(external)
    rows = find(abs(model.S(:, c)) > 0);
    if numel(rows) == 1
        exCols(end+1) = c;
        exRows(end+1) = rows;
    end
end
exSign = zeros(1, 0);
if ~isempty(exCols)
    exSign = -sign(full(model.S(sub2ind(size(model.S), exRows, exCols))));
end
exStem = metBase(exRows);
atoms = struct();
elements = {'C', 'N', 'S', 'P'};
for e = 1:numel(elements)
    counts = zeros(1, numel(exRows));
    for k = 1:numel(exRows)
        counts(k) = elementCount(model.metFormulas{exRows(k)}, elements{e});
    end
    atoms.(elements{e}) = counts;
end

elementNames = {'N', 'nitrogen'; 'S', 'sulfur'; 'P', 'phosphorus'};
sides = {'uptake', 'taken up', true; 'secretion', 'secreted', false};
for e = 1:size(elementNames, 1)
    symbol = elementNames{e, 1};
    label = elementNames{e, 2};
    weights = atoms.(symbol);
    present = any(weights > 0);
    reason = '';
    if ~present, reason = sprintf('no exchanged metabolite contains %s', label); end
    for s = 1:size(sides, 1)
        d(end+1) = mk(sprintf('elem_%s_%s', lower(symbol), sides{s, 1}), 'elemental', ...
            'external', 'flux', ...
            sprintf('%s atoms %s per unit time, from the formulas of exchanged metabolites', ...
                    label, sides{s, 2}), ...
            @(V, vals) sum(exchangeParts(V, exCols, exSign, weights, sides{s, 3}), 2), ...
            sprintf('an exchanged metabolite with %s in its formula', label), ...
            present, reason, 'conditional');
    end
end

hasCarbon = any(atoms.C > 0);
hasNitrogen = any(atoms.N > 0);
amino = cell(1, numel(exRows));
for k = 1:numel(exRows)
    amino{k} = aminoAcidOf(model.mets{exRows(k)}, model.metNames{exRows(k)});
end
aaCols = find(~cellfun(@isempty, amino));
whereStem = @(stems) find(ismember(exStem, stems));
shares = { ...
    'elem_c_from_glucose', 'share of carbon uptake supplied by glucose', ...
        whereStem({'glc__D', 'glc_D'}), true, 'C', 'a glucose exchange'; ...
    'elem_c_from_lactate_pyruvate', 'share of carbon uptake supplied by lactate and pyruvate', ...
        whereStem({'lac__L', 'lac_L', 'pyr'}), true, 'C', 'a lactate or pyruvate exchange'; ...
    'elem_c_from_amino_acids', 'share of carbon uptake supplied by amino acids', ...
        aaCols, true, 'C', 'an amino-acid exchange'; ...
    'elem_c_to_lactate', 'share of carbon secretion leaving as lactate', ...
        whereStem({'lac__L', 'lac_L'}), false, 'C', 'a lactate exchange'; ...
    'elem_c_to_co2', 'share of carbon secretion leaving as CO2', ...
        whereStem({'co2'}), false, 'C', 'a CO2 exchange'; ...
    'elem_n_from_amino_acids', 'share of nitrogen uptake supplied by amino acids', ...
        aaCols, true, 'N', 'an amino-acid exchange'};
for k = 1:size(shares, 1)
    columns = shares{k, 3};
    uptake = shares{k, 4};
    weights = atoms.(shares{k, 5});
    present = any(weights > 0) && ~isempty(columns);
    reason = '';
    if ~present, reason = ['no ' regexprep(shares{k, 6}, '^an? ', '')]; end
    d(end+1) = mk(shares{k, 1}, 'elemental', 'external', 'ratio', shares{k, 2}, ...
        @(V, vals) partShare(exchangeParts(V, exCols, exSign, weights, uptake), columns), ...
        shares{k, 6}, present, reason, 'conditional');
end

diversityKinds = { ...
    'c', 'carbon', 'source_entropy', 'entropy', 'nats', 'Shannon entropy of the carbon sources taken up'; ...
    'c', 'carbon', 'effective_sources', 'effective', 'dimensionless', 'effective number of carbon sources'; ...
    'c', 'carbon', 'dominant_source', 'dominant', 'ratio', 'share of carbon uptake from the largest single source'; ...
    'n', 'nitrogen', 'source_entropy', 'entropy', 'nats', 'Shannon entropy of the nitrogen sources taken up'; ...
    'n', 'nitrogen', 'effective_sources', 'effective', 'dimensionless', 'effective number of nitrogen sources'};
for k = 1:size(diversityKinds, 1)
    weights = atoms.(upper(diversityKinds{k, 1}));
    present = any(weights > 0);
    reason = '';
    if ~present, reason = sprintf('no exchanged metabolite contains %s', diversityKinds{k, 2}); end
    kind = diversityKinds{k, 4};
    d(end+1) = mk(sprintf('div_%s_%s', diversityKinds{k, 1}, diversityKinds{k, 3}), ...
        'source_diversity', 'external', diversityKinds{k, 5}, diversityKinds{k, 6}, ...
        @(V, vals) diversityOf(exchangeParts(V, exCols, exSign, weights, true), kind), ...
        sprintf('an exchanged metabolite with %s in its formula', diversityKinds{k, 2}), ...
        present, reason, 'conditional');
end

[acidTable, essentialAcids, trackedAcids] = aminoAcidTable();
acidNames = unique(amino(aaCols), 'stable');
byAcid = cell(1, numel(acidNames));
for k = 1:numel(acidNames)
    byAcid{k} = aaCols(strcmp(amino(aaCols), acidNames{k}));
end
for k = 1:size(acidTable, 1)
    acid = acidTable{k, 1};
    if any(strcmp(trackedAcids, acid)), continue, end
    at = find(strcmp(acidNames, acid), 1);
    present = ~isempty(at);
    columns = [];
    reason = '';
    if present, columns = byAcid{at}; else, reason = sprintf('no %s exchange', acid); end
    d(end+1) = mk(sprintf('aa_%s_net_uptake', acid), 'amino_acid', 'external', 'flux', ...
        sprintf('net uptake of %s (uptake minus secretion)', acid), ...
        @(V, vals) -sum(exportOf(V, exCols, exSign, columns), 2), ...
        sprintf('an exchange for %s', acid), present, reason, 'conditional');
end
found = ~isempty(acidNames);
essential = find(ismember(acidNames, essentialAcids));
aggregates = { ...
    'aa_total_uptake', 'flux', 'total uptake of amino acids', 'uptake'; ...
    'aa_total_secretion', 'flux', 'total secretion of amino acids', 'secretion'; ...
    'aa_essential_uptake', 'flux', 'uptake of the nine essential amino acids', 'essential'; ...
    'aa_uptake_entropy', 'nats', 'Shannon entropy of the amino-acid uptake profile', 'entropy'; ...
    'aa_n_taken_up', 'count', 'number of amino acids taken up', 'count'};
for k = 1:size(aggregates, 1)
    present = found && (~strcmp(aggregates{k, 1}, 'aa_essential_uptake') || ~isempty(essential));
    reason = '';
    if ~present, reason = 'no recognised amino-acid exchange'; end
    kind = aggregates{k, 4};
    d(end+1) = mk(aggregates{k, 1}, 'amino_acid', 'external', aggregates{k, 2}, aggregates{k, 3}, ...
        @(V, vals) aminoAggregate(acidMatrix(V, exCols, exSign, byAcid), kind, essential, tol), ...
        'an amino-acid exchange recognised by identifier or name', present, reason, 'conditional');
end

% carriers: the rows of every intracellular compartment, moves excluded
intCols = find(internal);
Sint = model.S(:, intCols);
h2o2 = carrierCoefficient(Sint, metBase, metComp, {'h2o2'}, opt.extracellular);
o2s = carrierCoefficient(Sint, metBase, metComp, {'o2s'}, opt.extracellular);
nadph = carrierCoefficient(Sint, metBase, metComp, {'nadph'}, opt.extracellular);
atp = carrierCoefficient(Sint, metBase, metComp, {'atp'}, opt.extracellular);
live = @(c) ~isempty(c) && any(c);
catalase = reactionSignature(Sint, metBase, {'h2o2'}, {'o2'}, true, '');
gpx = reactionSignature(Sint, metBase, {'h2o2', 'gthrd'}, {'gthox'}, false, '');
prdx = reactionSignature(Sint, metBase, {'h2o2'}, {}, false, 'prdx');
sod = reactionSignature(Sint, metBase, {'o2s'}, {'h2o2'}, false, '');
grNadph = reactionSignature(Sint, metBase, {'gthox', 'nadph'}, {'gthrd'}, false, '');
gr = union(grNadph, reactionSignature(Sint, metBase, {'gthox', 'nadh'}, {'gthrd'}, false, ''));
defence = { ...
    'ros_h2o2_production', 'flux', 'gross internal production of hydrogen peroxide', ...
        live(h2o2), 'no peroxide', ...
        @(V, vals) sum(madeBy(V, intCols, h2o2), 2); ...
    'ros_h2o2_removal', 'flux', 'gross internal removal of hydrogen peroxide', ...
        live(h2o2), 'no peroxide', ...
        @(V, vals) sum(usedBy(V, intCols, h2o2), 2); ...
    'ros_h2o2_removal_per_production', 'ratio', ...
        'model-based peroxide removal-to-production ratio (above 1 when peroxide is taken up)', ...
        live(h2o2), 'no peroxide', ...
        @(V, vals) safeDivide(sum(usedBy(V, intCols, h2o2), 2), sum(madeBy(V, intCols, h2o2), 2)); ...
    'ros_h2o2_removed_by_catalase', 'ratio', 'share of peroxide removal by catalase', ...
        live(h2o2) && ~isempty(catalase), 'no catalase signature (2 H2O2 -> O2)', ...
        @(V, vals) partShare(usedBy(V, intCols, h2o2), catalase); ...
    'ros_h2o2_removed_by_glutathione_peroxidase', 'ratio', ...
        'share of peroxide removal by glutathione peroxidase', ...
        live(h2o2) && ~isempty(gpx), 'no glutathione peroxidase signature', ...
        @(V, vals) partShare(usedBy(V, intCols, h2o2), gpx); ...
    'ros_h2o2_removed_by_peroxiredoxin', 'ratio', 'share of peroxide removal by peroxiredoxins', ...
        live(h2o2) && ~isempty(prdx), 'no peroxiredoxin signature', ...
        @(V, vals) partShare(usedBy(V, intCols, h2o2), prdx); ...
    'ros_superoxide_production', 'flux', 'gross internal production of superoxide', ...
        live(o2s), 'no superoxide', ...
        @(V, vals) sum(madeBy(V, intCols, o2s), 2); ...
    'ros_superoxide_dismutated', 'ratio', 'share of superoxide removal by superoxide dismutase', ...
        live(o2s) && ~isempty(sod), 'no superoxide dismutase signature', ...
        @(V, vals) partShare(usedBy(V, intCols, o2s), sod); ...
    'ros_glutathione_regeneration_cover', 'ratio', ...
        'NADPH produced per NADPH consumed by glutathione reductases', ...
        live(nadph) && ~isempty(grNadph), 'no NADPH-dependent glutathione reductase', ...
        @(V, vals) safeDivide(sum(madeBy(V, intCols, nadph), 2), ...
                              sumColumns(usedBy(V, intCols, nadph), grNadph)); ...
    'ros_glutathione_cycle_closure', 'ratio', ...
        'glutathione disulfide reduced by glutathione reductases per disulfide formed by peroxide removal', ...
        ~isempty(gr) && ~isempty(gpx), 'no glutathione reductase or peroxidase', ...
        @(V, vals) safeDivide(sum(max(V(:, intCols(gr)), 0), 2), ...
                              sum(max(V(:, intCols(gpx)), 0), 2))};
for k = 1:size(defence, 1)
    reason = '';
    if ~defence{k, 4}, reason = defence{k, 5}; end
    d(end+1) = mk(defence{k, 1}, 'oxidative_defense', 'internal', defence{k, 2}, defence{k, 3}, ...
        defence{k, 6}, 'peroxide or superoxide and the named signature', ...
        defence{k, 4}, reason, 'conditional');
end

% budget: NADPH consumed by glutathione reductases (the signature above), and
% NADPH made by the oxidative branch (the named G6PDH2 and GND) per NADPH so used
oxidative = [];
for label = {'g6pd', 'gnd'}
    at = find(strcmp(opt.namedReactions(:, 1), label{1}), 1);
    if isempty(at), continue, end
    column = find(strcmp(model.rxns, opt.namedReactions{at, 2}), 1);
    position = find(intCols == column, 1);
    if ~isempty(position), oxidative(end+1) = position; end
end
oxidative = sort(oxidative);
branch = numel(oxidative) == 2;
budget = { ...
    'bud_nadph_to_glutathione_reduction', ...
        'share of NADPH consumption used to reduce glutathione disulfide', ...
        live(nadph) && ~isempty(gr), 'no NADPH or no glutathione reductase signature', ...
        @(V, vals) partShare(usedBy(V, intCols, nadph), gr); ...
    'bud_oxppp_cover_of_glutathione_reduction', ...
        'NADPH made by the oxidative branch per NADPH used to reduce glutathione', ...
        live(nadph) && ~isempty(gr) && branch, ...
        'no NADPH, no glutathione reductase signature, or G6PDH2 and GND not both present', ...
        @(V, vals) safeDivide(sumColumns(madeBy(V, intCols, nadph), oxidative), ...
                              sumColumns(usedBy(V, intCols, nadph), gr))};
for k = 1:size(budget, 1)
    reason = '';
    if ~budget{k, 3}, reason = budget{k, 4}; end
    d(end+1) = mk(budget{k, 1}, 'budget', 'internal', 'ratio', budget{k, 2}, budget{k, 5}, ...
        'NADPH, the glutathione-reductase signature and the named reactions', ...
        budget{k, 3}, reason, 'conditional');
end

glucose = whereStem({'glc__D', 'glc_D'});
lactate = whereStem({'lac__L', 'lac_L'});
co2 = whereStem({'co2'});
carbon = atoms.C;
conversion = { ...
    'conv_atp_per_carbon', 'gross ATP production per carbon atom taken up', ...
        hasCarbon && ~isempty(atp), 'no ATP or no carbon exchange', ...
        @(V, vals) safeDivide(sum(madeBy(V, intCols, atp), 2), ...
                              sum(exchangeParts(V, exCols, exSign, carbon, true), 2)); ...
    'conv_nadph_per_carbon', 'gross NADPH production per carbon atom taken up', ...
        hasCarbon && ~isempty(nadph), 'no NADPH or no carbon exchange', ...
        @(V, vals) safeDivide(sum(madeBy(V, intCols, nadph), 2), ...
                              sum(exchangeParts(V, exCols, exSign, carbon, true), 2)); ...
    'conv_lactate_per_glucose_carbon', 'lactate carbon secreted per glucose carbon taken up', ...
        hasCarbon && ~isempty(glucose) && ~isempty(lactate), ...
        'no glucose or lactate exchange, or no formulas with carbon', ...
        @(V, vals) safeDivide(sumColumns(exchangeParts(V, exCols, exSign, carbon, false), lactate), ...
                              sumColumns(exchangeParts(V, exCols, exSign, carbon, true), glucose)); ...
    'conv_co2_per_carbon', 'CO2 carbon released per carbon atom taken up', ...
        hasCarbon && ~isempty(co2), 'no CO2 exchange or no carbon exchange', ...
        @(V, vals) safeDivide(sumColumns(exchangeParts(V, exCols, exSign, carbon, false), co2), ...
                              sum(exchangeParts(V, exCols, exSign, carbon, true), 2))};
for k = 1:size(conversion, 1)
    reason = '';
    if ~conversion{k, 3}, reason = conversion{k, 4}; end
    d(end+1) = mk(conversion{k, 1}, 'conversion', 'combined', 'ratio', conversion{k, 2}, ...
        conversion{k, 5}, 'metabolite formulas, the carrier and the named exchanges', ...
        conversion{k, 3}, reason, 'conditional');
end

% -- capacity: the only family that depends on the constraints ---------
finite = isfinite(model.lb)' & isfinite(model.ub)' & (model.ub' > model.lb');
for b = 1:size(blocks, 1)
    block = blocks{b, 1}; mask = blocks{b, 2}; tag = blocks{b, 3};
    usable = finite & mask;
    present = any(usable);
    reason = '';
    if ~present, reason = sprintf('no %s reaction has finite bounds', block); end
    lo = model.lb(usable)'; hi = model.ub(usable)';
    d(end+1) = mk(['cap_' tag '_mean_saturation'], 'capacity', block, 'ratio', ...
        sprintf(['mean position within the admissible interval, over bounded ' ...
                 '%s reactions'], block), ...
        @(V, vals) mean((V(:, usable) - lo) ./ (hi - lo), 2), ...
        'finite bounds', present, reason, 'model_specific');
    d(end+1) = mk(['cap_' tag '_at_bound'], 'capacity', block, 'ratio', ...
        sprintf(['fraction of bounded %s reactions within tolerance of a bound'], block), ...
        @(V, vals) mean(V(:, usable) <= lo + tol | V(:, usable) >= hi - tol, 2), ...
        'finite bounds', present, reason, 'model_specific');
end

% -- the boundary ------------------------------------------------------
if any(external)
    d(end+1) = mk('ex_total_uptake', 'exchange_burden', 'external', 'flux', ...
        'total magnitude taken up across the boundary', ...
        @(V, vals) sum(max(-V(:, external), 0), 2), 'an exchange reaction');
    d(end+1) = mk('ex_total_secretion', 'exchange_burden', 'external', 'flux', ...
        'total magnitude secreted across the boundary', ...
        @(V, vals) sum(max(V(:, external), 0), 2), 'an exchange reaction');
    d(end+1) = mk('ex_n_uptaken', 'exchange_burden', 'external', 'count', ...
        'number of metabolites taken up', ...
        @(V, vals) sum(V(:, external) < -tol, 2), 'an exchange reaction');
    d(end+1) = mk('ex_n_secreted', 'exchange_burden', 'external', 'count', ...
        'number of metabolites secreted', ...
        @(V, vals) sum(V(:, external) > tol, 2), 'an exchange reaction');
end

boundaryOf = containers.Map('KeyType', 'char', 'ValueType', 'double');
exColumns = find(external);
for k = 1:numel(exColumns)
    rows = find(abs(model.S(:, exColumns(k))) > 0);
    if numel(rows) == 1
        boundaryOf(model.mets{rows}) = exColumns(k);
    end
end
for k = 1:size(opt.trackedExchanges, 1)
    met = opt.trackedExchanges{k, 1};
    label = opt.trackedExchanges{k, 2};
    present = isKey(boundaryOf, met);
    column = [];
    reason = '';
    if present, column = boundaryOf(met);
    else, reason = sprintf('no exchange reaction for %s', met); end
    d(end+1) = mk(['ex_' label '_net'], 'boundary_metabolite', 'external', 'flux', ...
        sprintf('net exchange of %s; negative is uptake', label), ...
        @(V, vals) V(:, column), ...
        sprintf('an exchange for %s', met), present, reason, 'conditional');
end

% -- coupling: only meaningful because the blocks are kept apart -------
both = any(internal) && any(external);
reasonBoth = '';
if ~both, reasonBoth = 'one of the two blocks is empty'; end
d(end+1) = mk('cpl_boundary_share', 'coupling', 'combined', 'ratio', ...
    'share of the state''s total flux magnitude that crosses the boundary', ...
    @(V, vals) safeDivide(sum(abs(V(:, external)), 2), sum(abs(V), 2)), ...
    'both blocks non-empty', both, reasonBoth, 'conditional');
d(end+1) = mk('cpl_internal_per_boundary', 'coupling', 'combined', 'ratio', ...
    'internal activity supported per unit of boundary exchange', ...
    @(V, vals) safeDivide(sum(abs(V(:, internal)), 2), sum(abs(V(:, external)), 2)), ...
    'both blocks non-empty', both, reasonBoth, 'conditional');
reasonEx = '';
if ~any(external), reasonEx = 'no exchange reactions'; end
d(end+1) = mk('cpl_secretion_over_uptake', 'coupling', 'combined', 'ratio', ...
    'secreted magnitude relative to taken-up magnitude', ...
    @(V, vals) safeDivide(sum(max(V(:, external), 0), 2), ...
                          sum(max(-V(:, external), 0), 2)), ...
    'an exchange reaction', any(external), reasonEx, 'conditional');
end


function [base, comp] = splitMetabolites(mets)
% atp_m or atp[m] -> 'atp' and 'm'; no recognisable suffix -> the id and ''.
base = mets(:)';
comp = repmat({''}, 1, numel(mets));
for k = 1:numel(mets)
    tokens = regexp(mets{k}, '^(.*?)(?:_([A-Za-z0-9]{1,2})|\[([A-Za-z0-9]{1,2})\])$', ...
                    'tokens', 'once');
    if isempty(tokens), continue, end
    base{k} = tokens{1};
    if numel(tokens) >= 2 && ~isempty(tokens{2})
        comp{k} = tokens{2};
    elseif numel(tokens) >= 3
        comp{k} = tokens{3};
    end
end
end


function share = compartmentShare(V, S, rows, metComp, internalColumns, code)
% Chemical production of a metabolite in compartment CODE over its production
% in every intracellular compartment; transport of the metabolite is excluded.
Sm = S(rows, internalColumns);
moving = full(sum(abs(Sm) > 0, 1)) >= 2;
total = zeros(size(V, 1), 1);
mine = zeros(size(V, 1), 1);
for r = 1:numel(rows)
    weights = full(Sm(r, :)) .* ~moving;
    made = sum(max(V(:, internalColumns) .* weights, 0), 2);
    total = total + made;
    if strcmp(metComp{rows(r)}, code), mine = mine + made; end
end
share = safeDivide(mine, total);
end


function tf = isComputable(d, name)
% Has this descriptor been declared already, and can this model compute it?
tf = false;
index = find(strcmp({d.name}, name), 1);
if ~isempty(index), tf = d(index).computable; end
end


function d = emptyDescriptor()
d = struct('name', {}, 'family', {}, 'block', {}, 'units', {}, ...
           'description', {}, 'requires', {}, 'computable', {}, 'reason', {}, ...
           'scope', {}, 'compute', {});
end


function s = mk(name, family, block, units, description, compute, requires, ...
                computable, reason, scope)
if nargin < 7,  requires = 'nothing named'; end
if nargin < 8,  computable = true; end
if nargin < 9,  reason = ''; end
if nargin < 10, scope = 'common'; end
s = struct('name', name, 'family', family, 'block', block, 'units', units, ...
           'description', description, 'requires', requires, ...
           'computable', logical(computable), 'reason', reason, ...
           'scope', scope, 'compute', compute);
end


function keep = selectDescriptors(declared, which)
% Which descriptors to return: tiers of the release, 'all', 'computable' or names.
tierNames = {'core', 'module', 'generic', 'alias'};
if isstring(which)
    which = cellstr(which);
    if numel(which) == 1
        which = which{1};
    end
end
if ischar(which) && any(strcmpi(which, tierNames))
    which = {lower(which)};
end
if iscell(which) && ~isempty(which) && all(ismember(lower(which), tierNames))
    if all(strcmp({declared.tier}, 'unclassified'))
        error('flux2feature:status', ['this copy of flux2feature carries no descriptor ' ...
            'status; use ''Which'', ''all'' or a list of names']);
    end
    keep = ismember({declared.tier}, lower(which));
elseif iscell(which)
    keep = ismember({declared.name}, which);
    unknown = setdiff(which, {declared.name});
    if ~isempty(unknown)
        error('flux2feature:unknown', 'not declared: %s', strjoin(unknown, ', '));
    end
elseif strcmpi(which, 'computable')
    keep = [declared.computable];
elseif strcmpi(which, 'all')
    keep = true(1, numel(declared));
else
    error('flux2feature:which', ['Which must be a tier (''core'', ''module'', ' ...
        '''generic'', ''alias''), a cell array of tiers, ''computable'', ''all'' ' ...
        'or a cell array of names']);
end
end


function declared = attachStatus(declared)
% Gives every declared descriptor its tier of the release and, for an alias,
% the canonical descriptor it repeats (see descriptorStatus).
%
% USAGE:
%
%    declared = attachStatus(declared)
status = descriptorStatus();
tiers = repmat({'unclassified'}, 1, numel(declared));
canonical = repmat({''}, 1, numel(declared));
if ~isempty(status)
    [found, at] = ismember({declared.name}, status(:, 1)');
    tiers(found) = status(at(found), 2)';
    canonical(found) = status(at(found), 3)';
end
[declared.tier] = tiers{:};
[declared.canonical] = canonical{:};
end


function status = descriptorStatus()
% The tier of every declared descriptor in the release and, for an alias,
% its canonical descriptor: one row {name, tier, canonical} per descriptor.
% Written by the analysis notebook (stage 28) from the audit of the study;
% do not edit by hand.
% BEGIN DESCRIPTOR STATUS
status = { ...
    'int_L1', 'generic', ''; ...
    'int_L2', 'generic', ''; ...
    'int_Linf', 'generic', ''; ...
    'int_active_fraction', 'generic', ''; ...
    'int_gini', 'alias', 'int_shannon_entropy'; ...
    'int_shannon_entropy', 'generic', ''; ...
    'int_directionality', 'generic', ''; ...
    'int_cv_abs', 'generic', ''; ...
    'int_top10pct_share', 'generic', ''; ...
    'ext_L1', 'generic', ''; ...
    'ext_L2', 'generic', ''; ...
    'ext_Linf', 'generic', ''; ...
    'ext_active_fraction', 'generic', ''; ...
    'ext_gini', 'alias', 'ext_cv_abs'; ...
    'ext_shannon_entropy', 'alias', 'ext_cv_abs'; ...
    'ext_directionality', 'generic', ''; ...
    'ext_cv_abs', 'generic', ''; ...
    'ext_top10pct_share', 'generic', ''; ...
    'net_active_metabolites', 'core', ''; ...
    'net_top_decile_turnover_share', 'core', ''; ...
    'path_glycolysis_activity', 'core', ''; ...
    'path_glycolysis_share', 'core', ''; ...
    'path_glycolysis_engagement', 'module', ''; ...
    'path_glycolysis_directionality', 'core', ''; ...
    'path_pentose_phosphate_activity', 'core', ''; ...
    'path_pentose_phosphate_share', 'core', ''; ...
    'path_pentose_phosphate_engagement', 'module', ''; ...
    'path_pentose_phosphate_directionality', 'core', ''; ...
    'path_glutathione_activity', 'core', ''; ...
    'path_glutathione_share', 'alias', 'path_glutathione_activity'; ...
    'path_glutathione_engagement', 'module', ''; ...
    'path_glutathione_directionality', 'core', ''; ...
    'path_redox_detox_activity', 'core', ''; ...
    'path_redox_detox_share', 'alias', 'path_redox_detox_activity'; ...
    'path_redox_detox_engagement', 'module', ''; ...
    'path_redox_detox_directionality', 'core', ''; ...
    'path_purine_activity', 'core', ''; ...
    'path_purine_share', 'alias', 'path_purine_activity'; ...
    'path_purine_engagement', 'module', ''; ...
    'path_purine_directionality', 'core', ''; ...
    'path_pyrimidine_activity', 'core', ''; ...
    'path_pyrimidine_share', 'alias', 'path_pyrimidine_activity'; ...
    'path_pyrimidine_engagement', 'module', ''; ...
    'path_pyrimidine_directionality', 'core', ''; ...
    'path_transport_activity', 'core', ''; ...
    'path_transport_share', 'core', ''; ...
    'path_transport_engagement', 'module', ''; ...
    'path_transport_directionality', 'core', ''; ...
    'path_rapoport_luebering_activity', 'core', ''; ...
    'path_rapoport_luebering_share', 'alias', 'path_rapoport_luebering_activity'; ...
    'path_rapoport_luebering_engagement', 'module', ''; ...
    'path_rapoport_luebering_directionality', 'core', ''; ...
    'turnover_atp_c', 'core', ''; ...
    'turnover_adp_c', 'core', ''; ...
    'turnover_amp_c', 'core', ''; ...
    'turnover_nadh_c', 'core', ''; ...
    'turnover_nad_c', 'core', ''; ...
    'turnover_nadph_c', 'core', ''; ...
    'turnover_nadp_c', 'alias', 'turnover_nadph_c'; ...
    'turnover_gthrd_c', 'core', ''; ...
    'turnover_gthox_c', 'core', ''; ...
    'turnover_glc__D_c', 'core', ''; ...
    'turnover_g6p_c', 'core', ''; ...
    'turnover_f6p_c', 'core', ''; ...
    'turnover_pyr_c', 'core', ''; ...
    'turnover_lac__L_c', 'core', ''; ...
    'turnover_h2o2_c', 'core', ''; ...
    'turnover_o2_c', 'core', ''; ...
    'turnover_13dpg_c', 'core', ''; ...
    'turnover_r5p_c', 'core', ''; ...
    'turnover_ru5p__D_c', 'alias', 'turnover_f6p_c'; ...
    'turnover_6pgc_c', 'core', ''; ...
    'rxn_g6pd', 'alias', 'turnover_6pgc_c'; ...
    'rxn_pgl', 'alias', 'turnover_6pgc_c'; ...
    'rxn_gnd', 'alias', 'turnover_6pgc_c'; ...
    'rxn_catalase', 'alias', 'ros_h2o2_removed_by_catalase'; ...
    'rxn_gpx', 'core', ''; ...
    'rxn_sod', 'core', ''; ...
    'rxn_hexokinase', 'core', ''; ...
    'rxn_pfk', 'core', ''; ...
    'rxn_pyruvate_kinase', 'core', ''; ...
    'rxn_ldh', 'core', ''; ...
    'rxn_adenylate_kinase', 'core', ''; ...
    'eff_nadph_per_glucose', 'core', ''; ...
    'eff_nadph_per_atp', 'core', ''; ...
    'eff_atp_per_glucose', 'core', ''; ...
    'eff_nadph_per_h2o2', 'core', ''; ...
    'eff_redox_buffer_ratio', 'alias', 'cof_glutathione_share'; ...
    'eff_ppp_per_glycolysis', 'core', ''; ...
    'eff_ppp_flux_split', 'core', ''; ...
    'cof_nadph_nadp_share', 'core', ''; ...
    'cof_nadh_nad_share', 'core', ''; ...
    'cof_glutathione_share', 'core', ''; ...
    'cmp_atp_made_in_mitochondrion', 'core', ''; ...
    'cmp_atp_made_in_cytosol', 'core', ''; ...
    'cmp_nadh_made_in_mitochondrion', 'core', ''; ...
    'cmp_nadh_made_in_cytosol', 'core', ''; ...
    'cmp_nadph_made_in_mitochondrion', 'core', ''; ...
    'cmp_nadph_made_in_cytosol', 'core', ''; ...
    'cmp_co2_made_in_mitochondrion', 'core', ''; ...
    'cmp_co2_made_in_cytosol', 'core', ''; ...
    'cmp_flux_in_mitochondrion', 'core', ''; ...
    'cmp_flux_in_cytosol', 'core', ''; ...
    'cmp_flux_between_compartments', 'core', ''; ...
    'elem_n_uptake', 'module', ''; ...
    'elem_n_secretion', 'alias', 'elem_n_uptake'; ...
    'elem_s_uptake', 'module', ''; ...
    'elem_s_secretion', 'alias', 'elem_s_uptake'; ...
    'elem_p_uptake', 'module', ''; ...
    'elem_p_secretion', 'module', ''; ...
    'elem_c_from_glucose', 'module', ''; ...
    'elem_c_from_lactate_pyruvate', 'module', ''; ...
    'elem_c_from_amino_acids', 'module', ''; ...
    'elem_c_to_lactate', 'alias', 'rxn_ldh'; ...
    'elem_c_to_co2', 'module', ''; ...
    'elem_n_from_amino_acids', 'module', ''; ...
    'div_c_source_entropy', 'module', ''; ...
    'div_c_effective_sources', 'module', ''; ...
    'div_c_dominant_source', 'module', ''; ...
    'div_n_source_entropy', 'module', ''; ...
    'div_n_effective_sources', 'module', ''; ...
    'aa_alanine_net_uptake', 'module', ''; ...
    'aa_arginine_net_uptake', 'module', ''; ...
    'aa_asparagine_net_uptake', 'module', ''; ...
    'aa_aspartate_net_uptake', 'module', ''; ...
    'aa_glutamine_net_uptake', 'module', ''; ...
    'aa_histidine_net_uptake', 'module', ''; ...
    'aa_isoleucine_net_uptake', 'module', ''; ...
    'aa_leucine_net_uptake', 'module', ''; ...
    'aa_lysine_net_uptake', 'module', ''; ...
    'aa_methionine_net_uptake', 'module', ''; ...
    'aa_phenylalanine_net_uptake', 'module', ''; ...
    'aa_proline_net_uptake', 'module', ''; ...
    'aa_serine_net_uptake', 'module', ''; ...
    'aa_threonine_net_uptake', 'module', ''; ...
    'aa_tryptophan_net_uptake', 'module', ''; ...
    'aa_tyrosine_net_uptake', 'module', ''; ...
    'aa_valine_net_uptake', 'module', ''; ...
    'aa_total_uptake', 'module', ''; ...
    'aa_total_secretion', 'module', ''; ...
    'aa_essential_uptake', 'module', ''; ...
    'aa_uptake_entropy', 'module', ''; ...
    'aa_n_taken_up', 'module', ''; ...
    'ros_h2o2_production', 'alias', 'turnover_h2o2_c'; ...
    'ros_h2o2_removal', 'alias', 'turnover_h2o2_c'; ...
    'ros_h2o2_removal_per_production', 'core', ''; ...
    'ros_h2o2_removed_by_catalase', 'module', ''; ...
    'ros_h2o2_removed_by_glutathione_peroxidase', 'alias', 'rxn_gpx'; ...
    'ros_h2o2_removed_by_peroxiredoxin', 'module', ''; ...
    'ros_superoxide_production', 'module', ''; ...
    'ros_superoxide_dismutated', 'module', ''; ...
    'ros_glutathione_regeneration_cover', 'alias', 'bud_nadph_to_glutathione_reduction'; ...
    'ros_glutathione_cycle_closure', 'core', ''; ...
    'bud_nadph_to_glutathione_reduction', 'core', ''; ...
    'bud_oxppp_cover_of_glutathione_reduction', 'core', ''; ...
    'conv_atp_per_carbon', 'core', ''; ...
    'conv_nadph_per_carbon', 'core', ''; ...
    'conv_lactate_per_glucose_carbon', 'core', ''; ...
    'conv_co2_per_carbon', 'alias', 'elem_c_to_co2'; ...
    'cap_int_mean_saturation', 'module', ''; ...
    'cap_int_at_bound', 'module', ''; ...
    'cap_ext_mean_saturation', 'module', ''; ...
    'cap_ext_at_bound', 'module', ''; ...
    'ex_total_uptake', 'module', ''; ...
    'ex_total_secretion', 'generic', ''; ...
    'ex_n_uptaken', 'generic', ''; ...
    'ex_n_secreted', 'generic', ''; ...
    'ex_glucose_net', 'module', ''; ...
    'ex_lactate_net', 'alias', 'rxn_ldh'; ...
    'ex_oxygen_net', 'module', ''; ...
    'ex_carbon_dioxide_net', 'module', ''; ...
    'ex_hydrogen_peroxide_net', 'module', ''; ...
    'ex_reduced_glutathione_net', 'module', ''; ...
    'ex_pyruvate_net', 'module', ''; ...
    'ex_cysteine_net', 'module', ''; ...
    'ex_phosphate_net', 'module', ''; ...
    'cpl_boundary_share', 'core', ''; ...
    'cpl_internal_per_boundary', 'alias', 'cpl_boundary_share'; ...
    'cpl_secretion_over_uptake', 'alias', 'ext_directionality'; ...
    };
% END DESCRIPTOR STATUS
end


function out = keepFields(values, names)
out = struct();
for k = 1:numel(names)
    out.(names{k}) = values.(names{k});
end
end


function features = computeAll(declared, V)
% Compute in declaration order, so a ratio can read the descriptors it needs.
nStates = size(V, 1);
vals = struct();
for k = 1:numel(declared)
    if declared(k).computable
        value = declared(k).compute(V, vals);
        value = double(value(:))';
    else
        value = nan(1, nStates);
    end
    vals.(declared(k).name) = value;
end
features = vals;
end


% ==========================================================================
% the small numerical helpers, kept explicit rather than imported
% ==========================================================================
function out = safeDivide(numerator, denominator)
% Divide, returning NaN where the denominator vanishes. Zero is a
% measurement; a quotient by nothing is not.
shape = size(denominator);
numerator = numerator(:); denominator = denominator(:);
out = nan(size(denominator));
usable = abs(denominator) > 0;
out(usable) = numerator(usable) ./ denominator(usable);
out = reshape(out, shape);
end


function out = giniOf(A)
ordered = sort(A, 2);
n = size(ordered, 2);
index = 1:n;
total = sum(ordered, 2);
out = safeDivide(2 * (ordered * index') - (n + 1) * total, n * total);
end


function out = entropyOf(A)
total = sum(A, 2);
p = A ./ max(total, realmin);
p(total == 0, :) = 0;
terms = zeros(size(p));
positive = p > 0;
terms(positive) = p(positive) .* log(p(positive));
out = -sum(terms, 2);
end


function out = cvOf(A)
if size(A, 2) < 2
    out = nan(size(A, 1), 1);
    return
end
out = safeDivide(std(A, 0, 2), mean(A, 2));     % sample standard deviation
end


function out = topShare(A, fraction)
k = max(1, ceil(size(A, 2) * fraction));
ordered = sort(A, 2);
top = sum(ordered(:, end-k+1:end), 2);
out = safeDivide(top, sum(A, 2));
end


% ==========================================================================
% elements, amino acids, carriers and signatures
% ==========================================================================
function n = elementCount(formula, symbol)
% Atoms of SYMBOL in a formula such as C6H12O6; 0 if absent or unreadable.
n = 0;
if ~ischar(formula) || isempty(formula), return, end
tokens = regexp(formula, '([A-Z][a-z]?)(\d*(?:\.\d+)?)', 'tokens');
for k = 1:numel(tokens)
    if strcmp(tokens{k}{1}, symbol)
        if isempty(tokens{k}{2})
            n = n + 1;
        else
            n = n + str2double(tokens{k}{2});
        end
    end
end
end


function [acids, essential, tracked] = aminoAcidTable()
% The proteinogenic amino acids: canonical name and identifier stems (BiGG, VMH).
acids = { ...
    'alanine', {'ala__L', 'ala_L'}; 'arginine', {'arg__L', 'arg_L'}; ...
    'asparagine', {'asn__L', 'asn_L'}; 'aspartate', {'asp__L', 'asp_L'}; ...
    'cysteine', {'cys__L', 'cys_L'}; 'glutamine', {'gln__L', 'gln_L'}; ...
    'glutamate', {'glu__L', 'glu_L'}; 'glycine', {'gly'}; ...
    'histidine', {'his__L', 'his_L'}; 'isoleucine', {'ile__L', 'ile_L'}; ...
    'leucine', {'leu__L', 'leu_L'}; 'lysine', {'lys__L', 'lys_L'}; ...
    'methionine', {'met__L', 'met_L'}; 'phenylalanine', {'phe__L', 'phe_L'}; ...
    'proline', {'pro__L', 'pro_L'}; 'serine', {'ser__L', 'ser_L'}; ...
    'threonine', {'thr__L', 'thr_L'}; 'tryptophan', {'trp__L', 'trp_L'}; ...
    'tyrosine', {'tyr__L', 'tyr_L'}; 'valine', {'val__L', 'val_L'}};
essential = {'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', ...
             'phenylalanine', 'threonine', 'tryptophan', 'valine'};
% reported one by one as tracked exchanges, so left out of the per-acid list
tracked = {'cysteine', 'glycine', 'glutamate'};
end


function name = aminoAcidOf(identifier, metName)
% Canonical amino-acid name of a metabolite, by identifier and then by name;
% D forms excluded. '' if it is not an amino acid.
name = '';
acids = aminoAcidTable();
stem = splitMetabolites({identifier});
stem = stem{1};
for k = 1:size(acids, 1)
    if any(strcmp(stem, acids{k, 2}))
        name = acids{k, 1};
        return
    end
end
if ~ischar(metName), return, end
text = strtrim(regexprep(lower(strtrim(metName)), '\(.*?\)', ''));
if isempty(text) || strncmp(text, 'd-', 2), return, end
synonyms = {'aspartic acid', 'aspartate'; 'glutamic acid', 'glutamate'; ...
            'l-aspartate', 'aspartate'; 'l-glutamate', 'glutamate'};
hit = find(strcmp(synonyms(:, 1), text), 1);
if ~isempty(hit), text = synonyms{hit, 2}; end
if strncmp(text, 'l-', 2), text = text(3:end); end
if any(strcmp(acids(:, 1), text)), name = text; end
end


function X = exportOf(V, exCols, exSign, columns)
% Export of the exchanged metabolites (positive out of the cell); COLUMNS, if
% given, picks some of them.
X = V(:, exCols) .* exSign;
if nargin >= 4, X = X(:, columns); end
end


function P = exchangeParts(V, exCols, exSign, weights, uptake)
% Uptake (or secretion) of every exchanged metabolite, weighted by atoms.
X = exportOf(V, exCols, exSign);
if uptake
    P = max(-X, 0) .* weights;
else
    P = max(X, 0) .* weights;
end
end


function share = partShare(P, columns)
% Share of the row total carried by some columns.
share = safeDivide(sum(P(:, columns), 2), sum(P, 2));
end


function s = sumColumns(M, columns)
s = sum(M(:, columns), 2);
end


function value = diversityOf(parts, kind)
% Entropy, effective number (Hill number of order 2) or dominant share of the
% composition in each row; NaN where the row is empty.
total = sum(parts, 2);
f = zeros(size(parts));
usable = total > 0;
f(usable, :) = parts(usable, :) ./ total(usable);
switch kind
    case 'entropy'
        terms = zeros(size(f));
        positive = f > 0;
        terms(positive) = f(positive) .* log(f(positive));
        value = -sum(terms, 2);
    case 'effective'
        value = 1 ./ sum(f .^ 2, 2);
    case 'dominant'
        value = max(f, [], 2);
end
value(~usable) = NaN;
end


function A = acidMatrix(V, exCols, exSign, byAcid)
% Net export of every recognised amino acid, one column per acid.
X = exportOf(V, exCols, exSign);
A = zeros(size(V, 1), numel(byAcid));
for k = 1:numel(byAcid)
    A(:, k) = sum(X(:, byAcid{k}), 2);
end
end


function value = aminoAggregate(A, kind, essential, tol)
uptake = max(-A, 0);
switch kind
    case 'uptake',    value = sum(uptake, 2);
    case 'secretion', value = sum(max(A, 0), 2);
    case 'essential', value = sum(uptake(:, essential), 2);
    case 'entropy',   value = diversityOf(uptake, 'entropy');
    case 'count',     value = sum(uptake > tol, 2);
end
end


function coefficient = carrierCoefficient(Sint, metBase, metComp, stems, extracellular)
% Net coefficient of a carrier (all its intracellular forms) in every internal
% reaction, 0 where the reaction only moves it; [] if the model lacks it.
rows = find(ismember(metBase, stems) & ~ismember(metComp, extracellular));
coefficient = [];
if isempty(rows), return, end
coefficient = full(sum(Sint(rows, :), 1));
every = full(Sint(ismember(metBase, stems), :));
moves = any(every > 0, 1) & any(every < 0, 1);
coefficient(moves) = 0;
end


function hits = reactionSignature(Sint, metBase, consumes, produces, only, partner)
% Internal reactions (positions among the internal columns) that consume and
% produce the given metabolites; ONLY admits nothing else consumed but protons;
% PARTNER requires a consumed metabolite whose identifier starts with it.
hits = [];
for j = 1:size(Sint, 2)
    column = full(Sint(:, j));
    taken = unique(metBase(column < 0));
    given = unique(metBase(column > 0));
    if ~all(ismember(consumes, taken)) || ~all(ismember(produces, given)), continue, end
    if only && ~isempty(setdiff(taken, [consumes, {'h'}])), continue, end
    if ~isempty(partner) && ~any(strncmp(taken, partner, numel(partner))), continue, end
    hits(end+1) = j;
end
end


function M = madeBy(V, columns, coefficient)
M = max(V(:, columns) .* coefficient, 0);
end


function M = usedBy(V, columns, coefficient)
M = max(-V(:, columns) .* coefficient, 0);
end
