function [features, info, comparison] = flux2feature(model, fluxData, flux2Compare, varargin)
% FLUX2FEATURE  Interpretable descriptors from a metabolic flux vector.
%
%   One file, base MATLAB, no toolbox. Give it a model and a flux vector and
%   it returns named biochemical descriptors instead of a list of reaction
%   fluxes. Ask it first and it tells you which descriptors your model can
%   support and why the others cannot be computed.
%
%   INFO = FLUX2FEATURE(MODEL)
%       The audit, and nothing else: one entry per declared descriptor, with
%       its family, block, units, what it requires from a reconstruction,
%       whether this model has it and, when it does not, why. Call this
%       before computing anything. STRUCT2TABLE(INFO) prints it readably.
%
%   FEATURES = FLUX2FEATURE(MODEL, FLUXDATA)
%       Descriptors of one state, or of many. FLUXDATA is a column vector of
%       length numel(MODEL.rxns), or a matrix with the reactions down the
%       rows and one state per column — the shape SAMPLECBMODEL returns.
%       FEATURES is a struct whose fields are the descriptor names, each a
%       1-by-nStates row.
%
%   [FEATURES, INFO] = FLUX2FEATURE(...) also returns the audit.
%
%   [FEATURES, INFO, COMPARISON] = FLUX2FEATURE(MODEL, FLUXDATA, FLUX2COMPARE)
%       FLUX2COMPARE is a reference state — a healthy reference, a wild-type
%       solution, the median of another sample. COMPARISON reports, per
%       descriptor, the reference value, the difference and the log2 ratio,
%       which is usually what you wanted when you computed two sets of
%       features and started subtracting them by hand.
%
%   ... = FLUX2FEATURE(..., 'Name', VALUE) accepts:
%       'Which'      'computable' (default) returns only what this model
%                    supports; 'all' returns every declared descriptor, NaN
%                    where it cannot be computed — use this when comparing
%                    two reconstructions, so the fields line up and an
%                    absence stays visibly an absence; a cellstr returns
%                    exactly those descriptors.
%       'Tolerance'  below this a flux counts as zero (default 1e-6).
%       'Exchange'   logical vector marking the boundary reactions. By
%                    default a reaction is a boundary reaction if its
%                    identifier starts with EX_, and failing that if it
%                    touches exactly one metabolite — which also catches
%                    sinks and demand reactions, so pass this explicitly if
%                    your model distinguishes them.
%       'Pathways', 'PathwayReactions', 'TurnoverMetabolites',
%       'NamedReactions', 'TrackedExchanges', 'CofactorPairs', 'Currency'
%                    the declarations. The defaults describe human
%                    erythrocyte metabolism (RBC-GEM identifiers); pass your
%                    own for another reconstruction. Everything that does
%                    not depend on names keeps working untouched.
%
%   Three rules are enforced, and they are the point of the thing:
%     * every descriptor is DECLARED before it is computed — name, units,
%       block, what it needs from the model, what it means;
%     * a descriptor whose requirement is missing returns NaN, NEVER 0,
%       because zero is a measurement and absence is not;
%     * nothing is renamed, reordered or dropped silently: the audit is the
%       contract and the computation honours it.
%
%   Example:
%       info = flux2feature(model);
%       nnz([info.computable])                       % what this model supports
%       f = flux2feature(model, points);             % points: nRxns x nSamples
%       f.path_pentose_phosphate_activity            % one row, one value per state
%       [f, ~, cmp] = flux2feature(model, points, vHealthy);
%
%   See also TEST_FLUX2FEATURE.

%   Companion to the Python implementation flux2feature.py; the two agree on
%   every shared descriptor to 1e-5 relative on a real sampled state of
%   RBC-GEM. Version 1.0.

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
declared = declare(model, opt);

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
opt.which      = 'computable';
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
        case 'pathways',            opt.pathways = args{k+1};
        case 'pathwayreactions',    opt.pathwayReactions = args{k+1};
        case 'turnovermetabolites', opt.turnoverMetabolites = args{k+1};
        case 'namedreactions',      opt.namedReactions = args{k+1};
        case 'trackedexchanges',    opt.trackedExchanges = args{k+1};
        case 'cofactorpairs',       opt.cofactorPairs = args{k+1};
        case 'efficiencyratios',    opt.efficiencyRatios = args{k+1};
        otherwise, error('flux2feature:option', 'unknown option "%s"', args{k});
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
if iscell(which)
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
    error('flux2feature:which', ...
        'Which must be ''computable'', ''all'' or a cell array of names');
end
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
