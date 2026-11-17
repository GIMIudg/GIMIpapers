%% Interpretable descriptors of flux states with flux2feature
%% Author(s): Germán Preciat, Ramyet Sotelo Rodriguez, Carolina Elizabeth Gómez Márquez, Universidad de Guadalajara
%% Reviewer(s):
%% INTRODUCTION
% A flux vector of a constraint-based model, whether a flux balance optimum or
% a sampled state, is a list of reaction rates. |flux2feature| maps it to
% $\Phi_M(v)$, a vector of declared, interpretable descriptors: pathway
% activity and share, metabolite turnover, cofactor allocation, NADPH and
% glutathione budgets, conversion yields, boundary exchange and compartment
% shares. Three rules hold throughout: every descriptor is declared before it
% is computed (name, units, block, requirement, reading); a descriptor whose
% requirement the model lacks is |NaN|, never 0; nothing is renamed,
% reordered or dropped silently.
%
% Each descriptor carries a tier: |core| (the default output), |module|,
% |generic| and |alias|. The core is the set of interpretable coordinates
% released with the study; it is not a claim that these are the most
% important descriptors of every system.
%
% The defaults name human erythrocyte metabolism (RBC-GEM). This tutorial
% applies the tool to the _E. coli_ core model: it maps a few identifiers,
% audits what the model can support, describes a flux balance optimum, then
% reads the distribution of the descriptors over sampled feasible states, and
% compares a perturbed state with the reference.
%% EQUIPMENT SETUP
% Initialise the COBRA Toolbox and choose an LP solver.

initCobraToolbox(false)
changeCobraSolver('glpk', 'LP');
%% PROCEDURE
%% 1. Load the model and map its identifiers
% |flux2feature| reads metabolite identifiers in BiGG style (|atp_c|,
% |glc__D_e|). The distributed _E. coli_ core model writes |atp[c]| and
% |glc-D[e]|, so convert them first.

model = getDistributedModel('ecoli_core_model.mat');
model.mets = strrep(regexprep(model.mets, '\[(\w+)\]$', '_$1'), '-', '__');
%%
% The labelled declarations are merged with the defaults: a label given
% replaces its default identifiers and every other label is kept, so the
% coordinates stay the same whatever the model. Here the glucose-6-phosphate
% dehydrogenase is |G6PDH2r|, lactate dehydrogenase |LDH_D|, the pathways
% carry _E. coli_ core's subsystem names, and secreted lactate is D-lactate.

declaration = { ...
    'NamedReactions',   {'g6pd', 'G6PDH2r'; 'ldh', 'LDH_D'}, ...
    'Pathways',         {'glycolysis', {'Glycolysis/Gluconeogenesis'}; ...
                         'pentose_phosphate', {'Pentose Phosphate Pathway'}; ...
                         'transport', {'Transport, Extracellular'}}, ...
    'TrackedExchanges', {'lac__D_e', 'lactate'}};
%% 2. Audit first
% Called without fluxes, |flux2feature| returns the audit: one entry per
% declared descriptor, with its tier and, when the model cannot support it,
% the reason.

info = flux2feature(model, [], [], declaration{:});
fprintf('%d descriptors declared, %d computable on E. coli core\n', ...
        numel(info), nnz([info.computable]));
tiers = {'core', 'module', 'generic', 'alias'};
for k = 1:numel(tiers)
    inTier = strcmp({info.tier}, tiers{k});
    fprintf('  %-8s %3d of %3d computable\n', tiers{k}, ...
            nnz([info(inTier).computable]), nnz(inTier));
end
%%
% The core descriptors this model cannot support, and why:

missing = info(strcmp({info.tier}, 'core') & ~[info.computable]);
disp(table({missing.name}', {missing.reason}', 'VariableNames', {'descriptor', 'reason'}))
%%
% _E. coli_ core has no peroxide and no mitochondrion, so the oxidative
% defence and compartment descriptors are |NaN|: absent, not zero.
%% 3. Describe a flux balance optimum
% The parsimonious optimum: the smallest total flux at maximal growth.

solution = optimizeCbModel(model, 'max', 'one');
phi = flux2feature(model, solution.x, declaration{:});      % the core of Phi
names = {'path_glycolysis_share'; 'path_pentose_phosphate_share'; ...
         'turnover_atp_c'; 'turnover_nadph_c'; 'cof_nadph_nadp_share'; ...
         'conv_atp_per_carbon'; 'cpl_boundary_share'};
values = cellfun(@(n) phi.(n), names);
disp(table(names, values, 'VariableNames', {'descriptor', 'value'}))
%%
% Other tiers on request: |'Which', {'core', 'module'}|, |'Which', 'all'|, or
% a cell array of descriptor names.

every = flux2feature(model, solution.x, 'Which', 'all', declaration{:});
fprintf('%d descriptors, %d finite\n', numel(fieldnames(every)), ...
        nnz(isfinite(cell2mat(struct2cell(every)))));
%% 4. Read the distribution over feasible states
% An optimum is one state of many. Sample the feasible space at 90% of
% maximal growth and map every sampled state; each descriptor becomes a
% distribution. A narrow distribution marks a property the constraints fix;
% a wide one marks a property they leave free (metabolic reserve).

growth = model;
growth.lb(growth.c ~= 0) = 0.9 * solution.f;
options.nStepsPerPoint = 50;
options.nPointsReturned = 500;
[sampled, samples] = sampleCbModel(growth, [], 'CHRR', options);
Phi = flux2feature(sampled, samples, declaration{:});        % one value per state
%%
% The 5th and 95th percentiles are an empirical range of the sampled states:
% not a confidence interval and not a physiological limit.

names = {'path_glycolysis_share'; 'path_pentose_phosphate_share'; ...
         'turnover_nadph_c'; 'cof_nadph_nadp_share'; 'conv_atp_per_carbon'};
percentile = @(x, q) interp1(linspace(0, 100, nnz(isfinite(x))), sort(x(isfinite(x))), q);
low = cellfun(@(n) percentile(Phi.(n), 5), names);
high = cellfun(@(n) percentile(Phi.(n), 95), names);
disp(table(names, low, high, 'VariableNames', {'descriptor', 'p05', 'p95'}))
%% 5. Compare a perturbed state with the reference
% Limit G6PDH2r to 10% of its flux at the optimum, solve again, and compare
% the two states descriptor by descriptor. The third input is the reference;
% the third output gives the difference and the log2 ratio.

perturbed = model;
g6pd = strcmp(model.rxns, 'G6PDH2r');
perturbed.ub(g6pd) = 0.1 * solution.x(g6pd);
perturbedSolution = optimizeCbModel(perturbed, 'max', 'one');
[~, ~, comparison] = flux2feature(model, perturbedSolution.x, solution.x, declaration{:});
changed = comparison(abs([comparison.difference]) > 1e-6);
disp(table({changed.name}', [changed.reference]', [changed.value]', [changed.log2ratio]', ...
     'VariableNames', {'descriptor', 'reference', 'perturbed', 'log2ratio'}))
%% REFERENCES
% Preciat G, Sotelo Rodriguez R, Gómez Márquez CE. From metabolic fluxes to
% functional features: an interpretable representation of constraint-based
% model states (the article accompanying flux2feature).
%
% Haiman ZB, Key A, D'Alessandro A, Palsson BO. RBC-GEM: a genome-scale
% metabolic model for systems biology of the human red blood cell. _PLoS
% Comput Biol_ 21(3): e1012109 (2025).
