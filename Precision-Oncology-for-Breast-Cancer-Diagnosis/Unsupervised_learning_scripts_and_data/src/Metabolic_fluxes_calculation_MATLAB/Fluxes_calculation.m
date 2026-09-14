initCobraToolbox(false);
changeCobraSolver('gurobi','all');


folderPath = 'C:\Users\aleja\Documents\Modelos_actual';
files = dir(fullfile(folderPath, '*.mat'));
disp(length(files))
numModels = length(files);

featureMatrix = [];
modelNames    = {};
maxBiomass    = [];

% =========================================================
% === REACTION LISTS DEFINITION ===
% =========================================================

% --- Warburg ---
warburgRxnsIDs_Glycolysis = {
    'EX_glc__D_e',
    'PFK',
    'LDH_L',
    'EX_lac__L_e'
};
warburgRxnsIDs_OXPHOS = {
    'CSm',        % TCA entry
    'SUCD1m',     % complex II
    'ICDHxm',     % mitochondrial TCA
    'ATPS4m'      % mitochondrial ATP synthase
};
% --- CU / EA ---
cuRxnsIDs = {'EX_glc__D_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e', ...
    'EX_met__L_e','EX_phe__L_e','EX_thr__L_e','EX_trp__L_e', ...
    'EX_val__L_e','EX_ala__L_e','EX_asn__L_e','EX_gln__L_e', ...
    'EX_tyr__L_e','EX_arg__L_e','EX_gly_e','EX_ser__L_e','EX_glu__L_e','EX_pyr_e'};

eaRxnsIDs = {'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cys__L_e', ...
    'EX_gln__L_e','EX_glu__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e', ...
    'EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_phe__L_e','EX_pro__L_e', ...
    'EX_ser__L_e','EX_thr__L_e','EX_trp__L_e','EX_tyr__L_e','EX_val__L_e'};

kcatKM   = [1.4,140,3.2,4.5,0.7,1.9,2.3,0.35,1.1,2.6, ...
            2.6,0.8,3.0,4.2,1.5,0.95,1.2,0.6,1.3,2.6];
W        = ones(1, length(kcatKM)) * 10;
EAcoeffs = W ./ kcatKM;


% --- Redox NAD+/NADH ---
redoxRxns_NADplus = {'LDH_L','LDH_Lm','MDH','MDHm'};
redoxRxns_NADH    = {'PDHm','GAPD','ICDHxm','ICDHym'};
% --- Metabolic flexibility (MFI) ---
altSubstrateRxns = {
    'EX_gln__L_e',
    'EX_lac__L_e',
    'EX_pyr_e',
    'EX_glc__D_e'
};
MFI_threshold = 1e-6;
% --- Anabolism vs catabolism ---
anabolicRxns  = {'PGI','TALA','TKT1','RPI'};
catabolicRxns = {'PDHm','PCm','CSm','PYK'};

% --- Oncometabolites ---
oncometRxns  = {'EX_lac__L_e','EX_succ_e','EX_akg_e'};
oncometNames = {'Lactate','Succinate','AlphaKG'};
% --- NADPH ---
nadphRxns = {'MTHFD','MTHFDm','ICDHxm','ICDHym'};

% --- Lipids ---
lipidRxns_sat   = {'FAOXC140','ACACT7m'};
lipidRxns_unsat = {'C161CPT22','C30CPT1'};
lipidRxns_PL    = {'PSSA1_hs','CEPTC'};
% =========================================================
% === STEP 1: MASTER LIST OF REACTIONS ===
% =========================================================
nSols    = 5;
solNames = {'FBA','pFBA','L1w','L2','L2w'};
% =========================================================
% === STEP 1: MASTER LIST OF REACTIONS AND SUBSYSTEMS ===
% =========================================================

masterRxns = {};
masterSubsystems = {};

for modelIdx = 1:numModels

    modelName = files(modelIdx).name;
    fprintf('Scanning model %d/%d: %s\n', modelIdx, numModels, modelName);

    try
        model = prepareModel(fullfile(folderPath, modelName));
    catch ME
        warning('Error preparing model %s: %s', modelName, ME.message);
        continue;
    end

    % Reactions
    masterRxns = union(masterRxns, model.rxns);

    % Subsystems
    ss = getUniqueSubsystems(model);
    masterSubsystems = union(masterSubsystems, ss);
end


% =========================================================
% === STEP 2: MAIN ANALYSIS LOOP ===
% =========================================================

fprintf('\nModels detected: %d\n\n', numModels);

for modelIdx = 1:numModels

    modelName = files(modelIdx).name;   % FIX: was 'fileName', undefined variable
    fprintf('Processing model %d/%d: %s\n', modelIdx, numModels, modelName);

    try
        model = prepareModel(fullfile(folderPath, modelName));
    catch ME
        warning('Error preparing model %s: %s', modelName, ME.message);
        continue;
    end

    % ==============================
    % BIOMASS
    % ==============================
    biomassRxn = 'biomass_ob';
    bioIdx = find(strcmp(model.rxns, biomassRxn));

    if isempty(bioIdx)
        warning('⚠️  Biomass reaction not found in: %s', modelName);
        maxBiomass(end+1) = NaN;
        modelNames{end+1} = strrep(modelName, '.mat', '');
        continue
    end
    subSystemsToTrack = masterSubsystems;
    % -------------------------------------------------------
    % OPTIMIZATIONS
    % -------------------------------------------------------

    % SOL 1: FBA — maximize biomass
    sol_FBA = optimizeCbModel(changeObjective(model,biomassRxn));

    if sol_FBA.stat == 1
        bioVal = sol_FBA.f;
    else
        warning('⚠️  FBA without optimal solution for: %s (stat=%d)', modelName, sol_FBA.stat);
        bioVal = NaN;
    end
    model1 = model;
    % SOL 2: pFBA — minimize total flux (L1) at biomass optimum
    model_pFBA           = model1;
    model_pFBA.lb(bioIdx) = sol_FBA.f * 0.99;
    sol_pFBA             = optimizeCbModel(model_pFBA, 'min', 'one');

    % SOL 3: L1 weighted by gene expression
    model2 = model;
    model_L1w            = model2;
    model_L1w.lb(bioIdx)  = sol_FBA.f * 0.99;
    model_L1w.expressionRxns(isnan(model_L1w.expressionRxns)) = 0;
    weights_L1           = model_L1w.expressionRxns + abs(min(model_L1w.expressionRxns)) + 1e-6;
    sol_L1w              = optimizeCbModel(model_L1w, 'min', weights_L1);

    % SOL 4: L2 — minimize Euclidean norm of fluxes
    model3 = model;
    model_L2             = model3;
    model_L2.lb(bioIdx)   = sol_FBA.f * 0.99;
    sol_L2               = optimizeCbModel(model_L2, 'min', 1e-6);

    if sol_L2.stat ~= 1
        warning('⚠️  L2 without optimal solution for: %s (stat=%d)', modelName, sol_L2.stat);  % FIX: was fileName
        sol_L2.x = zeros(length(model.rxns), 1);
    end

    % SOL 5: L2 weighted by gene expression
    model4 = model;
    model_L2w            = model4;
    model_L2w.lb(bioIdx)  = sol_FBA.f * 0.99;
    model_L2w.expressionRxns(isnan(model_L2w.expressionRxns)) = 0;
    weights_L2           = model_L2w.expressionRxns + abs(min(model_L2w.expressionRxns)) + 1e-6;
    sol_L2w              = optimizeCbModel(model_L2w, 'min', weights_L2 * 1e-6);

    if sol_L2w.stat ~= 1
        warning('⚠️  L2w without optimal solution for: %s (stat=%d)', modelName, sol_L2w.stat);  % FIX: was fileName
        sol_L2w.x = zeros(length(model.rxns), 1);
    end

    fluxSolutions = {sol_FBA.x, sol_pFBA.x, sol_L1w.x, sol_L2.x, sol_L2w.x};
    solList       = {sol_FBA,   sol_pFBA,   sol_L1w,   sol_L2,   sol_L2w};

    % Biomass flux (for normalization)
    bioFlux = sol_FBA.f;
    if isempty(bioFlux) || bioFlux <= 0, bioFlux = 1; end

    % -------------------------------------------------------
    % FLUX VECTOR (masterRxns x nSols)
    % -------------------------------------------------------
    fluxVector = zeros(1, length(masterRxns) * nSols);
    [~, idxInMaster] = ismember(model.rxns, masterRxns);
    for s = 1:nSols
        fluxVector((s-1)*length(masterRxns) + idxInMaster) = fluxSolutions{s};
    end

    % -------------------------------------------------------
    % METRICS
    % -------------------------------------------------------

    % CU and EA
    CU = zeros(1, nSols);
    EA = zeros(1, nSols);
    for s = 1:nSols
        CU(s) = calcCU(solList{s}, model, cuRxnsIDs);
        EA(s) = calcEA(solList{s}, model, eaRxnsIDs, EAcoeffs);
    end

    % Warburg Index
    WarburgIndex = zeros(1, nSols);
    for s = 1:nSols
        gly = calcRobustSum(fluxSolutions{s}, model, warburgRxnsIDs_Glycolysis);
        oxp = calcRobustSum(fluxSolutions{s}, model, warburgRxnsIDs_OXPHOS);
        WarburgIndex(s) = gly.value - oxp.value;
    end

    % ATP
    ATPConsumption = zeros(1, nSols);
    ATPProduction  = zeros(1, nSols);
    atpMetID = find(strcmp(model.mets, 'atp_c'));
    if ~isempty(atpMetID)
        S_atp = model.S(atpMetID, :);
        for s = 1:nSols
            atp_bal = S_atp' .* fluxSolutions{s};
            ATPProduction(s)  = sum(atp_bal(atp_bal >  1e-9));
            ATPConsumption(s) = abs(sum(atp_bal(atp_bal < -1e-9)));
        end
    else
        ATPConsumption(:) = NaN;
        ATPProduction(:)  = NaN;
    end


SA_Matrix = zeros(length(subSystemsToTrack), nSols);

for s = 1:nSols
    flux = fluxSolutions{s};
    
    for ss_idx = 1:length(subSystemsToTrack)
        
        rxnsBool = strcmp(model.subSystems, subSystemsToTrack{ss_idx});
        numRxns  = sum(rxnsBool);
        
        if numRxns > 0
            SA_Matrix(ss_idx, s) = sum(abs(flux(rxnsBool))) / numRxns;
        end
        
    end
end

SubsystemActivity = reshape(SA_Matrix', 1, []);

    % Redox Index
    RedoxIndex = zeros(1, nSols);
    for s = 1:nSols
        rPlus = calcRobustSum(fluxSolutions{s}, model, redoxRxns_NADplus);
        rNADH = calcRobustSum(fluxSolutions{s}, model, redoxRxns_NADH);
        if (rPlus.coverage + rNADH.coverage) / 2 >= 0.5
            RedoxIndex(s) = rPlus.value - rNADH.value;
        else
            RedoxIndex(s) = NaN;
        end
    end

    % MFI — Metabolic flexibility
    MFI = zeros(1, nSols);
    for s = 1:nSols
        flux   = fluxSolutions{s};
        count  = 0;
        nFound = 0;
        for k = 1:length(altSubstrateRxns)
            idx = find(strcmp(model.rxns, altSubstrateRxns{k}));
            if ~isempty(idx)
                nFound = nFound + 1;
                if flux(idx) < -MFI_threshold
                    count = count + 1;
                end
            end
        end
        MFI(s) = ternary(nFound >= 2, count, NaN);
    end

    % Anabolism Score
    AnabolismScore = zeros(1, nSols);
    for s = 1:nSols
        rAna = calcRobustSum(fluxSolutions{s}, model, anabolicRxns);
        rCat = calcRobustSum(fluxSolutions{s}, model, catabolicRxns);
        if rAna.coverage >= 0.5 && rCat.coverage >= 0.5
            AnabolismScore(s) = rAna.value / (rCat.value + eps);
        else
            AnabolismScore(s) = NaN;
        end
    end

    % Oncometabolite Score
    OncometScore = zeros(length(oncometRxns), nSols);
    for s = 1:nSols
        flux = fluxSolutions{s};
        for k = 1:length(oncometRxns)
            idx = find(strcmp(model.rxns, oncometRxns{k}));
            if ~isempty(idx)
                OncometScore(k, s) = max(0, flux(idx));
            else
                OncometScore(k, s) = NaN;
            end
        end
    end
    OncometFlat = reshape(OncometScore', 1, []);

    % NADPH Demand
    NADPHdemand = zeros(1, nSols);
    for s = 1:nSols
        r = calcRobustSum(fluxSolutions{s}, model, nadphRxns);
        if r.coverage >= 0.5
            NADPHdemand(s) = r.value / bioFlux;
        else
            NADPHdemand(s) = NaN;
        end
    end

    % TCA Completeness
    TCA_completeness = zeros(1, nSols);
    tcaBool = strcmp(model.subSystems, 'Citric acid cycle');
    nTCA    = sum(tcaBool);
    for s = 1:nSols
        if nTCA > 0
            activeRxns = sum(abs(fluxSolutions{s}(tcaBool)) > 1e-6);
            TCA_completeness(s) = activeRxns / nTCA;
        else
            TCA_completeness(s) = NaN;
        end
    end

    % Lipid Profile
    LipidSat   = zeros(1, nSols);
    LipidUnsat = zeros(1, nSols);
    LipidPL    = zeros(1, nSols);
    for s = 1:nSols
        rSat   = calcRobustSum(fluxSolutions{s}, model, lipidRxns_sat);
        rUnsat = calcRobustSum(fluxSolutions{s}, model, lipidRxns_unsat);
        rPL    = calcRobustSum(fluxSolutions{s}, model, lipidRxns_PL);
        LipidSat(s)   = ternary(rSat.coverage   >= 0.5, rSat.value,   NaN);
        LipidUnsat(s) = ternary(rUnsat.coverage  >= 0.5, rUnsat.value, NaN);
        LipidPL(s)    = ternary(rPL.coverage     >= 0.5, rPL.value,    NaN);
    end

    % Glutamine Dependence
    GlnDependence = zeros(1, nSols);
    glnIdx = find(strcmp(model.rxns, 'EX_gln__L_e'));
    glcIdx = find(strcmp(model.rxns, 'EX_glc__D_e'));
    for s = 1:nSols
        flux   = fluxSolutions{s};
        hasGln = ~isempty(glnIdx);
        hasGlc = ~isempty(glcIdx);
        if hasGln && hasGlc
            fGln = abs(min(0, flux(glnIdx)));
            fGlc = abs(min(0, flux(glcIdx)));
            GlnDependence(s) = fGln / (fGlc + fGln + eps);
        else
            GlnDependence(s) = NaN;
        end
    end

    % -------------------------------------------------------
    % FINAL FEATURE VECTOR
    % -------------------------------------------------------

    features = [ ...
        fluxVector, ...
        CU, EA, ...
        WarburgIndex, ...
        ATPConsumption, ATPProduction, ...
        SubsystemActivity, ...
        RedoxIndex, ...
        MFI, ...
        AnabolismScore, ...
        OncometFlat, ...
        NADPHdemand, ...
        TCA_completeness, ...
        LipidSat, LipidUnsat, LipidPL, ...
        GlnDependence ...
    ];

    featureMatrix = [featureMatrix; features];
    modelNames{end+1, 1} = strrep(modelName, '.mat', '');   % FIX: was fileName
    maxBiomass(end+1)    = bioVal;

    fprintf('✅ [%d/%d] %s — Max biomass = %.6f\n', modelIdx, numModels, modelName, bioVal);  % FIX: was i, length(files), fileName

    if mod(modelIdx, 50) == 0   % FIX: was mod(i, 50)
        fprintf('Processed %d / %d models...\n', modelIdx, numModels);
    end
end

% =========================================================
% === STEP 3: COLUMN NAMES ===
% =========================================================
colNames = {};

% Fluxes
for s = 1:nSols
    for r = 1:length(masterRxns)
        colNames{end+1} = sprintf('Flux_%s_%s', masterRxns{r}, solNames{s});
    end
end

% Scalar metrics
for s = 1:nSols, colNames{end+1} = sprintf('CU_%s',                  solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('EA_%s',                  solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('WarburgIndex_%s',        solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('ATPConsumption_%s',      solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('ATPProduction_%s',       solNames{s}); end

% Subsystem Activity
for ss_idx = 1:length(subSystemsToTrack)
    ssName = strrep(subSystemsToTrack{ss_idx}, ' ', '_');
    for s = 1:nSols
        colNames{end+1} = sprintf('SA_%s_%s', ssName, solNames{s});
    end
end

% New scalar metrics
for s = 1:nSols, colNames{end+1} = sprintf('RedoxIndex_%s',     solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('MFI_%s',            solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('AnabolismScore_%s', solNames{s}); end

% Oncometabolites
for k = 1:length(oncometNames)
    for s = 1:nSols
        colNames{end+1} = sprintf('Oncomet_%s_%s', oncometNames{k}, solNames{s});
    end
end

for s = 1:nSols, colNames{end+1} = sprintf('NADPHdemand_%s',      solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('TCA_completeness_%s', solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('LipidSat_%s',         solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('LipidUnsat_%s',       solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('LipidPL_%s',          solNames{s}); end
for s = 1:nSols, colNames{end+1} = sprintf('GlnDependence_%s',    solNames{s}); end

% =========================================================
% === DIMENSION VERIFICATION ===
% =========================================================
nColsMatrix = size(featureMatrix, 2);
nColNames   = length(colNames);
if nColsMatrix ~= nColNames
    error('MISMATCH: featureMatrix has %d columns but colNames has %d entries.\n', ...
          nColsMatrix, nColNames);
else
    fprintf('OK: %d columns in featureMatrix == %d names in colNames\n', ...
            nColsMatrix, nColNames);
end

% =========================================================
% === EXPORT ===
% =========================================================
T_biomass = table(modelNames(:), maxBiomass(:), ...
                  'VariableNames', {'Model', 'MaxBiomass_FBA'});
writetable(T_biomass, fullfile(folderPath, 'MaxBiomass_perModel.csv'));
fprintf('\n✅ Maximum biomass saved.\n');

featureTable       = array2table(featureMatrix, 'VariableNames', matlab.lang.makeValidName(colNames));
featureTable.Model = modelNames;
outFile = fullfile(folderPath, 'FeatureMatrix_TumorPhenotype_norm2agregado.csv');
writetable(featureTable, outFile);
fprintf('\n✅ Complete feature matrix saved to:\n%s\n', outFile);
fprintf('Dimensions: %d models x %d features\n', size(featureMatrix, 1), size(featureMatrix, 2));

% =========================================================
% === HELPER FUNCTIONS ===
% =========================================================


function result = calcRobustSum(flux, model, rxnIDs)
    total  = 0;
    nFound = 0;
    for j = 1:length(rxnIDs)
        idx = find(strcmp(model.rxns, rxnIDs{j}));
        if ~isempty(idx)
            total  = total + abs(flux(idx));
            nFound = nFound + 1;
        end
    end
    result.value    = total;
    result.coverage = nFound / length(rxnIDs);
end

function CU = calcCU(sol, model, cuRxnsIDs)
    CU = 0;
    for j = 1:length(cuRxnsIDs)
        rxnIdx = find(strcmp(model.rxns, cuRxnsIDs{j}));
        if isempty(rxnIdx), continue; end
        flux = sol.x(rxnIdx);
        if flux < 0
            mets = model.mets(model.S(:, rxnIdx) < 0);
            for m = 1:length(mets)
                metIdx = find(strcmp(model.mets, mets{m}));
                if ~isempty(metIdx)
                    coef = -model.S(metIdx, rxnIdx);
                    nC   = countCarbons(model.metFormulas{metIdx});
                    CU   = CU + nC * coef * (-flux);
                end
            end
        end
    end
end

function EA = calcEA(sol, model, eaRxnsIDs, EAcoeffs)
    EA = 0;
    for k = 1:length(eaRxnsIDs)
        rxnIdx = find(strcmp(model.rxns, eaRxnsIDs{k}));
        if isempty(rxnIdx), continue; end
        EA = EA + EAcoeffs(k) * sol.x(rxnIdx);
    end
end

function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end

function subSystemsToTrack = getUniqueSubsystems(model)

    % --- Normalize subSystems ---
    subList = model.subSystems;

    for i = 1:length(subList)

        if iscell(subList{i})
            if isempty(subList{i})
                subList{i} = '';
            else
                subList{i} = subList{i}{1};
            end
        end

        if isempty(subList{i})
            subList{i} = '';
        end
    end

    % --- Remove empty entries ---
    subList = subList(~cellfun(@isempty, subList));

    % --- Get unique entries ---
    subSystemsToTrack = unique(subList);

end