%% =========================================================================
%% MULTI-OBJECTIVE METABOLIC MODELING: PARETO SURFACE & PHENOTYPE PROFILING
%% =========================================================================
% Purpose:
%   This pipeline evaluates trade-offs between Carbon Uptake (CU) and Enzyme 
%   Allocation (EA) across genome-scale metabolic models (GEMs) under varying 
%   demands for growth (biomass) and cellular maintenance (ATP hydrolysis).
%
% Methodology:
%   1. Latin Hypercube Sampling (LHS) samples feasible (Biomass, ATP) constraints.
%   2. Epsilon-constraint / stratified sampling explores the bi-objective front (min CU, min EA).
%   3. Multi-dimensional phenotypic metrics (Warburg index, redox balance, 
%      oncometabolite secretion, subsystem activities, etc.) are computed for each state.
%   4. Vectorized Pareto filtering isolates non-dominated metabolic operating states.
%   5. Results are exported to CSV and MAT files for downstream analysis.
% =========================================================================

%% =========================================================================
%% 1. COBRA TOOLBOX & SOLVER INITIALIZATION
%% =========================================================================
initCobraToolbox(false);
changeCobraSolver('gurobi','all');

% Global numerical precision tolerance and shared subsystem names
global EPSILON;
EPSILON = 1e-9;
global SA_ALL_NAMES;

%% =========================================================================
%% 2. FILE PATHS, RUNTIME PARAMETERS & TARGET REACTIONS
%% =========================================================================
modelsFolder = '/Users/eduardoruiz/Documents/MCBCI/MCBCI2/Sistemas metabólicos/Proyecto_Tesis/Modelos_actual';
outputFileCSV   = 'ParetoSurfacenew.csv';
outputFileMAT   = 'ParetoSolutionsnew.mat';
if exist(outputFileCSV,'file'), delete(outputFileCSV); end
if exist(outputFileMAT,'file'), delete(outputFileMAT); end

% Sampling parameters:
% - numSamples: Number of LHS samples for the (Biomass, ATP) space.
% - maxStrata:  Maximum number of grid partitions for Pareto frontier sampling.
numSamples = 100;
maxStrata  = 10;

% Key objective reactions in human GEMs
biomassRxn      = 'biomass_ob';
tpHydrolysisRxn = 'ATP_HYDROLYSIS';

%% =========================================================================
%% 3. RESOURCE CONSTRAINTS: CARBON UPTAKE (CU) & ENZYME ALLOCATION (EA)
%% =========================================================================
% Carbon Uptake (CU): Exchange reactions of carbon-bearing nutrients (glucose, amino acids).
cuRxnsIDs = {'EX_glc__D_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e', ...
    'EX_met__L_e','EX_phe__L_e','EX_thr__L_e','EX_trp__L_e', ...
    'EX_val__L_e','EX_ala__L_e','EX_asn__L_e','EX_gln__L_e', ...
    'EX_tyr__L_e','EX_arg__L_e','EX_gly_e','EX_ser__L_e','EX_glu__L_e','EX_pyr_e'};

% Enzyme Allocation (EA): Exchange reactions for amino acids reflecting proteomic investment.
eaRxnsIDs = {'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cys__L_e', ...
    'EX_gln__L_e','EX_glu__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e', ...
    'EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_phe__L_e','EX_pro__L_e', ...
    'EX_ser__L_e','EX_thr__L_e','EX_trp__L_e','EX_tyr__L_e','EX_val__L_e'};

% Enzyme molecular weights (Daltons) and catalytic efficiencies (kcat/KM):
% EA coefficient = Enzyme_Mass / (kcat / KM), representing protein burden
% per unit flux. kcatM obtained by 
enzymeMasses_Da = [39460,36275,64370,46530,44508,76715,165283,33566,74141,55010, ...
                   134466,68048,26132,51862,70911,45675,34625,47872,50399,140476];
kcatKM          = [1.4,140,3.2,4.5,0.7,1.9,2.3,0.35,1.1,2.6, ...
                   2.6,0.8,3.0,4.2,1.5,0.95,1.2,0.6,1.3,2.6];
EAcoeffs        = enzymeMasses_Da ./ kcatKM;

%% =========================================================================
%% 4. PHENOTYPIC & METABOLIC METRICS DEFINITION
%% =========================================================================
% --- Warburg Index ---
% Compares glycolytic flux against oxidative phosphorylation (OXPHOS) flux.
warburgRxnsIDs_Glycolysis = {'EX_glc__D_e','PFK','LDH_L','EX_lac__L_e'};
warburgRxnsIDs_OXPHOS     = {'NADH2_u10mi','NADHtru','SUCDim','CYOR_u10mi', ...
                              'CYOOm3i','CYOOm2i','ATPS4mi'};

% --- Nitrogen Metabolism & Growth Precursors ---
nitrogenRxnsIDs           = {'ASPTAm','ASNS1','EX_glu__L_e','EX_gln__L_e'};
growthRxnsIDs             = {'MTHFD2m','RNDR1'};

% --- Redox Balance: NAD+ vs NADH Regeneration ---
redoxRxns_NADplus = {'LDH_L','LDH_Lm','MDH','MDHm'};
redoxRxns_NADH    = {'PDHm','GAPD','ICDHxm','ICDHym'};

% --- Metabolic Flexibility Index (MFI) ---
% Evaluates uptake capacity of non-glucose alternative carbon sources.
altSubstrateRxns  = {'EX_gln__L_e','EX_lac__L_e','EX_pyr_e','EX_glc__D_e'};
MFI_threshold     = 1e-6;

% --- Anabolism vs Catabolism ---
anabolicRxns      = {'PGI','TALA','TKT1','RPI'};
catabolicRxns     = {'PDHm','PCm','CSm','PYK'};

% --- Oncometabolite Secretion ---
oncometRxns       = {'EX_lac__L_e','EX_succ_e','EX_akg_e'};
oncometNames      = {'Lactate','Succinate','AlphaKG'};

% --- NADPH Balance (Reductive Biosynthesis & Antioxidant Defense) ---
nadphRxns         = {'MTHFD','MTHFDm','ICDHxm','ICDHym'};

% --- Lipid Metabolism Profile ---
lipidRxns_sat     = {'FAOXC140','ACACT7m'}; % Saturated fatty acids
lipidRxns_unsat   = {'C161CPT22','C30CPT1'}; % Unsaturated fatty acids
lipidRxns_PL      = {'PSSA1_hs','CEPTC'};   % Phospholipids

%% =========================================================================
%% 5. CSV EXPORT HEADERS SETUP
%% =========================================================================
HEADER_SET = false;

base_header        = {'ModelName','CU_real','EA_real','Biomass','ATP'};
pareto_metrics_hdr = {'WarburgIndex','NitrogenAnaplerosis','GrowthMetabolism', ...
                      'ATPConsumption','ATPProduction','RatioCU_EA', ...
                      'TotalTurnOver','MeanSaturation','GrowthEfficiency', ...
                      'MeanFlowCentrality','TotalReversibleFlux'};
new_metrics_hdr    = {'RedoxIndex','MFI','AnabolismScore', ...
                      'Oncomet_Lactate','Oncomet_Succinate','Oncomet_AlphaKG', ...
                      'NADPHdemand','TCA_completeness', ...
                      'LipidSat','LipidUnsat','LipidPL','GlnDependence'};

ParetoSolutions = struct();
modelFiles      = dir(fullfile(modelsFolder,'*.mat'));
numModels       = length(modelFiles);

%% =========================================================================
%% 6. MAIN MODEL ITERATION & EXPLORATION LOOP
%% =========================================================================
for modelIdx = 1:numModels
    modelName = modelFiles(modelIdx).name;
    fprintf('Processing model %d/%d: %s\n', modelIdx, numModels, modelName);

    try
        model = prepareModel(fullfile(modelsFolder, modelName));
    catch ME
        warning('Error preparing model %s: %s', modelName, ME.message);
        continue;
    end

    % ── Initialize CSV Header using first processed model ──────────────────
    if ~HEADER_SET
        [~, SA_ALL_NAMES_temp] = calcSubsystemActivity( ...
            struct('x', zeros(length(model.rxns),1)), model);
        SA_ALL_NAMES = cellfun(@(x) ['SA_' strrep(x,' ','_')], ...
            SA_ALL_NAMES_temp, 'UniformOutput', false);
        SA_ALL_NAMES = SA_ALL_NAMES(:).';

        header = [base_header(:).', pareto_metrics_hdr(:).', ...
                  new_metrics_hdr(:).', SA_ALL_NAMES];
        fid = fopen(outputFileCSV,'w');
        fprintf(fid,'%s\n', strjoin(header,','));
        fclose(fid);
        HEADER_SET = true;
    end

    % ── Determine Maximum Attainable Capacities (ATP & Biomass) ───────────
    solATP     = optimizeCbModel(changeObjective(model, tpHydrolysisRxn));
    solBiomass = optimizeCbModel(changeObjective(model, biomassRxn));
    if solATP.stat ~= 1 || solBiomass.stat ~= 1
        warning('Failed to optimize ATP or Biomass in model: %s', modelName);
        continue;
    end
    FBiomassMax = solBiomass.f;
    FATPMax     = solATP.f;
    bioFlux     = FBiomassMax;
    if bioFlux <= 0, bioFlux = 1; end

    % ── Latin Hypercube Sampling (LHS) across Biomass & ATP requirements ──
    epsilonCombinations      = lhsdesign(numSamples, 2);
    epsilonCombinations(:,1) = epsilonCombinations(:,1) * FBiomassMax;
    epsilonCombinations(:,2) = epsilonCombinations(:,2) * FATPMax;

    % ── Initialize Result Accumulators ────────────────────────────────────
    CU_real_all = []; EA_real_all = []; Biomass_all = []; ATP_all = [];
    % Pareto topological & physiological metrics
    WBI_all=[]; NA_all=[]; GM_all=[]; ATPC_all=[]; ATPP_all=[];
    RatioCU_EA_all=[]; TurnOver_all=[]; Saturation_all=[];
    Efficiency_all=[]; MeanFlowCentrality_all=[]; TotalReversibleFlux_all=[];
    SA_all = zeros(0, length(SA_ALL_NAMES));
    % Phenotype profiling metrics
    Redox_all=[]; MFI_all=[]; Anab_all=[];
    OncLac_all=[]; OncSucc_all=[]; OncAKG_all=[];
    NADPH_all=[]; TCA_all=[];
    LipSat_all=[]; LipUnsat_all=[]; LipPL_all=[]; Gln_all=[];

    % ── Iterate Over Sampled (Biomass, ATP) Demands ───────────────────────
    for i = 1:numSamples
        % Impose sampled lower bounds for growth and maintenance
        model_sample = changeRxnBounds(model, biomassRxn,      epsilonCombinations(i,1), 'l');
        model_sample = changeRxnBounds(model_sample, tpHydrolysisRxn, epsilonCombinations(i,2), 'l');

        solTest = optimizeCbModel(model_sample);
        if solTest.stat ~= 1, continue; end

        % Add pseudo-reactions aggregating CU and EA
        [model_sample, CU_rxnName, CU_demandRxn] = createPseudoCU(model_sample, cuRxnsIDs);
        [model_sample, EA_rxnName, EA_demandRxn] = createPseudoEA(model_sample, eaRxnsIDs, EAcoeffs);

        % Compute theoretical bounds for the Pareto trade-off
        CU_max = optimizeCbModel(changeObjective(model_sample, CU_rxnName), 'max');
        CU_min = optimizeCbModel(changeObjective(model_sample, CU_rxnName), 'min');
        EA_max = optimizeCbModel(changeObjective(model_sample, EA_rxnName), 'max');
        EA_min = optimizeCbModel(changeObjective(model_sample, EA_rxnName), 'min');

        if any([CU_max.stat, CU_min.stat, EA_max.stat, EA_min.stat] ~= 1), continue; end

        m1 = min(max(round(CU_max.f - CU_min.f), 1), maxStrata);
        m2 = min(max(round(EA_max.f - EA_min.f), 1), maxStrata);

        % Stratified random sampling within bounds
        CU_edges  = linspace(CU_min.f, CU_max.f, m1+1);
        EA_edges  = linspace(EA_min.f, EA_max.f, m2+1);
        CU_samples = arrayfun(@(a,b) a + rand()*(b-a), CU_edges(1:end-1), CU_edges(2:end));
        EA_samples = arrayfun(@(a,b) a + rand()*(b-a), EA_edges(1:end-1), EA_edges(2:end));

        % Epsilon-constraint 1: Minimize EA while fixing CU upper bounds
        for j = 1:m1
            model_lp = changeRxnBounds(model_sample, CU_demandRxn, CU_samples(j), 'u');
            model_lp = changeObjective(model_lp, EA_rxnName);
            sol = optimizeCbModel(model_lp, 'min');
            if sol.stat ~= 1, continue; end

            [row, SA_mapped] = buildResultRow(sol, model_lp, ...
                cuRxnsIDs, eaRxnsIDs, EAcoeffs, biomassRxn, tpHydrolysisRxn, ...
                warburgRxnsIDs_Glycolysis, warburgRxnsIDs_OXPHOS, ...
                nitrogenRxnsIDs, growthRxnsIDs, SA_ALL_NAMES, ...
                redoxRxns_NADplus, redoxRxns_NADH, altSubstrateRxns, MFI_threshold, ...
                anabolicRxns, catabolicRxns, oncometRxns, nadphRxns, ...
                lipidRxns_sat, lipidRxns_unsat, lipidRxns_PL, bioFlux);

            [CU_real_all,EA_real_all,Biomass_all,ATP_all, ...
             WBI_all,NA_all,GM_all,ATPC_all,ATPP_all, ...
             RatioCU_EA_all,TurnOver_all,Saturation_all, ...
             Efficiency_all,MeanFlowCentrality_all,TotalReversibleFlux_all, ...
             Redox_all,MFI_all,Anab_all, ...
             OncLac_all,OncSucc_all,OncAKG_all, ...
             NADPH_all,TCA_all,LipSat_all,LipUnsat_all,LipPL_all,Gln_all, ...
             SA_all] = appendRow(row, SA_mapped, ...
             CU_real_all,EA_real_all,Biomass_all,ATP_all, ...
             WBI_all,NA_all,GM_all,ATPC_all,ATPP_all, ...
             RatioCU_EA_all,TurnOver_all,Saturation_all, ...
             Efficiency_all,MeanFlowCentrality_all,TotalReversibleFlux_all, ...
             Redox_all,MFI_all,Anab_all, ...
             OncLac_all,OncSucc_all,OncAKG_all, ...
             NADPH_all,TCA_all,LipSat_all,LipUnsat_all,LipPL_all,Gln_all, ...
             SA_all);
        end

        % Epsilon-constraint 2: Minimize CU while fixing EA upper bounds
        for j = 1:m2
            model_lp = changeRxnBounds(model_sample, EA_demandRxn, EA_samples(j), 'u');
            model_lp = changeObjective(model_lp, CU_rxnName);
            sol = optimizeCbModel(model_lp, 'min');
            if sol.stat ~= 1, continue; end

            [row, SA_mapped] = buildResultRow(sol, model_lp, ...
                cuRxnsIDs, eaRxnsIDs, EAcoeffs, biomassRxn, tpHydrolysisRxn, ...
                warburgRxnsIDs_Glycolysis, warburgRxnsIDs_OXPHOS, ...
                nitrogenRxnsIDs, growthRxnsIDs, SA_ALL_NAMES, ...
                redoxRxns_NADplus, redoxRxns_NADH, altSubstrateRxns, MFI_threshold, ...
                anabolicRxns, catabolicRxns, oncometRxns, nadphRxns, ...
                lipidRxns_sat, lipidRxns_unsat, lipidRxns_PL, bioFlux);

            [CU_real_all,EA_real_all,Biomass_all,ATP_all, ...
             WBI_all,NA_all,GM_all,ATPC_all,ATPP_all, ...
             RatioCU_EA_all,TurnOver_all,Saturation_all, ...
             Efficiency_all,MeanFlowCentrality_all,TotalReversibleFlux_all, ...
             Redox_all,MFI_all,Anab_all, ...
             OncLac_all,OncSucc_all,OncAKG_all, ...
             NADPH_all,TCA_all,LipSat_all,LipUnsat_all,LipPL_all,Gln_all, ...
             SA_all] = appendRow(row, SA_mapped, ...
             CU_real_all,EA_real_all,Biomass_all,ATP_all, ...
             WBI_all,NA_all,GM_all,ATPC_all,ATPP_all, ...
             RatioCU_EA_all,TurnOver_all,Saturation_all, ...
             Efficiency_all,MeanFlowCentrality_all,TotalReversibleFlux_all, ...
             Redox_all,MFI_all,Anab_all, ...
             OncLac_all,OncSucc_all,OncAKG_all, ...
             NADPH_all,TCA_all,LipSat_all,LipUnsat_all,LipPL_all,Gln_all, ...
             SA_all);
        end
    end

    if isempty(CU_real_all), continue; end

    %% ── Vectorized Pareto Filtering (Non-dominated sorting) ───────────────
    X        = [CU_real_all.', EA_real_all.'];
    isPareto = paretoFilter(X);

    %% ── Apply Pareto Filter across all metrics ────────────────────────────
    pX    = X(isPareto,:);
    pBio  = Biomass_all(isPareto);   pATP  = ATP_all(isPareto);
    pWBI  = WBI_all(isPareto);       pNA   = NA_all(isPareto);
    pGM   = GM_all(isPareto);        pATPC = ATPC_all(isPareto);
    pATTP = ATPP_all(isPareto);      pR    = RatioCU_EA_all(isPareto);
    pTO   = TurnOver_all(isPareto);  pSat  = Saturation_all(isPareto);
    pEff  = Efficiency_all(isPareto);pMFC  = MeanFlowCentrality_all(isPareto);
    pTRF  = TotalReversibleFlux_all(isPareto);
    pRedox= Redox_all(isPareto);     pMFI  = MFI_all(isPareto);
    pAnab = Anab_all(isPareto);
    pLac  = OncLac_all(isPareto);    pSucc = OncSucc_all(isPareto);
    pAKG  = OncAKG_all(isPareto);
    pNADPH= NADPH_all(isPareto);     pTCA  = TCA_all(isPareto);
    pLSat = LipSat_all(isPareto);    pLUns = LipUnsat_all(isPareto);
    pLPL  = LipPL_all(isPareto);     pGln  = Gln_all(isPareto);
    pSA   = SA_all(isPareto,:);

    %% ── Export Pareto Front to CSV ────────────────────────────────────────
    fid = fopen(outputFileCSV,'a');
    for k = 1:size(pX,1)
        % Base identifiers & objectives
        fprintf(fid,'%s,%.6f,%.6f,%.6f,%.6f,', ...
            modelName, pX(k,1), pX(k,2), pBio(k), pATP(k));
        % Pareto topological/flux metrics
        rCU_EA = ternary(pX(k,2) > EPSILON, pX(k,1)/pX(k,2), NaN);
        fprintf(fid,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,', ...
            pWBI(k), pNA(k), pGM(k), pATPC(k), pATTP(k), rCU_EA, ...
            pTO(k), pSat(k), pEff(k), pMFC(k), pTRF(k));
        % Phenotypic metrics
        fprintf(fid,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,', ...
            pRedox(k), pMFI(k), pAnab(k), ...
            pLac(k), pSucc(k), pAKG(k), ...
            pNADPH(k), pTCA(k), pLSat(k), pLUns(k), pLPL(k), pGln(k));
        % Subsystem activities
        sa_str = sprintf('%.6f,', pSA(k,:));
        fprintf(fid,'%s\n', sa_str(1:end-1));
    end
    fclose(fid);

    %% ── Store Pareto Front in MAT structure ──────────────────────────────
    ParetoSolutions.(matlab.lang.makeValidName(modelName)) = struct( ...
        'CU', pX(:,1), 'EA', pX(:,2), 'Biomass', pBio(:), 'ATP', pATP(:), ...
        'WarburgIndex', pWBI(:), 'NitrogenAnaplerosis', pNA(:), ...
        'GrowthMetabolism', pGM(:), 'ATPConsumption', pATPC(:), ...
        'ATPProduction', pATTP(:), 'RatioCU_EA', pR(:), ...
        'TotalTurnOver', pTO(:), 'MeanSaturation', pSat(:), ...
        'GrowthEfficiency', pEff(:), 'MeanFlowCentrality', pMFC(:), ...
        'TotalReversibleFlux', pTRF(:), ...
        'RedoxIndex', pRedox(:), 'MFI', pMFI(:), 'AnabolismScore', pAnab(:), ...
        'Oncomet_Lactate', pLac(:), 'Oncomet_Succinate', pSucc(:), ...
        'Oncomet_AlphaKG', pAKG(:), 'NADPHdemand', pNADPH(:), ...
        'TCA_completeness', pTCA(:), 'LipidSat', pLSat(:), ...
        'LipidUnsat', pLUns(:), 'LipidPL', pLPL(:), ...
        'GlnDependence', pGln(:), 'SubsystemActivity', pSA);

    %% ── Memory Management: Clear temporary model