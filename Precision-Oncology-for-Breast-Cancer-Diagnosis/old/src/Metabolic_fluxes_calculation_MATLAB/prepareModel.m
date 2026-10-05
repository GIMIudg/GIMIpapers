function preparedModel = prepareModel(model_path)

    biomassRxn = 'biomass_ob';
    tpHydrolysisRxn = 'ATP_HYDROLYSIS';
    newMet = 'M_cancer_biomass_c';
    
    preparedModel = readCbModel(model_path);
    
    % Add ATP reaction if it does not exist
    if ~ismember(tpHydrolysisRxn, preparedModel.rxns)
        preparedModel = addReaction(preparedModel, tpHydrolysisRxn, ...
            'reactionFormula', 'atp_c + h2o_c -> adp_c + pi_c + h_c', ...
            'lowerBound', -1000, 'upperBound', 1000);
    end
    
    % Add biomass metabolite if it does not exist
    if ~ismember(newMet, preparedModel.mets)
        preparedModel = addMetabolite(preparedModel, newMet, 'Cancer biomass intermediate');
    end
    
    % Verify that the biomass reaction exists
    if ~ismember(biomassRxn, preparedModel.rxns)
        error('Reaction %s was not found in model %s.', biomassRxn, model_path);
    end
    
    % Add metabolite to the biomass reaction
    metIdx = find(strcmp(preparedModel.mets, newMet));
    rxnIdx = find(strcmp(preparedModel.rxns, biomassRxn));
    preparedModel.S(metIdx, rxnIdx) = 1;
    
    % Add demand reaction if it does not exist
    demandRxnID = 'biomass_ex';
    if ~ismember(demandRxnID, preparedModel.rxns)
        metIndex = find(strcmp(preparedModel.mets, newMet));
        preparedModel.S(:, end+1) = 0;
        preparedModel.S(metIndex, end) = -1;
        preparedModel.rxns{end+1} = demandRxnID;
        preparedModel.rxnNames{end+1} = 'Export Cancer Biomass';
        preparedModel.lb(end+1) = 0;
        preparedModel.ub(end+1) = 1000;
        preparedModel.c(end+1) = 0;
        preparedModel.rules{end+1} = '';
        preparedModel.rxnGeneMat(end+1,:) = 0;
        preparedModel.grRules{end+1} = '';
    end
end
