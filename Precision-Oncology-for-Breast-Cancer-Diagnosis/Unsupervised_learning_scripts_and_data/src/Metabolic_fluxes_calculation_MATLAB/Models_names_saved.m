% Resolve paths relative to this script's location
% Hierarchy: Metabolic_fluxes_calculation_MATLAB/ -> src/ -> Unsupervised_learning_scripts_and_data/
scriptDir    = fileparts(mfilename('fullpath'));
dataRoot     = fullfile(scriptDir, '..', '..', 'Clinical_data_and_models_ids');
modelsFolder = fullfile(dataRoot, 'All_models_created');
outputDir    = fullfile(dataRoot, 'GEMs_Data_for_construction');
outputFileTXT = fullfile(outputDir, 'Model''s_ids.txt'); % Output text file path
if exist(outputFileTXT,'file'), delete(outputFileTXT); end % Delete if already exists

%% ===========================
initCobraToolbox(false);
changeCobraSolver('gurobi','all');
modelFiles = dir(fullfile(modelsFolder,'*.mat'));
numModels = length(modelFiles);

% 1. OPEN THE TEXT FILE FOR WRITING (outside the loop)
fileID = fopen(outputFileTXT, 'w'); % 'w' means write mode (creates or overwrites)

for modelIdx = 1:numModels
    modelName = modelFiles(modelIdx).name;
    
    fprintf('Processing model %d/%d: %s\n', modelIdx, numModels, modelName);
    
    % 2. WRITE MODEL NAME TO THE TXT FILE
    % Use '\n' to ensure each name is on a new line.
    fprintf(fileID, '%s\n', modelName); 
end

% 3. CLOSE THE FILE (IMPORTANT)
fclose(fileID); 

disp(['✅ Names of the ' num2str(numModels) ' models saved to ' outputFileTXT]);