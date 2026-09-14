modelsFolder = '/Users/eduardoruiz/Documents/MCBCI/MCBCI2/Sistemas metabólicos/Proyecto_Tesis/Modelos_actual';
outputFileTXT = 'Models_ids.txt'; % ⬅️ Name of the new text file
if exist(outputFileTXT,'file'), delete(outputFileTXT); end % ⬅️ Delete the .txt file if it already exists

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