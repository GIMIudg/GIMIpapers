% Generate iDopaNeuro models with different parameters using the
% function XomicsToMultipleModels. 
%
% The batch models are saved in 
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/multidimensionalModelGeneration


%% Set directory and specificData

clear

initCobraToolbox(false)
global CBTDIR

% specificData folder
dataFolder = uigetdir(pwd, 'Select the folder with the bibliomic data');

%filenames
bibliomicData = 'bibliomicData.xlsx';
files = dir([dataFolder filesep 'pythonResults' filesep 'processedData' filesep 'xomics_files' filesep '*.txt']);
transcriptomicData = {files.name};

% generic model
inputFolder = [CBTDIR filesep 'papers' filesep '2023_iDopaNeuro'];
genericModelName = 'Recon3DModel_301_xomics_input';

%% automatic processing from here on
%bibliomic data
specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

% specificData
for i = 1:length(transcriptomicData)
    rawName = regexprep(transcriptomicData{i}, '\.txt$', '');
    validFieldName = matlab.lang.makeValidName(rawName);
    modelGenerationConditions.specificData.(validFieldName) = specificData;
    modelGenerationConditions.specificData.(validFieldName).originalName = rawName;
    modelGenerationConditions.specificData.(validFieldName).transcriptomicData = readtable([dataFolder filesep 'transcriptomics' filesep transcriptomicData{i}]);
    modelGenerationConditions.specificData.(validFieldName).transcriptomicData.genes = string(modelGenerationConditions.specificData.(validFieldName).transcriptomicData.genes);
end

modelGenerationConditions.outputDir = [dataFolder filesep 'matlabResults' filesep 'models'];
modelGenerationConditions.cobraSolver = {'gurobi'}; 

%load generic model
load([inputFolder filesep genericModelName])
modelGenerationConditions.genericModel.model = model;


%% Fixed options
param.setObjective = ''; % No objective function
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.findThermoConsistentFluxSubset = 1;


%% Create models

% Remove expressionRxns
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Create models
directoriesWithModels = XomicsToMultipleModels(modelGenerationConditions, param);

disp('Models saved in:')
disp([dataFolder filesep 'models'])





%%

% 1. Definir directorios
baseDir = pwd
modelsDir =  'C:\Users\dange\Music\GIMI\01_posgrado\MCBCI_analisisCancer\resultados\newModels\models';
outputDir = fullfile(baseDir, 'models_renamed'); % Carpeta donde se moverán/guardarán los modelos

% Crear carpeta de destino si no existe
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Se creó la carpeta de destino: %s\n', outputDir);
end

% 2. Obtener lista de subcarpetas dentro de ./models
dirItems = dir(modelsDir);
% Filtrar solo directorios y excluir '.' y '..'
isSubFolder = [dirItems.isdir] & ~ismember({dirItems.name}, {'.', '..'});
modelFolders = dirItems(isSubFolder);

if isempty(modelFolders)
    disp('No se encontraron subcarpetas dentro del directorio ./models');
    return;
end

% 3. Recorrer carpeta por carpeta
for k = 1:length(modelFolders)
    folderName = modelFolders(k).name;
    modelFilePath = fullfile(modelsDir, folderName, 'Model.mat');

    fprintf('Procesando (%d/%d): Carpeta "%s"...\n', k, length(modelFolders), folderName);

    % Verificar si existe el archivo Model.mat en la subcarpeta
    if ~isfile(modelFilePath)
        warning('No se encontró "Model.mat" en %s. Omitiendo.', folderName);
        continue;
    end

    try
        % Cargar datos del archivo .mat
        matData = load(modelFilePath);

        % Verificar que contenga el struct 'Model'
        if isfield(matData, 'Model')
            Model = matData.Model;

            % Verificar existencia del campo con el nombre original
            if isfield(Model, 'XomicsToModelSpecificData') && ...
               isfield(Model.XomicsToModelSpecificData, 'originalName')

                originalName = Model.XomicsToModelSpecificData.originalName;

                % Limpiar la extensión si ya la traía el texto
                [~, cleanBaseName, ~] = fileparts(originalName);

                % Renombrar la variable a 'specificModel'
                specificModel = Model;

                % Definir la ruta completa del nuevo archivo .mat
                newMatPath = fullfile(outputDir, [cleanBaseName, '.mat']);

                % Guardar únicamente la variable 'specificModel'
                save(newMatPath, 'specificModel', '-v7.3');

                fprintf(' -> Exito: Guardado en "%s" como %s.mat\n', outputDir, cleanBaseName);
            else
                warning('El modelo en "%s" no contiene "XomicsToModelSpecificData.originalName".', folderName);
            end
        else
            warning('El archivo en "%s" no contiene la variable "Model".', folderName);
        end

    catch ME
        warning('Error al procesar el archivo en "%s": %s', folderName, ME.message);
    end
end

disp('¡Proceso completado exitosamente!');


% 
%%
% % 1. Definir el directorio objetivo (por defecto la carpeta del script actual)
% targetDir = 'C:\Users\dange\Music\GIMI\01_posgrado\MCBCI_analisisCancer\resultados\newModels\models'; % o especifica una ruta como 'C:\ruta\a\carpeta'
% 
% % Obtener la lista de subcarpetas en el directorio objetivo
% dirItems = dir(targetDir);
% isSubFolder = [dirItems.isdir] & ~ismember({dirItems.name}, {'.', '..'});
% subFolders = dirItems(isSubFolder);
% 
% if isempty(subFolders)
%     disp('No se encontraron subcarpetas para evaluar.');
%     return;
% end
% 
% fprintf('Evaluando %d subcarpetas en: %s\n\n', length(subFolders), targetDir);
% 
% % 2. Recorrer cada subcarpeta
% for k = 1:length(subFolders)
%     folderName = subFolders(k).name;
%     folderPath = fullfile(subFolders(k).folder, folderName);
% 
%     % Buscar archivos .mat dentro de la subcarpeta (incluyendo subcarpetas internas con '**')
%     matFiles = dir(fullfile(folderPath, '**', '*.mat'));
% 
%     if isempty(matFiles)
%         fprintf('[ELIMINANDO] "%s" no contiene archivos .mat.\n', folderName);
%         try
%             % La opción 's' elimina la carpeta y todo su contenido no-.mat que tenga dentro
%             rmdir(folderPath, 's');
%         catch ME
%             warning('No se pudo eliminar la carpeta "%s": %s', folderName, ME.message);
%         end
%     else
%         fprintf('[CONSERVADA] "%s" contiene %d archivo(s) .mat.\n', folderName, length(matFiles));
%     end
% end
% 
% disp(' ');
% disp('¡Proceso de limpieza completado!');











