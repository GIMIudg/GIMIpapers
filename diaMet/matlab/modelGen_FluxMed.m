% modelGen_FluxMed.m
% Genera modelos de islote pancreatico / celula beta humana (uno por muestra
% transcriptomica) con XomicsToMultipleModels -> XomicsToModel.
%
% Entradas (dentro de la carpeta del proyecto MCBCI_DiaMED):
%   datos/bibliomicData_pancreas.xlsx      bibliomica curada (islote/celula beta, IDs Recon3D)
%   resultados/transcriptomicData/*.txt    130 muestras de islotes, formato "genes,expVal"
%                                          (Entrez, TPM lineal). Se generan con
%                                          diaMet/python/csv2txt_transcriptomics.py a partir de
%                                          datos/DE_DONDE_SALIO_ISLOTESPANCREATICOS/*.csv
%
% Salidas:
%   resultados/modelosFluxMed/models/<muestra>/Model.mat   (un modelo por muestra)
%   resultados/modelosFluxMed/models_renamed/<muestra>.mat (variable specificModel)
%
% Requiere COBRA Toolbox (XomicsToModel) y Gurobi.

%% 0. Configuracion
clear
initCobraToolbox(false)
global CBTDIR

% Carpeta del proyecto (MCBCI_DiaMED). Cambiar por la ruta fija si se prefiere:
% projectDir = '~/work/Drive/GIMI/01_posgrado/MCBCI_DiaMED';
projectDir = uigetdir(pwd, 'Selecciona la carpeta MCBCI_DiaMED');
if isequal(projectDir, 0)
    error('No se selecciono la carpeta del proyecto.')
end

bibliomicFile     = fullfile(projectDir, 'datos', 'bibliomicData_pancreas.xlsx');
transcriptomicDir = fullfile(projectDir, 'resultados', 'transcriptomicData');
outputDir         = fullfile(projectDir, 'resultados', 'modelosFluxMed');

% Subconjunto opcional de muestras (vacio = las 130). Ej.: {'GSM4852033','GSM4852034'}
samplesToRun = {};

% Modelo generico
inputFolder      = fullfile(CBTDIR, 'papers', '2023_iDopaNeuro');
genericModelName = 'Recon3DModel_301_xomics_input.mat';

%% 1. Comprobaciones de entrada
if ~isfile(bibliomicFile)
    error('No existe %s', bibliomicFile)
end
files = dir(fullfile(transcriptomicDir, 'GSM*.txt'));
if isempty(files)
    error(['No hay archivos .txt en %s. Ejecuta primero:\n' ...
        '  python3 diaMet/python/csv2txt_transcriptomics.py'], transcriptomicDir)
end
transcriptomicData = {files.name};
if ~isempty(samplesToRun)
    transcriptomicData = transcriptomicData(ismember(regexprep(transcriptomicData, '\.txt$', ''), samplesToRun));
end
fprintf('%d muestras transcriptomicas encontradas.\n', numel(transcriptomicData))

%% 2. Modelo generico (Recon3D)
load(fullfile(inputFolder, genericModelName), 'model')
% Quitar expressionRxns ANTES de asignarlo a modelGenerationConditions
% (en modelGen.m se quitaba despues y no tenia efecto)
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end
modelGenerationConditions.genericModel.model = model;

%% 3. Datos bibliomicos curados (islote/celula beta)
specificData = preprocessingOmicsModel(bibliomicFile, 1, 1);

% Verificacion rapida de IDs contra Recon3D (las reacciones de rxns2add aun no existen en el generico)
newRxns = specificData.rxns2add.rxns;
newMets = {'zn2[g]', 'insulin[g]', 'inshex[g]', 'cpeptide[g]'};
missingRxns = setdiff([specificData.activeReactions; specificData.rxns2remove.rxns; specificData.rxns2constrain.rxns], ...
    [model.rxns; newRxns]);
missingMets = setdiff(specificData.presentMetabolites.mets, [model.mets; newMets']);
if ~isempty(missingRxns) || ~isempty(missingMets)
    disp(missingRxns); disp(missingMets)
    error('Hay IDs en la bibliomica que no estan en Recon3D ni en rxns2add.')
end

%% 4. Datos transcriptomicos: una condicion (specificData) por muestra
for i = 1:numel(transcriptomicData)
    rawName = regexprep(transcriptomicData{i}, '\.txt$', '');
    validFieldName = matlab.lang.makeValidName(rawName);

    T = readtable(fullfile(transcriptomicDir, transcriptomicData{i}), ...
        'Delimiter', ',', 'ReadVariableNames', true, 'FileType', 'text');
    if ~all(ismember({'genes', 'expVal'}, T.Properties.VariableNames))
        error('%s no tiene las columnas genes,expVal', transcriptomicData{i})
    end
    if any(T.expVal < 0) || any(isnan(T.expVal))
        error('%s contiene valores negativos o NaN', transcriptomicData{i})
    end
    T.genes = string(T.genes);

    modelGenerationConditions.specificData.(validFieldName) = specificData;
    modelGenerationConditions.specificData.(validFieldName).originalName = rawName;
    modelGenerationConditions.specificData.(validFieldName).transcriptomicData = T;
end

%% 5. Condiciones de generacion
modelGenerationConditions.outputDir   = fullfile(outputDir, 'models');
modelGenerationConditions.cobraSolver = {'gurobi'};
modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'};
modelGenerationConditions.activeGenesApproach  = {'oneRxnPerActiveGene'};

% expVal esta en TPM lineal; XomicsToModel compara log(expVal) (logaritmo natural)
% contra el umbral. 0 -> activo si TPM >= 1. (XomicsToMultipleModels usa 2 por
% defecto, es decir TPM >= 7.4, demasiado estricto porque INS/GCG acaparan los TPM.)
modelGenerationConditions.transcriptomicThreshold = 0;
% No borrar genes por baja expresion (PC, GCK, KCNJ11 tienen TPM bajos en islote completo)
modelGenerationConditions.inactiveGenesTranscriptomics = false;
% La curacion beta-celular prevalece sobre el RNA-seq de islote (contaminacion acinar)
modelGenerationConditions.curationOverOmics  = true;
% Las reacciones activas curadas (p.ej. PYRt2mMPC, INSsyn) no se eliminan por genes inactivos
modelGenerationConditions.activeOverInactive = true;
modelGenerationConditions.closeIons = false;

%% 6. Parametros fijos
param.setObjective = '';            % sin funcion objetivo
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
param.sinkDMinactive = 1;           % sinks y demandas no core inactivas
param.addCoupledRxns = 1;           % (no hay hoja coupledRxns para islote; no tiene efecto)
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = false;         % no hay mediaData/cellCultureData para islotes
param.findThermoConsistentFluxSubset = 1;
param.printLevel = 1;

%% 7. Generar modelos
directoriesWithModels = XomicsToMultipleModels(modelGenerationConditions, param);
disp('Modelos guardados en:')
disp(modelGenerationConditions.outputDir)

%% 8. Copiar cada Model.mat como <muestra>.mat (variable specificModel)
modelsDir = modelGenerationConditions.outputDir;
renamedDir = fullfile(outputDir, 'models_renamed');
if ~exist(renamedDir, 'dir')
    mkdir(renamedDir);
end

dirItems = dir(modelsDir);
modelFolders = dirItems([dirItems.isdir] & ~ismember({dirItems.name}, {'.', '..'}));

for k = 1:numel(modelFolders)
    folderName = modelFolders(k).name;
    modelFilePath = fullfile(modelsDir, folderName, 'Model.mat');
    if ~isfile(modelFilePath)
        warning('No se encontro Model.mat en %s. Se omite.', folderName);
        continue
    end
    try
        matData = load(modelFilePath, 'Model');
        Model = matData.Model;
        if isfield(Model, 'XomicsToModelSpecificData') && isfield(Model.XomicsToModelSpecificData, 'originalName')
            [~, cleanBaseName] = fileparts(Model.XomicsToModelSpecificData.originalName);
        else
            cleanBaseName = folderName;
        end
        specificModel = Model;
        save(fullfile(renamedDir, [cleanBaseName '.mat']), 'specificModel', '-v7.3');
        fprintf(' -> %s.mat (%d rxns, insulina: %d)\n', cleanBaseName, numel(Model.rxns), ...
            any(strcmp(Model.rxns, 'EX_inshex[e]')));
    catch ME
        warning('Error al procesar %s: %s', folderName, ME.message);
    end
end
disp('Proceso completado.')
