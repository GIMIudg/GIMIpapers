% 1. Inicializar COBRA
cd('C:\Users\ami_0\cobratoolbox')
initCobraToolbox

%% FBA + norma 2 con pesos % este es el que tiene los pesos de german

%%
%carpetaModelos = 'C:\Users\ami_0\OneDrive\Documentos\Tareas\Tesis de maestria\fbaL2\modelosGERMAN';
%isfolder(carpetaModelos)                          % debe dar 1
%numel(dir(fullfile(carpetaModelos, '*.mat')))     % debe dar el número de modelos

%%
S = load('C:\Users\ami_0\OneDrive\Documentos\Tareas\Tesis de maestria\fbaL2\modelosGERMAN\GSM4852033.mat');
fieldnames(S)
%"C:\Users\ami_0\OneDrive\Documentos\Tareas\Tesis de maestria\fbaL2\modelosGERMAN\GSM4852033.mat"
model = S.specificModel;    % ajusta si tu variable se llama distinto


%%
global CBT_LP_SOLVER CBT_QP_SOLVER
CBT_LP_SOLVER
CBT_QP_SOLVER


%% 0. Configuracion
carpetaModelos = 'C:\Users\ami_0\OneDrive\Documentos\Tareas\Tesis de maestria\fbaL2\modelosGERMAN';   % carpeta con tus .mat (ajusta la ruta)
expressionRxns= 'expressionRxns' % VERIFICA este nombre en tu modelo

% La norma 2 es un problema cuadratico (QP): necesitas un solver QP.
% Descomenta el que tengas instalado:
% changeCobraSolver('gurobi', 'QP');
% changeCobraSolver('ibm_cplex', 'QP');
% changeCobraSolver('tomlab_cplex', 'QP');
% como ya esta pdco ya no descongelamos las de arriba

%% 1. Listar los modelos
archivos = dir(fullfile(carpetaModelos, '*.mat'));
n = numel(archivos);
fprintf('Se encontraron %d archivos .mat\n', n);
 

%% Tabla de resultados
nombre   = strings(n, 1);
%biomasa  = nan(n, 1);     % valor objetivo (f)
%norma2   = nan(n, 1);     % norma 2 de TODOS los flujos
%norma2Int = nan(n, 1);    % norma 2 solo de flujos internos
estado   = nan(n, 1);     % 1 = optimo
etiqueta = 'expresion';

% Guardamos por modelo: sus reacciones y sus flujos (los modelos pueden
% tener reacciones distintas, asi que juntamos todo despues)
rxnsCelda   = cell(n, 1);
flujosCelda = cell(n, 1);
nombresCelda = cell(n, 1);   % rxnNames
subsisCelda  = cell(n, 1);

%% 2. Bucle sobre modelos
for i = 1:n
    ruta = fullfile(carpetaModelos, archivos(i).name);
    nombre(i) = string(erase(archivos(i).name, '.mat'));

    % Cargar el .mat y localizar la variable del modelo
    S = load(ruta);
    campos = fieldnames(S);
   
    if isfield(S, 'specificModel')
        model = S.specificModel;
    elseif isfield(S, 'model')
        model = S.model;
    elseif numel(campos) == 1
        model = S.(campos{1});
    else
        warning('%s: no se encontro la variable del modelo, se omite.', nombre(i));
        continue
    end
    
     % Fijar la biomasa como objetivo
    %if ~ismember(rxnBiomasa, model.rxns)
     %   warning('%s: la reaccion %s no existe en el modelo, se omite.', nombre(i), rxnBiomasa);
      %  continue
    %end
    %model = changeObjective(model, rxnBiomasa);

    % Optimizar: maximiza biomasa y, entre soluciones empatadas,
    % elige la de menor norma 2 de los flujos internos

    % --- sin objetivo ---
    model.c(:) = 0;

    % --- penalizacion por expresion ---PESOS DE GERMANA
    expr = model.expressionRxns(:);
    expr(expr <= 0) = NaN;                       % evita log(0) = -Inf
    pen = -(log(expr) - max(log(expr)));
    pen(isnan(pen)) = mean(pen, 'all', 'omitnan');   % sin dato -> promedio
    pen(~model.SConsistentRxnBool(:)) = 0;           % intercambios -> 0
    model.expressionRxns2norm = pen;

    pesos = pen + abs(min(pen)) + 1e-6;
    sol = optimizeCbModel(model, 'min', pesos);      % <-- guarda el resultado

    estado(i) = sol.stat;
    if sol.stat ~= 1
        warning('%s: solucion no optima (stat = %d).', nombre(i), sol.stat);
        continue
    end
  

    % --- Guardar el vector de flujos (el "2" y el "3" de la suma) ---
    rxnsCelda{i}    = model.rxns;
    flujosCelda{i}  = sol.v;
    if isfield(model, 'rxnNames'),   nombresCelda{i} = model.rxnNames;   end
    if isfield(model, 'subSystems'), subsisCelda{i}  = model.subSystems; end

    %fprintf('%s | objetivo = %.6g | ||v||2 = %.4f | ||v_int||2 = %.4f\n', ...
     %   nombre(i), biomasa(i), norma2(i), norma2Int(i));
end

%% 3. Matriz de flujos: reacciones (union de todos los modelos) x muestras
validos = find(~cellfun(@isempty, rxnsCelda));
rxnUnion = {};
for k = validos'
    rxnUnion = union(rxnUnion, rxnsCelda{k}, 'stable');
end

V = nan(numel(rxnUnion), n);          % NaN = esa reaccion no existe en ese modelo
nombreRxn = repmat({''}, numel(rxnUnion), 1);
subsistema = repmat({''}, numel(rxnUnion), 1);

for k = validos'
    [~, loc] = ismember(rxnsCelda{k}, rxnUnion);
    V(loc, k) = flujosCelda{k};

    % Completar nombre y subsistema la primera vez que aparece cada reaccion
    if ~isempty(nombresCelda{k})
        vacio = cellfun(@isempty, nombreRxn(loc));
        nombreRxn(loc(vacio)) = nombresCelda{k}(vacio);
    end
    if ~isempty(subsisCelda{k})
        ss = subsisCelda{k};
        for q = 1:numel(ss)                       % subSystems puede ser celda de celdas
            if iscell(ss{q}) && ~isempty(ss{q}), ss{q} = ss{q}{1}; end
            if iscell(ss{q}), ss{q} = ''; end
        end
        vacio = cellfun(@isempty, subsistema(loc));
        subsistema(loc(vacio)) = ss(vacio);
    end
end

%% 4. Exportar
% 4a) Resumen (objetivo y normas)
resumenlvl = table(nombre, estado, ...
    'VariableNames', {'Muestra', 'Estado'});
disp(resumenlvl)
writetable(resumenlvl, fullfile(carpetaModelos, ['resumen_' etiqueta '.csv']));


% 4b) Flujos por reaccion (una columna por muestra)
nombresColumnas = matlab.lang.makeValidName(cellstr(nombre));
tablaFlujos = array2table(V, 'VariableNames', nombresColumnas);
tablaFlujos = addvars(tablaFlujos, rxnUnion, nombreRxn, subsistema, ...
    'Before', 1, 'NewVariableNames', {'Reaccion', 'Nombre', 'Subsistema'});
writetable(tablaFlujos, fullfile(carpetaModelos, ['flujos_' etiqueta '.csv']));

% lo podria cambiar a esta forma pero no recuerdo que era la etiqueta, buscar para cambiar a la otra forma 
% writetable(resumenlvl, fullfile(carpetaModelos, 'resumenlvl_norma2.csv'));

% 4c) Version .mat para seguir trabajando en MATLAB
save(fullfile(carpetaModelos, ['flujos_' etiqueta '.mat']), ...
    'V', 'rxnUnion', 'nombreRxn', 'subsistema', 'nombre', 'resumenlvl');

%% 5. Ejemplo: ver los flujos de una muestra y comprobar la norma
% j = 1;                                   % indice de la muestra (
% paciente o reaccion??)
% v = V(:, j);  v = v(~isnan(v));          % solo reacciones que existen en ese modelo
% fprintf('Norma 2 recalculada: %.4f | guardada: %.4f\n', norm(v, 2), norma2(j))
% % Las 10 reacciones con mayor flujo (en valor absoluto):
% [~, orden] = sort(abs(V(:, j)), 'descend', 'MissingPlacement', 'last');
% tablaFlujos(orden(1:10), {'Reaccion', 'Nombre', nombresColumnas{j}})
