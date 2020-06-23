function [axon_collection, model, zoomed_model, output] = createOne2DWMModel(model_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Input required
%  axon_dictionary_path   : input axon dictionary path (provided in the
%  data folder)
 
% %%%%%%%%%% Input optional 
% model_params is a structure which contain all the algorithm options describe below
% The given values are the one set by defauts
%
% %%%% White matter model 
% model_params.number_of_axons = 400;
% % number of axons put on the original grid
% model_params.dims = estimated from the number of axons as well as the axon size 
% % dims of the original axon grid 
% model_params.mask = estimated from the original fiber volume fraction (FVF)
% % Represents the final WM model area after axon packing where the actual
% FVF, g-ratio are estimated
% 
% %%%% Axons packing 
% model_params.max_FVF = 0.85;
% % FVF value where the axon packing is stopped (it cannot be much higher than
% % 0.85 )
% model_params.max_iteration = 5000;
% % Max number of iteration of the axon packing (stop the axon packing if it
% % cannot reach the max_FVF)
% model_params.packing_speed = 0.5;
% % The packing speed weights the attraction/repulsion of the axons. A higher
% % value accelerates the packing process but can create axons overlap
% 
% %%%% Axons dispersion
% model_params.expected_FVF = 0.7; 
% %  FVF of the final WM model after remove or spread the axons
% model_params.dispersion_mode = 'spread'; 
% % Dispersion mode can be 
% % - remove where the axons are randomly remove from the packed model
% % - spread where the axons are spread from the grid center
% model_params.tolerance = 0.001;
% % Tolerance between expected FVF and actual FVF of the model
% 
% %%%% Change g-ratio
% model_params.expected_g_ratio = 0.6;
% % Change the myelin thickness to reach an expected g-ratio
% 
% %%%% Plot / save model
% model_params.plot_model = 1;
% model_params.save_model = '/project/3015069.04/WM_Models/toto.mat';
%
% %%%%%%%%%% Outputs
% %%%% axon_collection is a structure of the white matter model where each
% element represent an axon with one required field
% % - data which corresponds to the myelin sheath
% % - Centroid (required to run this function but not to simulate field
% perturbation)
% optional fields: gRatio, axonEquiDiameter, myelinThickness
%
% %%%% Model is the final WM model where 0 is extra axonal, 1 myelin, 0.5
% intra axonal
% %%%% Zoomed model is the model within the mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default parameters  
model_params = check_and_set_default(model_params);
% verbose = model_params.verbose;

if ~islogical(model_params.save_model)
    if  isfield(model_params, 'save_model')
        % makes sure directory exists
        a=max([find(model_params.save_model=='/'), find(model_params.save_model=='\')]);
        if ~isfolder(model_params.save_model(1:a))
            mkdir(model_params.save_model(1:a));
        end    
    end
end

%%% added by KC, 20200621
% for reproducibility
rng(model_params.seed)

if model_params.DEBUG
    disp('DEBUG mode on')
    save('DEBUG_tmp_output','model_params');
end
%%%

% Load pre-exist dictionary
% default: data/axonMediumDict.mat
fprintf('Load axon dictionary ...')
% load(axon_dictionary_path);
load(model_params.dictionary_fn) 
fprintf('done\n')

fprintf('Randomly select %d axon shapes ...',model_params.number_of_axons);
list_axons = randi(length(axonDico), model_params.number_of_axons, 1);
original_axon_collection = axonDico(list_axons);
for k = 1:length(original_axon_collection)
    original_axon_collection(k).data = double(original_axon_collection(k).data);
end
fprintf('done\n')

% Setup axons grid
fprintf('Setup axons on a grid ...')
if ~isfield(model_params, 'dims')
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection);
else
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection, model_params.dims);
end
fprintf('done\n')

if ~isfield(model_params, 'mask')
    model_params.mask = createAdaptedMask(axon_collection, dims);
end

% initial state
[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, model_params.mask, model_params.plot_model);

% Axon packing
fprintf('Process packing ...')
[axon_collection, FVF_packed_model] = packAxons(axon_collection, model_params.mask, model_params.max_iteration, model_params.max_FVF, model_params.packing_speed, model_params.plot_model);
fprintf('done\n')
disp(['FVF packed model : ' num2str(FVF_packed_model)]);

% Axon dispersion (optional)
if isfield(model_params, 'expected_FVF')
    if  strcmpi(model_params.dispersion_mode,'remove')
        disp('Remove axons ...')
        [axon_collection, FVF] = removeAxons(axon_collection, model_params.expected_FVF, model_params.tolerance, model_params.mask, model_params.plot_model);
    elseif strcmpi(model_params.dispersion_mode,'spread')
        disp('Spread axons ...')
        [axon_collection, FVF] = repulseAxons(axon_collection, model_params.expected_FVF, model_params.tolerance, model_params.mask, model_params.plot_model);
    else
        error('dispersion mode should be remove or spread');
    end
    
    disp(['current FVF : ' num2str(FVF)]);
    disp('done')
end
% Avoid axon overlap
axon_collection = avoidAxonOverlap(axon_collection, dims);
axon_collection = convertAxonDataToRoundValues(axon_collection);

% Change g-ratio (optional)
if isfield(model_params, 'expected_g_ratio')
    fprintf('Change g ratio ...')
    axon_collection = changeGRatio(axon_collection, model_params.expected_g_ratio, model_params.mask);
    fprintf('done\n')
end

% End state
[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, model_params.mask, model_params.plot_model);

% Save model
mask = model_params.mask;
if model_params.save_model
    fprintf('Save model ...')
    disp(model_params.save_model);
    try
        save(model_params.save_model, 'model', 'zoomed_model', 'FVF', 'g_ratio', 'axon_collection', 'dims', 'mask', 'model_params')
    catch
        disp ('failed to save...')
    end    
    fprintf('done\n')
end

%%% KC added 20200621
output.model        = model;
output.zoomed_model = zoomed_model;
output.FVF          = FVF;
output.g_ratio      = g_ratio;
output.model_params = model_params;
output.mask = mask;
output.dims = dims;
%%%

end

%% parser for essential input
function model_param2 = check_and_set_default(model_param)

model_param2 = model_param; 

curr_dir        = fileparts(mfilename('fullpath'));
Whist_HOME      = curr_dir(1:end-13);
Whist_data_dir  = fullfile(Whist_HOME,'data');

dictionary_fn = fullfile(Whist_data_dir,'axonMediumDict.mat');

try model_param2.max_FVF            = model_param.max_FVF;        	catch; model_param2.max_FVF = 0.8; end
try model_param2.max_iteration      = model_param.max_iteration;   	catch; model_param2.max_FVF = 5000; end
try model_param2.packing_speed      = model_param.packing_speed;  	catch; model_param2.packing_speed = 0.5; end
try model_param2.dispersion_mode    = model_param.dispersion_mode;	catch; model_param2.dispersion_mode = 'spread'; end
try model_param2.number_of_axons    = model_param.number_of_axons;  catch; model_param2.number_of_axons = 400; end
try model_param2.plot_model         = model_param.plot_model;     	catch; model_param2.plot_model = true; end
try model_param2.tolerance          = model_param.tolerance;      	catch; model_param2.tolerance = 0.1; end
try model_param2.seed               = model_param.seed;             catch; model_param2.seed = rng; end
% for faster accessibility
try model_param2.dictionary_fn   	= model_param.dictionary_fn; 	catch; model_param2.dictionary_fn = dictionary_fn; end

try model_param2.save_model       	= model_param.save_model;     	catch; model_param2.save_model = false; end
try model_param2.verbose            = model_param.verbose;          catch; model_param2.verbose = true; end
try model_param2.DEBUG              = model_param.DEBUG;            catch; model_param2.DEBUG = false; end


end