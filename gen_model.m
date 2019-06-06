% function that can be run as a job to construct model and generate initial
% values and store them to file
function model = gen_model(file_names, nival)

% model_files - comma seperated list of model, mechanism and data file
all_files = strsplit(file_names, ',');

if length(all_files) < 4
    all_files{4} = fullfile(pwd, 'models', 'kfitmodel');
end
if length(all_files) < 5
    all_files{5} = fullfile(pwd, 'initial_values', 'init_val');
end

% construct matlab model and save model as .mat file
model = modcompile(all_files{1}, all_files{2}, all_files{3});
save(all_files{4}, 'model');

% generate rnadom and feasible initial values and save them to file
generateival(model, nival, all_files{5});

