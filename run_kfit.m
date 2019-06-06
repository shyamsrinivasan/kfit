% function to run k-fit - function form of kinetic estimate (capable of
% being called from a shell script)
function run_kfit(mi_files, ivaltype, ivalid)

if nargin < 2
    % feasible
    ivaltype = 1;
end
if nargin < 3
    % choice of ival_id
    ivalid = 1;
end

% collect model and relevant initial values from file
all_fnames = strsplit(mi_files, ',');
load(all_fnames{1}, 'model');

% load initial values from file
ivalfile = fopen(all_fnames{2});
% scan header
headers = textscan(ivalfile, '%s', 4, 'Delimiter', ',');
% scan values
ivals = textscan(ivalfile, '%f %f %f %s', 'Delimiter', ',');
fclose(ivalfile);

% collect all  initial values
initial = 1;
if ivaltype == 1 
    x0 = ivals{1}(ivals{2} == ivalid);  
    initial = 0;
elseif ivaltype == 2
    x0 = ivals{3}(ivals{2} == ivalid);
end

% run k-fit for selected initial value


% store result to file

   

