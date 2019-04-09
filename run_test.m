% run script for original k-fit (Sarat's matlab code)

% empty cell array for storing results
res={};

% get initial time
cpu_init=cputime;

% formulate model
model = modcompile(strcat(pwd,'data/toy_model.xlsx'),strcat(pwd,'/toy_mechanism.xlsx'),strcat(pwd,'/toy_data.xlsx'));

% get kinetic parameter estimates
res = kineticestimate(model, res);

run_time=cputime-cpu_init;

% save results (matlab workspace)
save('toy_model.mat','res','model','run_time');