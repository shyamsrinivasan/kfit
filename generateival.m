% generate initial values and write them to file 
function generateival(model, nival, fname)

if nargin < 2
	nival = 1;
end

% computing total number of variables
ne = cell2mat(model.ensemble.ne);
nei = cell2mat(model.ensemble.nei);
nvr = cell2mat(model.ensemble.nvr);
nf = ne-1;
np = nf+nei+nvr;
nvar = sum(np);

% create nival random values
pos = 0;
inds = zeros(1,nvar);
vrind = inds;
for i = 1:length(ne)
    inds(1,pos+1:pos+nf(i)+nei(i)) = 1;
    vrind(1,pos+nf(i)+nei(i)+1:pos+np(i)) = model.d.flx{1}(i);
    pos = pos+np(i);
end
rng shuffle
xrand = rand(nvar, nival);
xrand(~inds, :) = 100*xrand(~inds, :);

% get feasible initial values for each random value
xfeas = cell(1, nival);
parfor i = 1:nival	
	x0 = xrand(:, i);
	opts = optimset('Display','iter','MaxFunEvals',10000000);
	xfeas{i} = fmincon(@(x) 1, x0, A, b, [], [], xlb, xub, [], opts);
end

% write all random and feasible values to file
colnames = {'feasible', 'id', 'random', 'status'};
tabinfo = cell(nvar * nival, 4);
for j = 1:nival
	tabinfo(nvar*(j-1) + 1:nvar*j, 1) = mat2cell(xfeas{i}, ones(1, nvar));	
	tabinfo(nvar*(j-1) + 1:nvar*j, 2) =...
	 mat2cell(repmat(j, nvar, 1);, ones(1, nvar));
	tabinfo(nvar*(j-1) + 1:nvar*j, 3) = mat2cell(xrand(:, i), ones(1, nvar));
end
tabinfo(:, 4) = {'optimal'};

table = cell2table(tabinfo, 'VariableNames', colnames);
writetable(table, fname, 'Delimiter', ',');




