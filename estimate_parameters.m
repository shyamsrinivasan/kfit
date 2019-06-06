% estimate parameters through k-fit for a given initial value x0 for a give
% model 'model'
function results = estimate_parameters(model, x0, initial, reinit)

if nargin < 3
    initial = 1;
end

if nargin < 4
    reinit = 0;
end

fmin = Inf;
model.options.dfbase = 1e-6;

itime = tic;

% problem setup for calculating active constraints and
% get feasible initial value from random x0 input
[x0, A, b, actcon] = initialize(x0, model, initial, reinit);

% run k-fit with lsqsolve
[x, f, ~, actcon] = lsqsolve(x0, model, A, b, actcon);
if f < fmin
    fmin = f;
    xinit = x;
end

% re-estimation w/ machine precision dfbase
model.options.dfbase = eps;
iter = 0;
fail = true;
x = xint;
xopt = xinit;

while fail || iter <= 1
    [x, f, fail, actcon] = lsqsolve(x, model, A, b, actcon);
    if f < fmin
        fmin = f;
        xopt = x;
    end
    iter = iter + 1;
end

etime = toc(itime);

% compile results and write to file
result = collectresults(model, xopt, etime);


fprintf('\n End time = %f', toc(itime));
end

function [x0, A, b, actcon] = initialize(x0, model, initial, reinit)

% lp setup
[A, b, lb, ub] = constraints(model);

if initial && ~reinit
    % get feasible values
    x0 = fmincon(@(x) 1, x0, A, b, [], [], lb, ub, [],...
                 optimset('Display','iter','MaxFunEvals',10000000));
end

% get active constraints
A = -A;
b = -b;
nx = length(x0);
A = [A;eye(nx);-eye(nx)];
b = [b;lb;-ub];
actcon = find(A*x0<=b);

end




