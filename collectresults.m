% mimics original compileresult.m and integrates elements from
% pycompileresult.m
function res = collectresults(model, xopt, time)

nr = length(model.rid);
ne = cell2mat(model.ensemble.ne);
nmet = length(model.metprop);

res.fluxid = cell(nr, 1);
res.metid = cell(nmet, 1);

res.kf = cell(nr, 1);
res.kr = cell(nr, 1);
res.kfstd = cell(nr, 1);
res.krstd = cell(nr, 1);
res.Ki = cell(nr, 1);
res.Kistd = cell(nr, 1);

res.ef = cell(nr, 1);

% calculate fluxes, concentrations, objective func value and kinetic
% parameters
% fmin
[r,W,J,v,css] = rescalc(xopt,model);
res.fmin = r'*W*r;
H = J'*W*J;

% kinetic parameters and reversibilities
[k,Jk,rev] = calc_k(model,xopt);

% Uncertainties in kinetic parameter estimation
H = (H+H')/2;
s = svds(H,1);
tol = max(size(H,1),size(Jk,1))*eps(s);
Hinv = pinv(H,tol);
Covmat = Jk*Hinv*Jk';
k_sd = sqrt(abs(diag(Covmat)));

ncond = length(v(1,:));

res.flux = cell(nr, ncond);
res.vpert = model.d.pert;
res.conc = zeros(nmet, ncond);

% data vs model comparison measures
res.mrxnid = cell(nr, ncond);
res.flxmodel = zeros(nr, ncond);
res.flxdata = zeros(nr, ncond);
res.flxWRES = zeros(nr, ncond);
res.flxSRES = zeros(nr, ncond);
res.mutant = cell(1, ncond);

% data vs model reversibility
nrevs = length(model.d.ridx);	
res.revid = cell(nrevs, 1);
res.revmodel = zeros(nrevs, 1);
res.revdata = zeros(nrevs, 1);
res.revWRES = zeros(nrevs, 1);
res.revSRES = zeros(nrevs, 1);

% information for reinitialization of kfit
res.x = [];

% standard deviations ?
res.x_sd = [];



kin = struct('id','','kf',[],'kr',[],'Ki',[],'kact',[],'kf_sd',[],'kr_sd',[],'Ki_sd',[],'kact_sd',[]);


[res.kinetic_params(1:nr)] = deal(kin);
%x1 = xopt;
%nu = model.ensemble.nu;
%u = x1(1:nu);
%x1(1:nu) = [];

%nei = cell2mat(model.ensemble.nei);
%nvr = cell2mat(model.ensemble.nvr);
%nf = ne-1;
%np = nf+nei+nvr;
for i = 1:nr
    res.fluxid{i} = model.rid(i);
    res.kinetic_params(i).id = model.rid(i);
    
    eind = 1:ne(i);
    krevind = 2*eind;    
    kfwdind = krevind-1;
    
    krxn = k(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    k_sd_rxn = k_sd(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    
    kcatal = krxn(1:2*ne(i));
    k_sd_catal = k_sd_rxn(1:2*ne(i));
    
    res.kf{i} = kcatal(kfwdind);
    res.kr{i} = kcatal(krevind);    
    res.kfstd{i} = k_sd_catal(kfwdind);
    res.krstd{i} = k_sd_catal(krevind);    
    
    krxn(1:2*ne(i)) = [];
    k_sd_rxn(1:2*ne(i)) = [];
    
    res.Ki{i}= krxn;
    res.Kistd = k_sd_rxn;
end

for i = 1:ncond
    pert = model.d.vpert(:,i);
    if ~any(~pert)
        res.mutant{i} = {'WT'};
    else
%         mut_id = model.rid(~pert);
        res.mutant{i} = model.rid(~pert);
    end
end

% Predicted fluxes and concentrations
flx = struct('id','','val',[]);
met = struct('id','','val',[]);
mut = struct('gene_KO','','fluxes',[],'concentrations',[]);


[m1(1:ncond)] = deal(mut);
for i = 1:ncond
    %m1 = mut;
    pert = model.d.vpert(:,i);
    if ~any(~pert)
        m1(i).gene_KO = {'WT'};
    else
        m1(i).gene_KO = model.rid(~pert);
    end
    [f(1:nr)] = deal(flx);
    for j = 1:nr
        f(j).id = model.rid(j);
        f(j).val = v(j,i);
    end
    [m(1:nmet)] = deal(met);
    for j = 1:nmet
        m(j).id = model.metprop(j).metid;
        m(j).val = css(j,i);
    end
    m1(i).fluxes = f;
    m1(i).concentrations = m;
end
res.predictions = m1;

% Lack-of-Fit

rever = struct('flxid','','val',[],'data',[],'WRES',[],'SRES',[]);
%flxft = struct('flxid','','val',[],'data',[],'WRES',[],'SRES',[]);

%reversibilities
nrevs = length(model.d.ridx);
[rx(1:nrevs)] = deal(rever);
for i = 1:nrevs
    rx(i).flxid = model.rid(model.d.ridx(i));
    rx(i).val = rev(model.d.ridx(i));
    rx(i).data = model.d.revs(i);
    rx(i).WRES = (rx(i).val-rx(i).data)/model.d.rerr(i);
    rx(i).SRES = (rx(i).WRES).^2;
end
res.residuals.reversibility = rx;

%fluxes
flx = struct('measid','','combination','','mut','','data',[],'val',[],'WRES',[],'SRES',[]);
%nflx = length(model.d.id);
nflx = 0;
for i = 1:ncond
    nflx = nflx+length(model.d.id{i});
end
[flxs(1:nflx)] = deal(flx);
ctr = 0;
%{
idx = model.d.idx;
fdat = model.d.flx;
erdat = model.d.err;
for i = 1:ncond
    relid = idx(idx<=length(model.rid));
    fd1 = fdat(relid);
    err1 = erdat(relid);
    pert = model.d.vpert(:,i);
    for j = 1:length(relid)
        ctr = ctr+1;
        if ~any(~pert)
            flxs(ctr).mut = {'WT'};
        else
            flxs(ctr).mut = model.rid(~pert);
        end
        flxs(ctr).measid = model.rid(relid(j));
        flxs(ctr).combination = model.rid(relid(j));
        flxs(ctr).val = v(relid(j),i);
        flxs(ctr).data = fd1(j);
        flxs(ctr).WRES = (flxs(ctr).val - flxs(ctr).data)/err1(j);
        flxs(ctr).SRES = (flxs(ctr).WRES).^2;
    end
    fdat(idx<=length(model.rid)) = [];
    erdat(idx<=length(model.rid)) = [];
    idx(idx<=length(model.rid)) = [];
    idx = idx-length(model.rid);
end
%}
for i = 1:ncond
    pert = model.d.vpert(:,i);
    for j = 1:length(model.d.id{i})
        ctr = ctr+1;
        if ~any(~pert)
            flxs(ctr).mut = {'WT'};
        else
            flxs(ctr).mut = model.rid(~pert);
        end
        flxs(ctr).measid = model.d.id{i}(j);
        flxs(ctr).combination = model.rid(find(model.d.rmap{i}(j,:)));
        flxs(ctr).val = model.d.rmap{i}(j,:)*v(:,i);
        flxs(ctr).data = model.d.flx{i}(j);
        flxs(ctr).WRES = (flxs(ctr).val - flxs(ctr).data)/model.d.err{i}(j);
        flxs(ctr).SRES = (flxs(ctr).WRES).^2;
    end
end


res.residuals.fluxes = flxs;        


%reinitialization
res.reinit_data = xopt;
res.xbest_sd = sqrt(abs(diag(Hinv)));





end