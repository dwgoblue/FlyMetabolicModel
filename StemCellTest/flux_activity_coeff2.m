function [dynamicmodel,grate_wt, geneko_growthrate, rxnko_growthrate,slope,solverobj,geneko_growthrate_obj,rxnko_growthrate_obj] = flux_activity_coeff2(model, timecourse_metabolomics_datafile,sheetname,kappa,kappa2,genedelflag,rxndelflag)

if (~exist('kappa','var')) || (isempty(kappa))
    kappa = 1;
end

if (~exist('kappa2','var')) || (isempty(kappa2))
    kappa2 = 1E-3;
end

if (~exist('sheetname','var')) || (isempty(sheetname))
    sheetname = 'Sheet1';
end

if (~exist('genedelflag','var')) || (isempty(genedelflag))
    genedelflag = 0;
end

if (~exist('rxndelflag','var')) || (isempty(rxndelflag))
    rxndelflag = 0;
end



model_rpmi = model;

%kappa = 10;
%kappa = 1.5;

model1 = model_rpmi;
model = model1;
model.A = model1.S;
model.obj = model1.c;
model.rhs = model1.b;
model.sense =repmat( '=',[size(model1.S,1),1]);
model.lb = model1.lb;
model.ub = model1.ub;
model.vtype = repmat('C',size(model1.S,2),1);
model.modelsense = 'max';
m2 = model;

% load metabolomics data
[num txt] = xlsread(timecourse_metabolomics_datafile,sheetname);
manual_matchpos = num(2:end,1:3); % coresponding position in model
%timevec =  num(1, 4:end); % time points
timevec =  num(1, 4:9); % time points
%maty = num(2:end,4:end); % metabolomics data
maty = num(2:end,4:9); % metabolomics data
maty = knnimpute(maty);
slope = zeros(size(maty,1),2);
for i = 1:size(maty,1)
    m1 = mean(maty(i,:));
    [p S] = polyfit(timevec, maty(i,:),1);
    slope(i,:) = p(1)/p(2);
    %slope(i,:) = p(1);
end
% normalize the data so that the F.A.C is a small number less than 1
%slope(:,1) = slope(:,1)./max(abs(slope(:,1))).*1; 
%slope(:,1) = slope(:,1);
%slope(1:154,1) = slope(1:154,1)./max(abs(slope(1:154,1))).*1; 
%slope(155:326,1) = slope(155:326,1)./max(abs(slope(155:326,1))).*1; 

%slope(1:154,1) = slope(1:154,1)./0.25; 
%slope(155:326,1) = slope(155:326,1)./0.2; 
%slope(slope(:,1)>1,1)=log10(slope(slope(:,1)>1,1))+1;
%slope(slope(:,1)<-1,1)=-log10(-slope(slope(:,1)<-1,1))-1;

for i = 1:size(maty,1)
    ix0 = (manual_matchpos(:,1) ~= 0);
    
    tol11 = zeros(size(slope(:,1)));
    if ix0(i)
        
        u3pos = manual_matchpos(i,:);
        u3pos(u3pos == 0) = '';
        tol11(i) = slope(i,1); % FLUX ACTIVITY COEFFICIENT
        epsilon2 = tol11(i);
        
        rowpos = u3pos;
        colpos = size(m2.A,2) + 1;
        colpos1 = colpos + 1;
        m2.rhs(rowpos) = -epsilon2;
        
        
        % metrow = tol - alpha + beta
        
        m2.A(rowpos,colpos) = 1;
        m2.A(rowpos,colpos1) = -1;
        m2.vtype(colpos) = 'C';
        
        % set si to be positive
        m2.lb(colpos1) = 0;
        m2.ub(colpos1) = 1000;
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos1) = -1*kappa; % minimized
        m2.obj(colpos) = -1*kappa; % minimized               
    end
end

epsilon = 0;
%kappa2 = 5.5E-3;
for jj = 1:length(model.rxns)
    if model.c(jj)==0
        rxnpos = jj;
        %        xi + si >= -eps2
        %     si >= 0
        %     rho(ri + si)
        % constraint 1
        rowpos = size(m2.A,1) + 1;
        colpos = size(m2.A,2) + 1;
        m2.A(rowpos,rxnpos) = 1;
        m2.A(rowpos,colpos) = 1;
        m2.rhs(rowpos) = -epsilon;
        m2.sense(rowpos) = '>';
        % set si to be positive
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos) = -1*kappa2; % minimized
    
        % constraint 2
        %     xi - ri <= eps2
        %     ri >= 0
        % new row and column
        rowpos = size(m2.A,1) + 1;
        colpos = size(m2.A,2) + 1;
        m2.A(rowpos,rxnpos) = 1;
        m2.A(rowpos,colpos) = -1;
        m2.rhs(rowpos) = epsilon;
        m2.sense(rowpos) = '<';
        % set ri to be positive
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos) = -1*kappa2; % minimized
    end
end

m2.vtype = repmat('C',size(m2.A,2),1);
params.outputflag = 0;
solg1 = gurobi(m2,params);
grate_wt =  solg1.x(find(model_rpmi.c));
solverobj = solg1;
dynamicmodel = m2;

% optional. do gene deletion analysis
if genedelflag
    unqgenes = unique(m2.genes);
    for kk = 1:length(unqgenes),
        model2 = deleteModelGenes(m2,unqgenes(kk));
        solg1 = gurobi(model2,params);
        geneko_growthrate(kk,1) = solg1.x(find(model_rpmi.c));
        geneko_growthrate_obj(kk,1) = solg1.objval;
        
        disp(kk)
        
    end
else
    geneko_growthrate = [];
    geneko_growthrate_obj = [];
end

if rxndelflag
    
    for kk = 1:length(model_rpmi.rxns),
        model2 = m2; model2.lb(kk) = 0; model2.ub(kk) = 0;
        solg1 = gurobi(model2,params);
        rxnko_growthrate(kk,1) = solg1.x(find(model_rpmi.c));
        rxnko_growthrate_obj(kk,1) = solg1.objval;
        
        disp(kk)
    end
    
else
    rxnko_growthrate = 0;
    rxnko_growthrate_obj = [];
end

close('all')

end