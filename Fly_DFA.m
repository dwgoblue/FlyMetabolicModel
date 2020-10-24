fly = load('FlySilico.mat');
fly = fly.FlySilico_v1;
path = 'Mets_G3_PLOS_GR.xlsx';
sheetname_early = 'PLOS_FirstHalf';
sheetname_late = 'PLOS_SecondHalf';
sheetname = 'PLOS_BestMatch';
[num txt] = xlsread(path,sheetname);


% Find genes.
uk = unique(fly.rxns);
uk(2:length(uk)+1, 1) = uk(:);
uk{1, 1} = 'mets';
uk{length(uk)+1, 1} = 'WT';

% Change the lower and upper bounds of oxygen uptake.
l_v = 0:-0.1:-2;
%%%%%%%%%%
%l_v = -2*ones(length(l_v));
%%%%%%%%%%%
V1 = zeros(length(uk), length(l_v));
V2 = zeros(length(uk), length(l_v));
S1 = zeros(length(txt(2:end, 1))+1, length(l_v));
S2 = zeros(length(txt(2:end, 1))+1, length(l_v));

for i=1:length(l_v);
    model = changeRxnBounds(fly,'EX_o2',l_v(i),'l'); % Set the maximum uptake rate.

    % Early stage
    [dynamicmodel,grate_wt,~,rxnko_growthrate,slope,~,~,~] = flux_activity_coeff2(model,path,sheetname_early,1,1E-3,0,1);

    % Find growth rates for gene deletions.
    gh = rxnko_growthrate;
    % Append growth rate for the wild-type strain.
    gh(end+1, 1) = grate_wt;
    % Save result in tb.
    V1(1, i) = l_v(i);
    V1(2:end, i) = gh(:);
    % Save result in tb.
    S1(1, i) = l_v(i);
    S1(2:end, i) = slope(:, 1);

    [dynamicmodel,grate_wt,~,rxnko_growthrate,slope,~,~,~] = flux_activity_coeff2(model,path,sheetname_late,1,1E-3,0,1);
    rxnko_growthrate
    % Late stage
    % Find growth rates for gene deletions.
    gh = rxnko_growthrate;
    % Append growth rate for the wild-type strain.
    gh(end+1, 1) = grate_wt;
    % Save result in tb.
    V2(1, i) = l_v(i);
    V2(2:end, i) = gh(:);
    % Save result in tb.
    S2(1, i) = l_v(i);
    S2(2:end, i) = slope(:, 1);
end

% Initialize tables for saving data from early and late stages.
T1 = table(uk, V1);
T2 = table(uk, V2);
T1S = table(txt(1:end, 1), S1);
T2S = table(txt(1:end, 1), S2);
% Save results in CSV.
writetable(T2,'[hypoxia]PLOS_FlyLateDFA.csv','Delimiter',';','QuoteStrings',true)
writetable(T1,'[hypoxia]PLOS_FlyEarlyDFA.csv','Delimiter',';','QuoteStrings',true)
writetable(T2S,'[hypoxia]PLOS_Fly_late_slope.csv','Delimiter',';','QuoteStrings',true)
writetable(T1S,'[hypoxia]PLOS_Fly_early_slope.csv','Delimiter',';','QuoteStrings',true)