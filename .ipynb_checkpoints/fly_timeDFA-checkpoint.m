% Load metabolic network model
fly = load('FlySilico.mat');
fly = fly.FlySilico_v1;

% Load essential files as constraints
path = 'FlyData.xlsx';

% Access the data from the G3 paper
sheetname_early = 'G3_FirstHalf';
sheetname_late = 'G3_SecondHalf';
sheetname = 'G3_BestMatch';

% Access the data from the PLOS ONE paper
sheetname_early = 'PLOS_FirstHalf';
sheetname_late = 'PLOS_SecondHalf';
sheetname = 'PLOS_BestMatch';

% Input mode.
sheetname = input('Please type the sheetname you are going to analyze:', 's')
[num txt] = xlsread(path,sheetname);


% Find genes from the model
uk = unique(fly.rxns);
uk{length(uk)+1, 1} = 'WT';


% Change the lower and upper bounds of oxygen uptake.
o2s = -1;%0:-0.1:-1;
% Time steps
t_step = num(1, 4:end);
t_points = unique(t_step);

%%%%%%%%%%%
%V = zeros(length(uk), length(t_step), length(o2));
%S = zeros(length(txt(2:end, 1))+1, length(t_step), length(o2));
V = {};
S = {};
start_p = [4, 10, 17, 23, 29, 36, 43, 50, 57, 64];
end_p = [16, 22, 28, 35, 42, 49, 56, 63, 69, 76];
for j = 1:length(o2s);
    
    % Change lower bound for oxygen uptake reaction.
    model = changeRxnBounds(fly,'EX_o2',0,'l');
    
    % Create a table for saving results.
    tp1 = "string"; tp2 = "double";
    vn = {}; vn{1} = 'rxns'; vn(2:length(t_points)) = sprintfc('hr%d', t_points(2:end));
    tp = [repelem([tp1, tp2], [1 length(t_points)-1])];
    Size = [length(uk) length(t_points)];
    V = table('Size',Size,'VariableTypes',tp,'VariableNames',vn);
    vn = {}; vn{1} = 'mets'; vn(2:length(t_points)) = sprintfc('hr%d', t_points(2:end));
    SizeSlope = [length(txt(2:end, 1)) length(t_points)];
    S = table('Size',SizeSlope,'VariableTypes',tp,'VariableNames',vn);
    
    % Save mets name in the first columns
    V.rxns = uk;
    S.mets = txt(2:end, 1);

    for i=1:length(t_points)-1;
        % DFA with one time step
        cols = [];
        for tmp_i=1:length(t_step);
            if t_step(tmp_i)==t_points(i) | t_step(tmp_i)==t_points(i+1);
                cols = [cols tmp_i+3];
            end
        end
        % DFA with one time step
        [dynamicmodel,grate_wt,~,rxnko_growthrate,slope,~,~,~] = time_flux_analysis(model,path,sheetname,1,1E-3,0,1,cols);

        % Find growth rates for gene deletions.
        gh = rxnko_growthrate;
        % Append growth rate for the wild-type strain.
        gh(end+1, 1) = grate_wt;
        % Save result in tb.
        V.(sprintf('hr%d', t_points(i+1))) = gh;
        % Save result in tb.
        S.(sprintf('hr%d', t_points(i+1))) = slope(:, 1);
    end
    excelname = 'fly_slope.xlsx';
    excelsheet = sprintf('o2%d', j)
    writetable(V,excelname,'FileType','spreadsheet','Sheet', excelsheet);
    excelsheet = sprintf('slope-o2%d', j)
    writetable(S,excelname,'FileType','spreadsheet','Sheet', excelsheet);
end
% % Initialize tables for saving data from early and late stages.
% T1 = table(uk, V1);
% T1S = table(txt(1:end, 1), S1);
% % Save results in CSV.
% writetable(T1,'[hypoxia]PLOS_FlyLateDFA.csv','Delimiter',';','QuoteStrings',true)
% writetable(T1S,'[hypoxia]PLOS_Fly_late_slope.csv','Delimiter',';','QuoteStrings',true)