%function [check] = CustomDFA(celltype, mode)

% load metabolomics data.
HumanModel = load('model_MGSA.mat');
% this file has the human metabolicmodel.
MetModel = HumanModel.model_MGSA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load dataset.
path = 'Temporal_Transcriptome.xlsx';
sheetname = 'Processed Transcriptome Data';
[num txt] = xlsread(path,sheetname);

% Save info
gene = txt(2:end, 2);
timepoints = txt(2, 3:end);
bound = size(num);

geneup = {}
genedw = {}

for j=1:bound(2)-1; %col
    u = 1;
    d = 1;
    for i=1:bound(1); %row
        geneup(1, j) = timepoints(j+1);
        genedw(1, j) = timepoints(j+1);
        if num(i, j+1)>0;
            u = u+1;
            geneup(u, j) = gene(i+1);
        elseif num(i, j+1)<0;
            d = d+1;
            genedw(d, j) = gene(i+1);
        else;
            disp('Neutral');
        end
    end
end

for ind=1:length(timepoints)-1;
    up = geneup(:, ind);
    dw = genedw(:, ind);
    
    model = MetModel.genes;
    uk = unique(model);
    uplist = {}

    count=0
    for i=2:length(uk);
        for j=1:length(dw);
            if isempty(dw{j})==0 && upper(string(dw(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                uplist{count, 1} = uk{i};
            end
        end
    end
    uplist = uplist(~cellfun('isempty',uplist));

    uk = unique(model);
    dwlist = {}

    count=0
    for i=1:length(uk);
        for j=1:length(up);
            if isempty(up{j})==0 && upper(string(up(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                dwlist{count, 1} = uk{i};
            end
        end
    end
    dwlist = dwlist(~cellfun('isempty',dwlist));
    
    [fluxstate,grate_naive,solverobj_naive]=constrain_flux_regulation(MetModel, uplist, dwlist, 1, 1, 1E-3, 0);

    T = table(MetModel.rxns,fluxstate);
    time = timepoints(ind)
    filename=sprintf('Naive_transcriptome_%s.csv',string(time))
    writetable(T,filename,'Delimiter',';','QuoteStrings',true)
    
end



for ind=1:length(timepoints)-1;
    dw = geneup(:, ind);
    up = genedw(:, ind);
    
    model = MetModel.genes;
    uk = unique(model);
    uplist = {}

    count=0
    for i=2:length(uk);
        for j=1:length(dw);
            if isempty(dw{j})==0 && upper(string(dw(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                uplist{count, 1} = uk{i};
            end
        end
    end
    uplist = uplist(~cellfun('isempty',uplist));

    uk = unique(model);
    dwlist = {}

    count=0
    for i=1:length(uk);
        for j=1:length(up);
            if isempty(up{j})==0 && upper(string(up(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                dwlist{count, 1} = uk{i};
            end
        end
    end
    dwlist = dwlist(~cellfun('isempty',dwlist));
    
    [fluxstate,grate_naive,solverobj_naive]=constrain_flux_regulation(MetModel, uplist, dwlist, 1, 1, 1E-3, 0);

    T = table(MetModel.rxns,fluxstate);
    time = timepoints(ind)
    filename=sprintf('Primed_transcriptome_%s.csv',string(time))
    writetable(T,filename,'Delimiter',';','QuoteStrings',true)
    
end


% load dataset.
path = 'Temporal_Proteomic.xlsx';
sheetname = 'Processed Proteome Data';
[num txt] = xlsread(path,sheetname);


% Save info
gene = txt(2:end, 2);
timepoints = txt(2, 3:end);
bound = size(num);

geneup = {}
genedw = {}

for j=1:bound(2)-1; %col
    u = 1;
    d = 1;
    for i=1:bound(1); %row
        geneup(1, j) = timepoints(j+1);
        genedw(1, j) = timepoints(j+1);
        if num(i, j+1)>0;
            u = u+1;
            geneup(u, j) = gene(i+1);
        elseif num(i, j+1)<0;
            d = d+1;
            genedw(d, j) = gene(i+1);
        else;
            disp('Neutral');
        end
    end
end

        
    
for ind=1:length(timepoints)-1;
    up = geneup(:, ind);
    dw = genedw(:, ind);
    
    model = MetModel.genes;
    uk = unique(model);
    uplist = {}

    count=0
    for i=2:length(uk);
        for j=1:length(dw);
            if isempty(dw{j})==0 && upper(string(dw(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                uplist{count, 1} = uk{i};
            end
        end
    end
    uplist = uplist(~cellfun('isempty',uplist));

    uk = unique(model);
    dwlist = {}

    count=0
    for i=1:length(uk);
        for j=1:length(up);
            if isempty(up{j})==0 && upper(string(up(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                dwlist{count, 1} = uk{i};
            end
        end
    end
    dwlist = dwlist(~cellfun('isempty',dwlist));
    
    [fluxstate,grate_naive,solverobj_naive]=constrain_flux_regulation(MetModel, uplist, dwlist, 1, 1, 1E-3, 0);

    T = table(MetModel.rxns,fluxstate);
    time = timepoints(ind)
    filename=sprintf('Naive_proteome_%s.csv',string(time))
    writetable(T,filename,'Delimiter',';','QuoteStrings',true)
    
end



for ind=1:length(timepoints)-1;
    dw = geneup(:, ind);
    up = genedw(:, ind);
    
    model = MetModel.genes;
    uk = unique(model);
    uplist = {}

    count=0
    for i=2:length(uk);
        for j=1:length(dw);
            if isempty(dw{j})==0 && upper(string(dw(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                uplist{count, 1} = uk{i};
            end
        end
    end
    uplist = uplist(~cellfun('isempty',uplist));

    uk = unique(model);
    dwlist = {}

    count=0
    for i=1:length(uk);
        for j=1:length(up);
            if isempty(up{j})==0 && upper(string(up(j)))==upper(string(uk(i)));
                count=count+1;
                %disp(count)
                %disp(string(uk(i)))
                dwlist{count, 1} = uk{i};
            end
        end
    end
    dwlist = dwlist(~cellfun('isempty',dwlist));
    
    [fluxstate,grate_naive,solverobj_naive]=constrain_flux_regulation(MetModel, uplist, dwlist, 1, 1, 1E-3, 0);

    T = table(MetModel.rxns,fluxstate);
    time = timepoints(ind)
    filename=sprintf('Primed_proteome_%s.csv',string(time))
    writetable(T,filename,'Delimiter',';','QuoteStrings',true)
    
end
