SC = load('model_MGSA.mat');
SC = SC.model_MGSA;
path = 'naive_prime_input_data_metabolomics.xlsx';
sheetname = 'naive';
[num txt] = xlsread(path,sheetname);

[dynamicmodel,grate_wt,~, rxnko_growthrate,slope,~,~,~] = flux_activity_coeff(SC,path,sheetname,0,1)

%{
[dynamicmodel,grate_wt,~, rxnko_growthrate,slope,~,~,~] = flux_activity_coeff(fly,path,sheetname,0,1)


[naivemodel, naive_growthrate, naive_genedeletion_growthrate, naive_reactiondeletion_growthrate, ~, ~, naive_genedeletion_objective, naive_reactiondeletion_objective] = flux_activity_coeff(detail, 'naive_prime_input_data_modeling_sep15.xlsx','naive',0, 1);

[primemodel, prime_growthrate, prime_genedeletion_growthrate, prime_reactiondeletion_growthrate, ~, ~, prime_genedeletion_objective, prime_reactiondeletion_objective] = flux_activity_coeff(detail, 'naive_prime_input_data_modeling_sep15.xlsx','prime',0, 1);


[dynamicmodel,grate_wt, geneko_growthrate, rxnko_growthrate,slope,solverobj,geneko_growthrate_obj,rxnko_growthrate_obj] = flux_activity_coeff2(MetModel, 'naive_prime_input_data_metabolomics.xlsx','prime',1,1E-3,0,1)

naivemodel
 
s = struct;
uk = unique(detail.rxns);
for i = 2:length(uk);
    word = string(uk(i));
    ff = strfind(word,',');
    if length(ff)>0
        disp('dot inside')
        disp(word)
        word = split(word,",")
        word = word(end)
    else
        word = word;
    end
    tmp = ['n', word]
    tmp = strjoin(tmp, '')
    s.(tmp) = naive_reactiondeletion_growthrate(i);
end

encodedJSON = jsonencode(s);
fid=fopen('Naive_0_05hr.json','w');
fprintf(fid, encodedJSON);
fclose('all');





uk = unique(detail.rxns);
Rxns = {}
for i = 2:length(uk);
    word = string(uk(i));
    Rxns{i-1, 1} = word;
    values(i-1, 1) = Prime_reactiondeletion_growthrate(i);
end

T = table(Rxns,values);
writetable(T,'Prime_time[1-2hr].csv','Delimiter',';','QuoteStrings',true)

writetable(T,'Naive_constraint.csv','Delimiter',';','QuoteStrings',true)
%}