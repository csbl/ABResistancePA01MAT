clear all

load tmodel.mat;
load Data.mat;

cutoff = 75;

express = cell2mat(Data(:,2));
names = Data(:,1);

set_solver('gurobi');

minim = min(express);
maxim = max(express);
threshold = prctile(express,cutoff);
disp('Threshold')
disp(threshold)

errorCount = 0;
i =1;

for i = 1:length(express)
    if ismember(names(i),tmodel.varnames) == 0% && express(i) > threshold
        %disp(express(i))
        express(i) = -1*abs(3*minim);%[];% minim;
        %names(i) = [];
        errorCount = errorCount + 1;
    end  
    %i = i+1;
end

A = express ~= -1*abs(3*minim);
express = express(A);
names = names(A);

disp('ErrorCount')
disp(errorCount)

threshold2 = prctile(express,cutoff);


[states,genes,sol,tiger] =  gimme(tmodel,express,threshold2,'gene_names',names);

newModel = extract_cobra(tiger)

