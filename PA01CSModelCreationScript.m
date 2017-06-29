%Example Script that uses the function createCSModelTiger to create the CS
%models and save the models as .mat cobra models for analysis in python.
%%
clear all
load tmodel.mat;

assession = cellstr(['GDS4244 ';'GDS3572 ';'GSE90620';'GSE65870';'GSE30021']);
width = 10;
nsamples = [12,6,12,6,9];
normalFact = 75;

[names,expressor,ave,sd,names2] = dataProcessor(assession,nsamples,width,tmodel,normalFact);

%%

j = 1;
states = cell(5,1);
genes = cell(5,1); 
allScores = [];
    
for i = 1:length(nsamples)
    
    disp(assession(i))

    [CSCobraModels,errorCount,names3,states{j},genes{j},totalScore] = createCSModelTiger(names{i},expressor{i},nsamples(i),tmodel,sd,'width',width,'percentile',ave,'ThreshNames',names2);

    allScores = [allScores , totalScore];
    for k = 1:nsamples(i)
        tempModel = CSCobraModels(k);
        fil = sprintf('%s_%d',assession{i},k);
        %save(strcat(fil,'.mat'),'tempModel')
        file = fopen(strcat(fil,'.txt'),'w');
        for l = 1:length(states{j})
            fprintf(file,'%s %d\n',genes{j}{l,k},round(states{j}(l,k)));
        end
        %fprintf(file,'%s  %d\n',genes{j}{:,k},round(states{j}(:,k)));
        fclose(file);
    end
    
    j= j+1;
end

%%
controls = [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1];
numE = sum(controls);
numC = length(controls) - numE;
k = 1;
genieE = [];
for j = 1:length(genes)
    A = states{j} == 0;
    [numrows , numcols] = size(states{j});
    for i = 1 : numcols
       if(controls(k) == 1)
           genieE = [genieE ; genes{j}(A(:,i))]; 
       end
       k = k +1;
    end
end

hitsE = zeros(length(tmodel.genes),1);
for i = 1 : length(tmodel.genes)
    hitsE(i) = sum(sum(genieE == tmodel.genes(i)));
end
geneDeletionsE = tmodel.genes(hitsE > 0);
hitsE = hitsE(hitsE > 0);


k = 1;
genieC = [];
for j = 1:length(genes)
    A = states{j} == 0;
    [numrows , numcols] = size(states{j});
    for i = 1 : numcols
       if(controls(k) == 0)
           genieC = [genieC ; genes{j}(A(:,i))]; 
       end
       k = k +1;
    end
end

hitsC = zeros(length(tmodel.genes),1);
for i = 1 : length(tmodel.genes)
    hitsC(i) = sum(sum(genieC == tmodel.genes(i)));
end
geneDeletionsC = tmodel.genes(hitsC > 0);
hitsC = hitsC(hitsC > 0);

summary{1} = hitsC;
summary{2} = geneDeletionsC;
summary{3} = hitsE;
summary{4} = geneDeletionsE;

file = [fopen('deletedGenesC.csv','w'),fopen('deletedGenesE.csv','w')];
minim = min([length(geneDeletionsC),length(geneDeletionsE)]);
[maxim,in] = max([length(geneDeletionsC),length(geneDeletionsE)]);
diff = maxim - minim;
for i = 1:minim
    for j = 1:length(file)
        fprintf(file(j),'%s,%d\n',summary{(j-1)*2+2}{i},summary{(j-1)*2+1}(i));
    end
end

for i = minim+1:maxim
    fprintf(file(in),'%s,%d\n',summary{2+2*(in-1)}{i},summary{1+2*(in-1)}(i));
end



%NEED TO SAVE THE NEW GENE STATES NO LONGER THE MODEL AND USE THE GENE
%STATES BACK IN PYTHON TO CREATE NEW MODELS WHICH WILL BE USED TO RUN
%SINGLE GENE DELETIONS AGAIN