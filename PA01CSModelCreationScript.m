%Example Script that uses the function createCSModelTiger to create the CS
%models and save the models as .mat cobra models for analysis in python.

clear all
load tmodel.mat;
assession = string(['GDS4244';'GDS3572']);
width = 10;
nsamples = [12,6];
percent = 75;
j = 1;
states = cell(5,1);
genes = cell(5,1);
for i = 1:length(nsamples)
    
    disp(assession(i))

    [CSCobraModels,errorCount,names2,states{j},genes{j}] = createCSModelTiger(assession(i),nsamples(i),tmodel,'width',width,'percentile',percent);

    for k = 1:nsamples(i)
        tempModel = CSCobraModels(k);
        fil = sprintf('%s_%d.mat',assession(i),k);
        save(fil,'tempModel')
    end
    
    j= j+1;
end

assession = string(['GSE90620';'GSE65870';'GSE30021']);
fileName = strcat(assession,'_series_matrix.txt');
width = 10;
nsamples = [12,6,9];

for i = 1:length(nsamples)
    [CSCobraModels,errorCount,names2,states{j},genes{j}] = createCSModelTiger(fileName(i),nsamples(i),tmodel,'width',width,'percentile',percent);

    for k = 1:nsamples(i)
        tempModel = CSCobraModels(k);
        fil = sprintf('%s_%d.mat',assession(i),k);
        save(fil,'tempModel')
    end
    j=j+1
end


%NEED TO SAVE THE NEW GENE STATES NO LONGER THE MODEL AND USE THE GENE
%STATES BACK IN PYTHON TO CREATE NEW MODELS WHICH WILL BE USED TO RUN
%SINGLE GENE DELETIONS AGAIN