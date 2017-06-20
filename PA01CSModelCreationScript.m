%Example Script that uses the function createCSModelTiger to create the CS
%models and save the models as .mat cobra models for analysis in python.

clear all
load tmodel.mat;
assession = string(['GDS4244';'GDS3572']);
width = 10;
nsamples = [12,6];
percent = 50;

for i = 1:length(nsamples)
    
    disp(assession(i))

    [CSCobraModels,errorCount,names2] = createCSModelTiger(assession(i),nsamples(i),tmodel,'width',width,'percentile',percent);

    for j = 1:nsamples(i)
        tempModel = CSCobraModels(j);
        fil = sprintf('%s_%d.mat',assession(i),j);
        save(fil,'tempModel')
    end
end

assession = string(['GSE90620';'GSE65870';'GSE30021']);
fileName = strcat(assession,'_series_matrix.txt');
width = 10;
nsamples = [12,6,9];

for i = 1:length(nsamples)
    [CSCobraModels,errorCount,names2] = createCSModelTiger(fileName(i),nsamples(i),tmodel,'width',width,'percentile',percent);

    for j = 1:nsamples(i)
        tempModel = CSCobraModels(j);
        fil = sprintf('%s_%d.mat',assession(i),j);
        save(fil,'tempModel')
    end
end