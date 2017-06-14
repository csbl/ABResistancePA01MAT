%Example Script that uses the function createCSModelTiger to create the CS
%models and save the models as .mat cobra models for analysis in python.

clear all
load tmodel.mat;
%assession = 'GDS4244';
%width = 10;
%nsamples = 12;

%[CSCobraModels,errorCount,names2] = createCSModelTiger(assession,nsamples,tmodel,'width',width);

%for i = 1:nsamples
    %tempModel = CSCobraModels(i);
    %fil = sprintf('%s_%d.mat',assession,i);
    %save(fil,'tempModel')
%end

assession = 'GSE90620';
fileName = strcat(assession,'_series_matrix.txt');
width = 10;
nsamples = 12;

[CSCobraModels,errorCount,names2] = createCSModelTiger(fileName,nsamples,tmodel,'width',width);

for i = 1:nsamples
    tempModel = CSCobraModels(i);
    fil = sprintf('%s_%d.mat',assession,i);
    save(fil,'tempModel')
end