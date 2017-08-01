function [names,express] = SOFTDataProcessor(fileName,numSamp,width,tmodel,normalFact)
%Function to process a .soft file given by fileName with numSamp number of
%samples. Will return normalized and processed gene names and expression values for the sample
%Output:
%   names : cell array of the names of the genes in the file
%   express: matrix of expression data where the columns repressent
%   the different samples and the rows the expression level for the gene
%   corresponding to that given in geneNames
%Input:
%   fileName   name of file .soft file
%   numSamp    number of samples in file
%   width      width for finding matching gene labels 
%   tmodel     tiger model to which the gene names in the sample should be
%              matched
%   normalFact percentile to normalize the data 

SoftData = geosoftread(char(fileName));
geneNames = SoftData.Data(:,2);
ExpressionData = cell2mat(SoftData.Data(:,3:numSamp+2));
geneNames2 =  string();
condition = 1;
numAddP = 0;
numAdd = 0;

% Match sample names with model gene names assuming in a sorted order
while(condition)
    for i =2:length(geneNames)
        if startsWith(geneNames(i-1),'PA') == 1 && max(startsWith(geneNames(i+1:i+width),'PA')) == 1 && startsWith(geneNames(i),'PA') ==0 
            temp = char(geneNames(i-1));
            value = str2num(temp(3:length(temp)))+1;
            value = sprintf('%04d',value);
            geneNames2(i,1) = string(['P','A',value]); 
            numAdd = numAdd+1;
        else
            geneNames2(i,1) = string(geneNames(i));
        end

    end
    geneNames = cellstr(geneNames2);
    if(numAdd - numAddP == 0)
        condition = 0;
    end
    numAddP = numAdd;
    numAdd = 0;
end

for j = 1:numSamp
    for i = 1 : length(geneNames)
         if ismember(geneNames(i),tmodel.varnames) == 0
            ExpressionData(i,j) = -1;
         end 
    end
end

A = ExpressionData ~= -1;
expresser = ExpressionData(A);
express = vec2mat(ExpressionData(A),numSamp);%pull out expression values for matched genes
for j = 1 : numSamp%normalize, set mean to 0, and transform variance to .1
    express(:,j) = ((express(:,j) - min(express(:,j))) ./ (prctile(express(:,j),normalFact) - min(express(:,j))));
    express(:,j) = express(:,j) - mean(express(:,j));
    express(:,j) = express(:,j) .* tranformVariance(express(:,j),.1) ;
end
figure;%plot the sample data
for i = 1 : numSamp
    scatter(1:length(express(:,i)),express(:,i))
end
names = geneNames(A(:,1));
end