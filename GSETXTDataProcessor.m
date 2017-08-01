function [names,express] = GSETXTDataProcessor(fileName,width,tmodel,normalFact)
%Function to process a GSE GEO dataset .txt file given by fileName. Will return normalized and processed gene names and expression values for the sample
%Output:
%   names : cell array of the names of the genes in the file
%   express: matrix of expression data where the columns repressent
%   the different samples and the rows the expression level for the gene
%   corresponding to that given in geneNames
%Input:
%   fileName   name of file .soft file
%   width      width for finding matching gene labels 
%   tmodel     tiger model to which the gene names in the sample should be
%              matched
%   normalFact percentile to normalize the data 

geoData = geoseriesread(char(fileName));%read data

geneNames = string();

for i = 1:length(geoData.Data.RowNames)
   temp = char(geoData.Data.RowNames(i));
   geneNames(i,1) = temp(1,1:6);%gather geneNames
end

condition = 1;
numAddP = 0;
numAdd = 0;
j=0;
while(condition)%match samples names with model names
    for i =2:length(geneNames)
        if startsWith(geneNames(i-1),'PA') == 1 && PADownStream(geneNames,i,width) == 1 && startsWith(geneNames(i),'PA') ==0 
            value = str2num(temp(3:length(temp)))+1;
            value = sprintf('%04d',value);
            geneNames(i,1) = string(['P','A',value]); 
            numAdd = numAdd+1;
            if(j == 0)
                
                
            end
            j = j+1;
        end
    end
    if(numAdd - numAddP == 0)
        condition = 0;
    end
    numAddP = numAdd;
    numAdd = 0;
end


Expressor = double(geoData.Data);%get data
[~,numSamp] = size(geoData.Data);%get number of samples


for j = 1:numSamp
    for i = 1 : length(geneNames)
         if ismember(geneNames(i),tmodel.varnames) == 0
            Expressor(i,j) = -1;
         end 
    end
end

A = Expressor ~= -1;
expresser = Expressor(A);%pull out matching expression values
express = vec2mat(expresser,numSamp);
for j = 1 : numSamp%normalize, set mean to 0, transform variance to .1
    express(:,j) = ((express(:,j) - min(express(:,j))) ./ (prctile(express(:,j),normalFact) - min(express(:,j))));
    express(:,j) = express(:,j) - mean(express(:,j));
    express(:,j) = express(:,j) .* tranformVariance(express(:,j),.1) ;

end
names = geneNames(A(:,1));
figure;
for i = 1 : numSamp
    scatter(1:length(express(:,i)),express(:,i))
end
end


function [result] = PADownStream(geneNames,i,width)

if (i + width) > length(geneNames)
    result = max(startsWith(geneNames(i+1:length(geneNames)),'PA'));
else 
    result = max(startsWith(geneNames(i+1:i+width),'PA'));
end

if result == 1
    result = true;
else
    result = false;
end
end