function [ names,expressor,med,sd,names2 ] = dataProcessor( accession , nsamples , width, tmodel , normalFact)
%dataProcessor: Processes the gene expression data identified by the
%accession numbers in accession. Outputs the normalized and processed gene expression data with 0 mean and .1 variance. Normalized
%to the mean. Along with the median and sd of each genes expression value.
%Inputs:
%   accession   vector of GEO accession numbers corresponding to saved GE data in .soft or .txt format
%   nsamples    vector of number of samples corresponding to the accession
%               numbers above
%   width       parameter for matching GE data gene names with model gene
%               names, 10 is a good number
%   tmodel      tiger model of which the gene names will be used
%   normalFact  percentile to normalize expression
%Outputs:
%   names       gene names for each sample
%   expressor   procesed gene expression values corresponding to names
%   med         median expression level for all genes in tmodel.genes
%   sd          standard deviation of expression for all genes in
%               tmodel.genes
%   names2      names corresponding to tmodel.genes
    %process data according to type
    for i = 1 : length(accession)
        if(startsWith(accession{i},'GDS') == 1)
            filename = strcat(accession{i}, '_full.soft');
            [names{i},expressor{i}] = SOFTDataProcessor(filename,nsamples(i),width,tmodel,normalFact);%Process the soft file and fill gaps that result in naming errors
        elseif(startsWith(accession{i},'GSE')) %Will only work with matrix series files
            fileName = strcat(accession{i},'_series_matrix.txt');
            [names{i},expressor{i}] = GSETXTDataProcessor(fileName,width,tmodel,normalFact);
        else
        end
    end
    
    for i = 1:length(expressor)
    [numRows(i),numCol(i)] = size(expressor{i}) ;
    end
    %calculate median and sd
    [~,I] = max(numRows);
    names2 = names{I};
    sums = zeros(length(names2),1);
    numbers = zeros(length(names2),1);
    sd = zeros(length(names2),1);
    med = zeros(length(names2),1);
    for i = 1 : length(names2)
       tempSD = [];
       z = 1;
       for j = 1 : length(expressor)
           k = find(strcmp(names{j}(:,1),names2(i))==1,1);
           if ~isempty(k)
              for l = 1:numCol(j)
                  temp = expressor{j}(:,l);
                  tempSD(z) = temp(k);
                  sums(i) = temp(k) + sums(i);
                  numbers(i) = numbers(i) + 1;
                  z= z+1;
              end
           end
       end
       sd(i) = std(tempSD);
       med(i) = median(tempSD);
    end
end