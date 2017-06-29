function [ names,expressor,med,sd,names2 ] = dataProcessor( accession , nsamples , width, tmodel , normalFact)
%dataProcessor: Processes the gene expression data identified by the
%accession numbers in accession. Pulls the normalized data out. Normalized
%to the mean. 
%   
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