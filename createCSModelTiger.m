%Trial Script for creating context specific metabolic models for
%Psuedonomas Aeruginosa given normalized gene expression datas from a full
%.soft file for GEO. The script implements TIGER to form the reconstruction
%and the cobra model is then extracted and saved as a .mat file. 

function [CSCobraModels,errorCount,names2] = createCSModelTiger(accession,nsamples,tmodel,varargin)
    p = inputParser;
    p.addParamValue('width',15);
    p.parse(varargin{:});

    width = p.Results.width;

    if(startsWith(accession,'GDS') == 1)
        filename = strcat(accession, '_full.soft');

        [names,expressor] = SOFTDataProcessor(filename,nsamples,width);%Process the soft file and fill gaps that result in naming errors
    elseif(startsWith(accession,'GSE')) %Will only work with matrix series files
        [names,expressor] = GSETXTDataProcessor(accession,width);
    else
        disp('Error Encountered, file was not found')
    end
    cutoff = 75;

    set_solver('gurobi');

    minim = min(expressor);
    maxim = max(expressor);
    threshold = prctile(expressor,cutoff);
    disp('Threshold')
    disp(threshold)

    errorCount = ones(nsamples,1);
    for j = 1:nsamples
        express = expressor(:,j);
        for i = 1:length(express)
            if ismember(names(i),tmodel.varnames) == 0
                express(i) = -1*abs(3*minim(j));
            end  
        end

        A = express ~= -1*abs(3*minim(j));
        express = express(A);
        names2 = names(A);

        errorCount(j) = 1-(double(length(express)) / length(tmodel.genes));

        disp('ErrorCount')
        disp(errorCount)

        threshold2(j) = prctile(express,cutoff);


        [states,genes,sol,tiger] =  gimme(tmodel,express,threshold2(j),'gene_names',names2);

        newModel = extract_cobra(tiger);
        
        CSCobraModels(j) = newModel;

    end

end
