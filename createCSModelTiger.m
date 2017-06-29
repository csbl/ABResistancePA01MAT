%Trial Script for creating context specific metabolic models for
%Psuedonomas Aeruginosa given normalized gene expression datas from a full
%.soft file for GEO. The script implements TIGER to form the reconstruction
%and the cobra model is then extracted and saved as a .mat file. 

function [CSCobraModels,errorCount,names2,GeneStates,GeneNames,totalScore] = createCSModelTiger(names,expressor,nsamples,tmodel,sd,varargin)
    p = inputParser;
    p.addParamValue('width',15);
    p.addParamValue('percentile',75);
    p.addParamValue('ThreshNames',tmodel.genes);
    p.parse(varargin{:});
    percentile = p.Results.percentile;
    threshName = p.Results.ThreshNames;
    GeneStates = zeros(length(tmodel.genes),nsamples);
    GeneNames = strings(length(tmodel.genes),nsamples);

    cutoff = percentile;

    set_solver('gurobi');

    errorCount = ones(nsamples,1);
    for j = 1:nsamples
        express = expressor(:,j);
        
        names2 = names;
        errorCount(j) = 1-(double(length(express)) / length(tmodel.genes));

        disp('ErrorCount')
        disp(errorCount)
        
        [states,genes,sol,tiger,totalScore(j)] =  gimme(tmodel,express,cutoff,threshName,sd,'gene_names',names2);

        if(j == 1)
            GeneStates = zeros(length(genes),nsamples);
            GeneNames = strings(length(genes),nsamples);
        end
        newModel = extract_cobra(tiger);
        disp(size(states))
        disp(size(GeneStates))
        CSCobraModels(j) = newModel;
        GeneStates(:,j) = states;
        GeneNames(:,j) = genes;

    end

end
