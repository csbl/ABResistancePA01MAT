
function [CSCobraModels,errorCount,names2,GeneStates,GeneNames,totalScore] = createCSModelTiger(names,expressor,nsamples,tmodel,sd,parameter,cutoff,varargin)
%Function for creating CS models from a base tiger model and expression
%data NOTE: AS OF NOW EXTRACT COBRA DOES NOT WORK, SO THE GENE STATES MUST
%BE MANUALLY ENFORCED USING DELETEMODELGENES!!!!!!!!!!!!!
%Inputs: names names of genes 
%   names   names of genes corresponding to the expresssion levels in
%           expressor
%   expressor  processed gene expression data
%   nsamples   number of samples in gene expression data
%   tmodel     base tiger model to derive context-specific models
%   sd         standard deviation of genes corresponding to the genes in
%              tmodel.genes
%   parameter  parameters for use in gimmeESMod algorithim. See
%              documentation
%   cutoff     median expression value for genes in tmodel.genes
%Optional Inputs
%   ThreshNames names for gene corresponding to cutoff, use if different than tmodel.genes
%Outputs
%   CSCobraModels   vector of context-specific cobra models(NOT FUNCTIONAL)
%   errorCount      number of unmatched genes from expression data to
%                   tmodel
%   names2          names of genes in CS model
%   GeneStates      states of genes 1 = on, 0 = off
%   GeneNames       Names of genes corresponding to GeneStates
%   totalScore      number of genes added back into model by gimmeESMod

    p = inputParser;
    p.addParamValue('ThreshNames',tmodel.genes);
    p.parse(varargin{:});
    percentile = p.Results.percentile;
    threshName = p.Results.ThreshNames;
    GeneStates = zeros(length(tmodel.genes),nsamples);
    GeneNames = strings(length(tmodel.genes),nsamples);


    set_solver('gurobi');

    errorCount = ones(nsamples,1);
    for j = 1:nsamples
        express = expressor(:,j);
        
        names2 = names;
        errorCount(j) = 1-(double(length(express)) / length(tmodel.genes));

        disp('ErrorCount')
        disp(errorCount)
        
        [states,genes,sol,tiger,totalScore(j)] =  gimmeESMod(tmodel,express,parameter,cutoff,threshName,sd,'gene_names',names2);

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
