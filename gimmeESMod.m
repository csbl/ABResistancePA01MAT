function [gene_states,genes,sol,tiger,error] = ...
                                      gimmeESMod(tiger,express,parameter,thresh,threshNames,sd,varargin)
% GIMME  Gene Inactivity Moderated by Metabolism and Expression
%
%   [GENE_STATES,GENES,SOL,TIGER] = 
%       GIMME(TIGER,EXPRESS,PARAMETER,THRESH,THRESHNAMES,SD,...parameters...)
%
%   Integrate expression data with a metabolic model using a modified
%   GIMME algorithm [Becker & Palsson (2008), PLoS Comput Biol].  This
%   implementation uses an integrated model to map expression values to
%   the gene instead of averaging over each reaction.
%
%   Inputs
%   TIGER    TIGER model.  COBRA models will be converted to TIGER models
%            with a warning.
%   EXPRESS  Expression value for each gene.
%   Parameter Vector of length 5 storing parameter values to be used in the
%       algorithm. [p1,p2,p3,p4,p5]: p1 = percentile of sample expression to be
%       used in cutoffs. (p2-50)/100 = fraction of sd to be added to median value to be
%       used as cutoff in combination with the p1 cutoff. p3/100 = fraction of sd
%       to be subtracted from median to turn off gene. p4 = percentile of
%       sample to be used a cutoff for low expression. p5 = percentile of sd
%       vector to be used with p4 to turn of consistently lowly expressed genes
%   THRESH   Threshold for genes to be turned "on". Corresponding to the
%       genes names given in THRESHNAMES. Normally, is the median expression of
%       individual gene across all samples.
%   THRESHNAMES   Gene names corresponding to the rows of THRESH
%   SD    standard deviation of each genes expression level across all
%       samples.
%
%   Parameters
%   'gene_names'  Cell array of names for genes in EXPRESS.  If not given,
%                 the default is TIGER.genes.
%   'obj_frac'    Fraction of metabolic objective required in the 
%                 resulting model (v_obj >= frac*v_obj_max). 
%                 Default is 0.3.
%
%   Outputs
%   GENE_STATES  Binary expression states calculated by GIMME.
%   GENES        Cell of gene names corresponding to GENE_STATES.
%   SOL          Solution structure with details from the MILP solver.
%   TIGER        TIGER model with GENE_STATES applied.
%   ERROR        Number of genes that had to be turned back on after
%                initial classification

tiger = assert_tiger(tiger);

p = inputParser;
p.addParamValue('gene_names',tiger.genes);
p.addParamValue('obj_frac',0.1);
p.parse(varargin{:});

genes = p.Results.gene_names;
frac = p.Results.obj_frac;

assert(nargin >= 6,'GIMMEESMOD requires at least six inputs');
assert(length(express) == length(genes), ...
       'gene_names must match size of EXPRESS');
   
cut = prctile(express,parameter(1));

j = 1;
w = 0;
for i = 1:length(genes)
   k = find(strcmp(genes(i),threshNames) == 1,1);
   if ~isempty(k) %%%Turn offf gene is expression level is below the median AND is in the bottom half of expression values for the sample. OR expression is less than 1sd lower than the median
       %%ALTERNATIVE FORMAT USING THE MEAN(median) FROM ALL SAMPLES ES 6/23
       if (express(i) < thresh(k) + (parameter(2)-50)/100. *sd(k) ... %{+ (50-parameter)/100*sd(k)
           && express(i) < cut) || express(i) < thresh(k) - sd(k)*parameter(3)/100. %* parameter/100.
       %if (express(i) < thresh(k)+sd(k)/2.)% && express(i) < cut) || express(i) < thresh(k) - sd(k)

          off_locs(i) = 1;
          w(j) = thresh(k) - express(i);
          j = j +1;
       elseif express(i) < prctile(express,parameter(4)) && sd(k) < prctile(sd,parameter(5))
           off_locs(i) = 1;
           w(j) = thresh(k) - express(i);
           j = j +1;
       else
           off_locs(i) = 0;
       end
   else
       off_locs(i) = 0;
   end
  
end

off_locs = logical(off_locs);
totalWeights = zeros(length(threshNames),1);
totalWeights(off_locs) = w; 


off_genes = genes(off_locs); 
[~,off_idxs] = convert_ids(tiger.varnames,off_genes);
if frac > 0
    model = add_growth_constraint(tiger,frac);
else
    model = tiger;
end
model.obj(:) = 0;
model.obj(off_idxs>0) = w; 
sol = cmpi.solve_mip(model);

if cmpi.is_acceptable_exit(sol)
    gene_states = ones(size(express));
    gene_states(off_locs) = round(sol.x(off_idxs));
    tiger = set_var(tiger,genes,gene_states);
    error = 0;
    for i = 1:length(gene_states)
       k = find(strcmp(threshNames,genes(i))==1,1);
       if ~isempty(k) && gene_states(i) == 1 && totalWeights(k) ~= 0
           error = error + 1;%totalWeights(k);
       end
    end
else
    gene_states = [];
end

