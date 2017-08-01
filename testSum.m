%%Objective function for finding the optimal parameter set the modified
%GIMME Algorithim
% [v,n,a] = testSum(allScores,states)
%Outputs:
%   v = objective value
%   n = number of deleted genes
%   a = mean number of genes added back into models
%Inputs:
%   allScores = vector of number of genes that had to be added back to make
%   functional model
%   states = cell array of matrices of gene states corresponding to the state of each gene
%   in the models. 
function [v,n,a] = testSum(allScores,states)

sums = 0;
for i = 1:5
    sums = sums + sum(sum(states{i}./45));%find mean number of genes in the model
end
sums %display the value 

a = mean(allScores)% find mean number added back into model
n = 1137 - sums % find number of genes deleted
if a < (sums)*2/5. %objective function
    v = n; 
else
    v = n - (a - sums*2/5.);
end
 v %display objective value
end