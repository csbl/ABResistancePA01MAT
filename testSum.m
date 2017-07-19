function [v,meanDeleted,meanAdded] = testSum(allScores,states)

sums = 0;
for i = 1:5
    sums = sums + sum(sum(states{i}./45));
end
sums
meanDeleted = 1146-sums;
meanAdded = mean(allScores);

a = mean(allScores);
n = sums;
if a < (1147-n)/3
    v = n; 
else
    v = n*(1-a/(1147/double(n)));
end

end