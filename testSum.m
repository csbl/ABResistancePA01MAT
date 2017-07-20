function [v,meanDeleted,meanAdded] = testSum(allScores,states)

sums = 0;
for i = 1:5
    sums = sums + sum(sum(states{i}./45));
end
sums
meanDeleted = 1137-sums;
meanAdded = mean(allScores)

a = mean(allScores);
n = 1137 - sums
if a < (1137-n)/3
    v = n; 
else
    v = n*(1-a/((1137.0-n)*2));
end

end