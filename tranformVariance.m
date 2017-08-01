function [ p ] = tranformVariance( data , v)
%transformVariance Transforms the variance of a mean 0 data set (data) to a
%specified value (v)
    ssx = sum(data.^2);
    n = length(data);
    p = sqrt(v*(n-1)/ssx);
end

