function [ p ] = tranformVariance( data , v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    ssx = sum(data.^2);
    n = length(data);
    p = sqrt(v*(n-1)/ssx);
end

