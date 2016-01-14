function [ R ] = crosscorr( vec1, vec2 )
%CROSSCORR Summary of this function goes here
%   Detailed explanation goes here

veclength = max([length(vec1) length(vec2)]);

for i = 0:veclength-1
    shiftVec = circshift(vec2,i);
    R(i+1) = vec1'*shiftVec;
end

