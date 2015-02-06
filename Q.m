function [ result ] = Q( x )
%Q Summary of this function goes here
%   Detailed explanation goes here

result = 0.5 * erfc(x/sqrt(2));

end

