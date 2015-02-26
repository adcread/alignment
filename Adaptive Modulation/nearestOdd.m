function [ out ] = nearestOdd( in )
%NEARESTODD Summary of this function goes here
%   Detailed explanation goes here

res = mod(in,2);
if res >= 1
    out = 2*(floor(in/2))+1;
else
    out = 2*(floor(in/2))-1;
end

end

