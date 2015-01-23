function [ X ] = gaussianCodeword( M,P )
%GAUSSCODEWORD Generates independently encoded Gaussian codewords of length M with power P
%   Detailed explanation goes here

X = sqrt(P/2) * (randn([M,1]) + 1i *randn([M,1]));

end


