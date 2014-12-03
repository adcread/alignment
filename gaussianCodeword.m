function [ X ] = gaussianCodeword( M,P )
%GAUSSCODEWORD Generates independently encoded Gaussian codewords of length M with power P
%   Detailed explanation goes here

X = zeros(M,1);
f = sqrt(P);

for m = 1:M
X(m) = f * (randn(1) + 1i *randn(1));

end

