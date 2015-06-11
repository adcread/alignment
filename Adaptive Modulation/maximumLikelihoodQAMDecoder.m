function [ index, estimate ] = maximumLikelihoodQAMDecoder( input, codebook )
%MAXIMUMLIKELIHOODQAMDECODER Summary of this function goes here
%   Detailed explanation goes here

% create the set of codewords to perform ML detection on

numberOfCodewords = length(codebook);

distance = zeros(1,numberOfCodewords);

for i = 1:numberOfCodewords
    distance(i) = abs(input-codebook(i));
end

[confidence, index] = min(distance);

estimate = codebook(index);

end

