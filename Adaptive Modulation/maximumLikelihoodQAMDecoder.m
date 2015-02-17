function [ index, estimate ] = maximumLikelihoodQAMDecoder( input, codingIndex )
%MAXIMUMLIKELIHOODQAMDECODER Summary of this function goes here
%   Detailed explanation goes here

% create the set of codewords to perform ML detection on
codebook = generateReceiveCodebook(codingIndex);

numberOfCodewords = 2^codingIndex;

distance = zeros(1,numberOfCodewords);

for i = 1:numberOfCodewords
    distance(i) = abs(input-codebook(i));
end

[confidence, index] = min(distance);

estimate = codebook(index);

end

