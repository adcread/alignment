function [ average ] = ensembleAve( ensemble )
%ENSEMBLEAVE Calculate matrix average over ensemble stored as cell
%   Detailed explanation goes here

ensembleSize = max(size(ensemble));

average = zeros(size(ensemble{1}));

for i = 1:ensembleSize
    average = average + ensemble{i};
end

average = average / ensembleSize;

end

