function [ average ] = ensembleAve( input )
%ENSEMBLEAVE Summary of this function goes here
%   Detailed explanation goes here

average = zeros(size(input{1}));

for i = 1:length(input)
    average = average + input{i};
end

average = average / length(input);

end

