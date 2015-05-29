function [ seq ] = trainingSequence( type, length )
%TRAININGSEQUENCE Generates a training sequence of codewords for MQAM
%constellations
%   Detailed explanation goes here

if (strcmp(type,'edge'))
        
elseif(strcmp(type,'onetwo'))
    seq = zeros(1,length);
    for i = 1:length
        seq(i) = mod(i+1,2) + 1;
    end
else
    error('Error: sequence type not recognised');
end

end

