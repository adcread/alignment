function [ mean ] = meanFrames( input, noFrames )
%MEANFRAMES Summary of this function goes here
%   Detailed explanation goes here

    % sanity check
    totalLength = length(input);
    
    if rem(totalLength,noFrames)
        error('source vector length is not a multiple of number of frames');   
    end
    
    sequenceLength = totalLength/noFrames;
    
    for i = 1:(sequenceLength)
        start = (i-1)*sequenceLength + 1;
        finish = i*sequenceLength - 1;
        

end

