function [ noise ] = circSymAWGN( rows, cols, var )
%CIRCSYMAWGN Generates a noise vector of length M with unit variance
%   Detailed explanation goes here

    amp = sqrt(var/2);
    
    noise = amp * (randn([rows,cols]) + 1i *randn([rows,cols]));
    
end

