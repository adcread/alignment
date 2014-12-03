
% Creates a complex Gaussian distribution and generates random channel
% coefficients from it
%
% Version 0.1 Chris Waters 18/11/14

% Change History
%
%   Version     Date                Comments
%   0.1         18/11/14            Initial version

function [ H ] = randomChannel( N, M )
%RANDOMCHANNEL Creates a random complex channel matrix of dimension N x M
%   Detailed explanation goes here

range = [-5:0.01:5];
distribution = normpdf(range,0,1);                                          % Creates a normal distribution with mean 0 and 
                                                                            % variance 1
for n = 1:N
    for m = 1:M
        H(n,m) = (1/sqrt(2)) * complex(normrnd(0,1),normrnd(0,1));
    end
end

end

