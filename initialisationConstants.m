% Initialise the channel model and pre-set some constants that define it.
% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         18/11/14            Initial version
%   1.0         23/12/14            First Formalised Version

addpath('Adaptive Modulation');
addpath('DoF Calculation');
addpath('Equalisation');

users = 2;                                      % K = number of users in network

txAntennas = [4 2];                             % M = number of transmit antennas
rxAntennas = [4 2];                             % N = number of receive antennas

maxIter = 100;

alpha = [1 .6; .6 1];

beta = zeros(users);

for i=1:users
    for j = 1:users
        beta(i,j) = max([0 (alpha(i,i)-alpha(i,j))]);
    end
end

baselinePower = 5000;
baselineNoise = 1;

H = generateChannel(users, txAntennas, rxAntennas, 'kronecker');            % Generates Kronecker channels with unit channel gain

scaleFactor = zeros(users);

for i=1:users
    for j = 1:users
        scaleFactor(i,j) = sqrt(10^(log10(baselinePower)*alpha(i,j)));      % scales the channel to achieve the desired SNR
        H{i,j} = scaleFactor(i,j) * H{i,j};
    end
end

