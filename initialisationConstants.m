% Initialise the channel model and pre-set some constants that define it.
% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         18/11/14            Initial version
%   1.0         23/12/14            First Formalised Version

users = 2;                                      % K = number of users in network

txAntennas = [3 2];                             % M = number of transmit antennas
rxAntennas = [3 2];                             % N = number of receive antennas

maxIter = 100;

alpha = [1 .4; .4 1];

beta = zeros(users);

for i=1:users
    for j = 1:users
        beta(i,j) = max([0 (alpha(i,i)-alpha(i,j))]);
    end
end

baselinePower = 2;



fading = 

H = generateChannel(users, txAntennas, rxAntennas,
