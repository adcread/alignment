% Initialise the channel model and pre-set some constants that define it.
% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         18/11/14            Initial version

users = 2;                                                                      % K = number of users in network

txAntennas = [3 3];                                                                  % M = number of transmit antennas
rxAntennas = [2 2];                                                                  % N = number of receive antennas

fading = zeros(users);

symmetricChannel = true;
    crossChannelFade = 1.6;
    directChannelFade = 2.0;
    
if (symmetricChannel == true)                                               % If the channel is symmetric (i.e. the cross-user links are all equal)
    fading = ones(users) * crossChannelFade;                                    % then generate the fading channel matrix. Otherwise explicitly define it
    for i = 1:users
        fading(i,i) = directChannelFade;
    end
else                                                                        % Explicit definition of the channel fade matrix.
end

transmitPower = ones(1,users);                                                  % Set all transmitters to equal, unity power output.

H = cell(users,users);

for i = 1:users                                                                 % Create the empty channel matrix
    for j = 1:users
        H{i,j} = zeros(rxAntennas(i),txAntennas(j));
    end
end                                                

