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

fading = zeros(users);

symmetricChannel = true;
    crossChannelFade = 63.09573445;             % To generate cross-user fading crossChannelFade = 2^(<desired alpha>*log2(directChannelFade))
    directChannelFade = 1000;
    
if (symmetricChannel == true)                   % If the channel is symmetric (i.e. the cross-user links are all equal)
    fading = ones(users) * crossChannelFade;    % then generate the fading channel matrix. Otherwise explicitly define it
    for i = 1:users
        fading(i,i) = directChannelFade;
    end
else      
end

% fading = [1 0.12411 ; 1.12653 1];             % Explicit definition of the channel fade matrix.

transmitPower = ones(1,users);                  % Set all transmitters to equal, unity power output.

H = cell(users,users);

for i = 1:users                                 % Create the empty channel matrix
    for j = 1:users
        H{i,j} = zeros(rxAntennas(i),txAntennas(j));
    end
end                                                

SNR = zeros(users);

for i = 1:users
    SNR(i,:) = sqrt(transmitPower(i)) * fading(i,:);
end

if (symmetricChannel == true)
    rho = directChannelFade;
else
    rho = 2;                                                                % Determine a baseline 'rho' to be used as the base for alpha calculations.
end

alpha = [1 .25; 1.25 1];
beta = zeros(users);

for i=1:users
    for j = 1:users
        beta(i,j) = max([0 (alpha(i,i)-alpha(i,j))]);
    end
end

