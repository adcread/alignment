% This script implements interference alignment as per Kamarkar and
% Varanasi (2011).

% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         18/11/14            Initial version


initialisationConstants;                                                     % Initialise the simulation parameters.

SNR = zeros(K);

for i = 1:K
    SNR(i,:) = sqrt(transmitPower(i)) * fading(i,:);
end

if (symmetricChannel == true)
    rho = directChannelFade;
else
    rho = 1.01;                                                                 % Determine a baseline 'rho' to be used as the base for alpha calculations.
end

alpha = log2(SNR)/log2(rho);

