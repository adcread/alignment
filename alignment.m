% This script implements interference alignment as per Kamarkar and
% Varanasi (2011).

% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         19/11/14            Initial version


initialisationConstants;                                                    % Initialise the simulation parameters.

SNR = zeros(users);

for i = 1:users
    SNR(i,:) = sqrt(transmitPower(i)) * fading(i,:);
end

if (symmetricChannel == true)
    rho = directChannelFade;
else
    rho = 2;                                                                % Determine a baseline 'rho' to be used as the base for alpha calculations.
end

alpha = log2(SNR)/log2(rho);
beta = zeros(users);

for i=1:users
    for j = 1:users
        beta(i,j) = max([0 (alpha(i,i)-alpha(i,j))]);
    end
end

[d1p, d1c, d2p, d2c] = degreesOfFreedom(txAntennas, rxAntennas, alpha, beta);



for i =1:users                                                                 
    for j = 1:users
        H{i,j} = randomChannel(rxAntennas(j),txAntennas(i));
    end
end

[U, sigma, V] = eigenchannel(H);

covariancePri = cell(1,users);
covariancePub = cell(1,users);

directionPri = cell(users,users);
directionPub = cell(users,users);

for i = 1:users
    for j = 1:users
        cardinality(i,j) = min([rxAntennas(j) txAntennas(i)]);
    end
end

eigenvalues = cell(users,users);

for i= 1:users
    for j = 1:users
        eigenvalues{i,j} = zeros(1,cardinality(i,j));
        eigenvalues{i,j} = eigs(H{i,j}*H{i,j}.');
    end
end

for i=1:users
    for j = 1:users
        if (i ~= j)
            [covariancePri{i}, covariancePub{i}] = covarianceHK(SNR(i,j), H{i,j});
        end
        for k1 = 1:txAntennas(i)
            for k2 = 1:txAntennas(i)
                if (k1 == k2)
                    if (k1 <= cardinality(i))
                        r = txAntennas(i) * (1 + SNR(i,j) * eigenvalues{i,j}(k1));
                        directionPri{i,j}(k1,k2) = 1/r;
                        directionPub{i,j}(k1,k2) = (1/txAntennas(i) - directionPri{1,j}(k1,k2));
                    else
                        directionPri{i,j}(k1,k2) = 1 / txAntennas(i);
                    end
                end
            end
        end
    end
end

% Create the transmitted symbols U, W, X for each transmitter

privateMessage = cell(1,users);
publicMessage = cell(1,users);

privateCodeword = cell(1,users);
publicCodeword = cell(1,users);

transmittedMessage = cell(1,users);
receivedMessage = cell(1,users);

privateCodeword{1} = gaussianCodeword(txAntennas(1),1);
privateCodeword{2} = gaussianCodeword(txAntennas(2),1);

publicCodeword{1} = gaussianCodeword(cardinality(1),1);
publicCodeword{2} = gaussianCodeword(cardinality(2),1);

privateMessage{1} = U{1,2} * sqrt(directionPri{1,2}) * privateCodeword{1};
privateMessage{2} = U{2,1} * sqrt(directionPri{2,1}) * privateCodeword{2};

publicMessage{1} = U{1,2}(:,1:cardinality) * sqrt(directionPub{1,2}) * publicCodeword{1};
publicMessage{2} = U{2,1}(:,1:cardinality) * sqrt(directionPub{1,2}) * publicCodeword{1};

transmittedMessage{1} = privateMessage{1} + publicMessage{1};
transmittedMessage{2} = privateMessage{2} + publicMessage{2};

for i=1:users
    receivedMessage{i} = zeros(rxAntennas(i),1);
end

        
for i = 1:users
    for j = 1:users
        receivedMessage{i} = receivedMessage{i} + ( H{j,i} * transmittedMessage{j} );
    end
end




