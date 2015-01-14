% This script implements interference alignment as per Kamarkar and
% Varanasi (2011).

% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         19/11/14            Initial version


initialisationConstants;                                                    % Initialise the simulation parameters.

d1p = 0;                                                                    % User 1 private DoF
d2p = 0;

d1c1 = 1;                                                                   % User 1 common DoF calculated by User 1
d1c2 = NaN;                                                                   % User 1 common DoF calculated by User 2
d2c1 = NaN;
d2c2 = 0;
d1c = 0;
d2c = 0;

% User 1 determines DoF using the absolute maximum of d1c as side
% information

iteration = 0;

while ((d1c1 >= 0) && ((d1c1 ~=d1c2) || (d2c1 ~=d2c2)))
    [d1p, d1c1, d2c1] = degreesOfFreedom2UserMIMOIC_LGP(1, 100, txAntennas, rxAntennas, alpha, beta, d2c2, 1);

    [d2p, d2c2, d1c2] = degreesOfFreedom2UserMIMOIC_LGP(2, 100, txAntennas, rxAntennas, alpha, beta, d1c1, 1);
    iteration = iteration + 1;
    d1c = min([d1c1 d1c2]);
    d2c = min([d2c1 d2c2]);
    d1c1 = d1c;
    d2c2 = d2c;
    sumDoF = [d1p d1c d2p d2c];
    string = ['Iteration ',num2str(iteration),': ',num2str(sumDoF)];
    disp(string);
end


[U, sigma, V] = eigenchannel(H);                                            % Eigenchannel decomposition

covariancePri = cell(1,users);                                              % set up cells for covariance matrices
covariancePub = cell(1,users);

directionPri = cell(users,users);                   
directionPub = cell(users,users);

cardinality = zeros(users);

for i = 1:users
    for j = 1:users
        cardinality(i,j) = min([rxAntennas(j) txAntennas(i)]);              % cardinality = maximum number of signalling dimensions (smallest number of antennas in channel)
    end
end

eigenvalues = cell(users,users);

for i= 1:users
    for j = 1:users
        eigenvalues{i,j} = zeros(1,cardinality(i,j));
        eigenvalues{i,j} = eigs(H{i,j}*H{i,j}.');                           % creates list of eigenvalues for each channel between users
    end
end

for i=1:users
    for j = 1:users
        if (i ~= j)
            [covariancePri{i}, covariancePub{i}] = covarianceHK(SNR(i,j), H{i,j}); % calculates covariance matrices for private and public messages 
        end
        for k1 = 1:txAntennas(i)
            for k2 = 1:txAntennas(i)
                if (k1 == k2)
                    if (k1 <= cardinality(i,j))
                        r = txAntennas(i) * (1 + SNR(i,j) * eigenvalues{i,j}(k1));  % [Dij]_kk = (1/M_i - 1/r_ik) 1 =< k =<m_ij
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

% privateCodeword{1} = qammod(privateMessageStream{1},4,0,'gray');
% privateCodeword{2} = qammod(privateMessageStream{2},4,0,'gray');

publicCodeword{1} = gaussianCodeword(cardinality(1,2),1);
publicCodeword{2} = gaussianCodeword(cardinality(2,1),1);

% publicCodeword{1} = qammod(publicMessageStream{1},4,0,'gray');
% publicCodeword{2} = qammod(publicMessageStream{2},4,0,'gray');

privateMessage{1} = U{1,2} * sqrt(directionPri{1,2}) * privateCodeword{1};
privateMessage{2} = U{2,1} * sqrt(directionPri{2,1}) * privateCodeword{2};

publicMessage{1} = U{1,2}(:,1:cardinality(1,2)) * sqrt(directionPub{1,2}) * publicCodeword{1};
publicMessage{2} = U{2,1}(:,1:cardinality(2,1)) * sqrt(directionPub{1,2}) * publicCodeword{1};

transmittedMessage{1} = privateMessage{1} + publicMessage{1};
transmittedMessage{2} = privateMessage{2} + publicMessage{2};

for i=1:users
    receivedMessage{i} = zeros(rxAntennas(i),1);
end

        
for i = 1:users
    for j = 1:users
        receivedMessage{i} = receivedMessage{i} + ( H{j,i} * transmittedMessage{j} );   % Yi = sum_j Hji * xj
    end
end





