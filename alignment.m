% This script implements interference alignment as per Kamarkar and
% Varanasi (2011).

% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   0.1         19/11/14            Initial version


initialisationConstants;                                                    % Initialise the simulation parameters.

SNR = zeros(users);

% Create a reciprocal channel from the channel H created in
% initialisationConstants script

H{2,1} = H{1,2}';

for i = 1:users
    for j = 1:users
        SNR(i,j) = trace(H{i,j}*H{i,j}') / baselineNoise;
    end
end

PublicDoF = zeros(1,users);
PrivateDoF = zeros(1,users);

d1p = 1; % User 1 private DoF
d2p = 0;

d1c1 = 0.4;       % User 1 common DoF calculated by User 1
d1c2 = NaN;     % User 1 common DoF calculated by User 2
d2c1 = NaN;     % User 2 common DoF calculated by User 1
d2c2 = 0;       % User 2 common DoF calculated by User 2
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

PrivateDoF(1) = d1p;
PublicDoF(1) = d1c;
PrivateDoF(2) = d2p;
PublicDoF(2) = d2c;

% Determines how DoF is split between the public & private symbols to be
% transmitted

dofSplitPri = cell(1,users);
dofSplitPub = cell(1,users);

for user=1:users
    dofSplitPri{user} = zeros(1,txAntennas(user));
    dofSplitPub{user} = zeros(1,txAntennas(user));
end

% Cardinality = maximum number of signalling dimensions (smallest number of antennas in channel)

cardinality = zeros(users);

for i = 1:users
    for j = 1:users
        cardinality(i,j) = min([rxAntennas(j) txAntennas(i)]);              
    end
end


spareDimensions(1) = pos(txAntennas(1) - cardinality(1,2));
spareDimensions(2) = pos(txAntennas(2) - cardinality(2,1));

% Allocate a whole DoF to null spaces of cross-user channel when available
% (i.e. spareDimension(x)>0)

for user = 1:users
    if spareDimensions(user)>0
        for dim =1:spareDimensions(user)
            dofSplitPri{user}(txAntennas(user)-(dim-1)) = 1;
            PrivateDoF(user) = PrivateDoF(user) - 1;
        end
    end
end

% Allocate the rest of the DoF assigned evenly between the remaining
% streams

for stream = 1:cardinality(1,2)
    dofSplitPri{1}(stream) = PrivateDoF(1)/cardinality(1,2);
    dofSplitPub{1}(stream) = PublicDoF(1)/cardinality(1,2);
end

for stream = 1:cardinality(2,1)
    dofSplitPri{2}(stream) = PrivateDoF(2)/cardinality(2,1);
    dofSplitPub{2}(stream) = PublicDoF(2)/cardinality(2,1);
end

dofSplitPri{1} = [0 0 1];
dofSplitPub{1} = [0 0 0];
dofSplitPri{2} = [0.4 0.4];
dofSplitPub{2} = [0.6 0.6];

for user = 1:users
    for stream = 1:txAntennas(user)
       [codebookPri{user}{stream},codebookIndexPri{user}(stream)] = generateTransmitCodebook(dofSplitPri{user}(stream), SNR(user,user),baselineNoise);
       [codebookPub{user}{stream},codebookIndexPub{user}(stream)] = generateTransmitCodebook(dofSplitPub{user}(stream), SNR(user,user),baselineNoise);        
    end
end

[U, sigma, V] = eigenchannel(H);  % Eigenchannel decomposition

covariancePri = cell(1,users);  % set up cells for covariance matrices
covariancePub = cell(1,users);

directionPri = cell(users,users);                   
directionPub = cell(users,users);

eigenvalues = cell(users,users);

% creates list of eigenvalues for each channel between users - not
% technically required now (21/01/15) since covariance matrices are used to
% construct the direction matrices D and D~

for i= 1:users
    for j = 1:users
        eigenvalues{i,j} = zeros(1,cardinality(i,j));
        eigenvalues{i,j} = eigs(H{i,j}*H{i,j}.');                           
    end
end

% Determine the direction matrices from the covariance matrices of the
% transmitted messages as per Gomadam, eq. 

for i=1:users
    for j = 1:users
        if (i ~= j)
            [covariancePri{i}, covariancePub{i}] = covarianceHK(SNR(i,j), H{i,j}); % calculates covariance matrices for private and public messages
            directionPri{i,j} = inv(V{i,j}) * covariancePri{i} * inv(V{i,j}');
            directionPub{i,j} = eye(txAntennas(i))/txAntennas(i) - directionPri{i,j};
        end
    end
end

% Create the transmitted symbols U, W, X for each transmitter

privateMessage = cell(1,users);
publicMessage = cell(1,users);

privateSymbol = cell(1,users);
publicSymbol = cell(1,users);

transmittedMessage = cell(1,users);
receivedMessage = cell(1,users);

% Preallocate memory for codewords created

for i = 1:users
    privateSymbol{i} = zeros(txAntennas(i),1);
    publicSymbol{i} = zeros(txAntennas(i),1);
end

% Create codewords with DoF splitting as dertermined by algorithm above

privateCodeword = zeros(users, max(txAntennas));
publicCodeword = zeros(users, max(txAntennas));

for user = 1:users
    for stream = 1:txAntennas(user)
        MPri = length(codebookPri{user}{stream});
        MPub = length(codebookPub{user}{stream});
        if (MPri > 1)
            privateCodeword(user,stream) = (randi(MPri));
            privateSymbol{user}(stream) = codebookPri{user}{stream}(privateCodeword(user,stream));
        end
        if (MPub > 1)
            publicCodeword(user,stream) = (randi(MPub));
            publicSymbol{user}(stream) = codebookPub{user}{stream}(publicCodeword(user,stream));
        end
    end
end

privateMessage{1} = V{1,2} * sqrt(directionPri{1,2}) * privateSymbol{1};
privateMessage{2} = V{2,1} * sqrt(directionPri{2,1}) * privateSymbol{2};

publicMessage{1} = V{1,2} * sqrt(directionPub{1,2}) * publicSymbol{1};
publicMessage{2} = V{2,1} * sqrt(directionPub{2,1}) * publicSymbol{2};

transmittedMessage{1} = privateMessage{1} + publicMessage{1};
transmittedMessage{2} = privateMessage{2} + publicMessage{2};

for i=1:users
    receivedMessage{i} = zeros(rxAntennas(i),1);
end

% form the received messages as the linear sum of Yi = sum_j(Hij*xj) + w
        
for i = 1:users
    for j = 1:users
        receivedMessage{i} = receivedMessage{i} + ( H{j,i} * transmittedMessage{j} );
    end
end

% for i = 1:users
%     receivedMessage{i} = receivedMessage{i} + circSymAWGN(rxAntennas(i),1,1);
% end

subspace = cell(1,2);
projection = cell(1,2);
orthogonal = cell(1,2);
equaliser = cell(1,2);



% User decoding

% for rxUser = 1:users
%     
%     dataPub{rxUser} = zeros(1,rxAntennas(user));
%     dataPri{rxUser} = zeros(1,rxAntennas(user));
%     
%     % Stores the status of if everything has been decoded from that user.
%     
%     decoded{rxUser} = zeros(1,users);
%         
%     % while not all paths have been decoded
%     while (min(decoded{rxUser})==0)
%         for txUser = 1:users
%             if (rxAntennas(rxUser)>txAntennas(txUser))
%         
%         
            

% User 1 decoding

projectedInt{1} = pseudoInverse(H{1,1}) * receivedMessage{1};

equalisedInt{1} = pseudoInverse(sqrt(directionPub{2,1})) * pseudoInverse(sigma{2,1}) * U{2,1}' * projectedInt{1};



% User 2 decoding

equalisedPub{2} = pseudoInverse(sqrt(directionPub{2,1})) * V{2,1}' * pseudoInverse(H{2,2}) * receivedMessage{2};

for i = 1:length(receivedMessage{2})
    [dataPub{2}(i), decodedPub{2}(i)] = maximumLikelihoodQAMDecoder(equalisedPub{2}(i),codebookIndexPub{2}(i));
end

isolatedPri{2} = receivedMessage{2} - (H{2,2} * V{2,1} * sqrt(directionPub{2,1}) * decodedPub{2}.');

equalisedPri{2} = pseudoInverse(sqrt(directionPri{2,1})) * V{2,1}' * pseudoInverse(H{2,2}) * isolatedPri{2};

for i = 1:length(receivedMessage{2})
    [dataPri{2}(i), decodedPri{2}(i)] = maximumLikelihoodQAMDecoder(equalisedPri{2}(i),codebookIndexPri{2}(i));
end

celldisp(dataPri)
privateCodeword

celldisp(dataPub)
publicCodeword










