%% blockalignment.m - performs alignment on a block of data, length determined by the input parameters.


%% Initialisation Phase

addpath('Adaptive Modulation');
addpath('DoF Calculation');
addpath('Equalisation');

dataSymbols = 128;                              % number of data symbols to transmit
trainingSymbols = 32;                           % number of training symbols tranmitted before the data.

totalSymbols = dataSymbols + trainingSymbols;

users = 2;                                      % K = number of users in network

txAntennas = [3 2];                             % M = number of transmit antennas
rxAntennas = [3 2];                             % N = number of receive antennas

power = [1 1];                                  % transmit power levels (per user)

alpha = [1 0.6 ; 0.6 1];

baselinePower = 1000;
baselineNoise = 1;

SNR = (baselinePower/baselineNoise) .^ alpha;   % work out SNR value for given alpha and baseline power levels

%% Creation of Channel matrices

H = generateChannel(users, txAntennas, rxAntennas, 'kronecker');

H{2,1} = H{1,2}';

%% Calculation of Degrees of Freedom in Network

privateDoF = [1 1];
publicDoF = [1 1];

cardinality = zeros(users);                                                 % Cardinality = maximum number of signalling dimensions 
                                                                            % (smallest number of antennas in channel)
for i = 1:users
    for j = 1:users
        cardinality(i,j) = min([rxAntennas(j) txAntennas(i)]);              
    end
end

spareDimensions(1) = pos(txAntennas(1) - cardinality(1,2));
spareDimensions(2) = pos(txAntennas(2) - cardinality(2,1));

for user = 1:users                                                          % Allocate a whole DoF to null spaces of cross-user channel when available
    if spareDimensions(user)>0                                              % (i.e. spareDimension(x)>0)
        for dim =1:spareDimensions(user)
            dofSplitPri{user}(txAntennas(user)-(dim-1)) = 1;
            privateDoF(user) = privateDoF(user) - 1;
        end
    end
end

for stream = 1:cardinality(1,2)
    dofSplitPri{1}(stream) = privateDoF(1)/cardinality(1,2);
    dofSplitPub{1}(stream) = publicDoF(1)/cardinality(1,2);
end

for stream = 1:cardinality(2,1)
    dofSplitPri{2}(stream) = privateDoF(2)/cardinality(2,1);
    dofSplitPub{2}(stream) = publicDoF(2)/cardinality(2,1);
end
                                                                            %%%%%%%%%%%%%%%%%%%%%%
dofSplitPri{1} = [0.4 0.4 1];                                               % PARAMETER TO CHANGE 
dofSplitPub{1} = [0.6 0.6 0];                                               %%%%%%%%%%%%%%%%%%%%%%
dofSplitPri{2} = [0.0 0.0];
dofSplitPub{2} = [0.0 0.0];

%% Creation of Source Alphabets

for user = 1:users
    for stream = 1:txAntennas(user)
       [codebookPri{user}{stream},codebookIndexPri{user}(stream)] = generateTransmitCodebook(dofSplitPri{user}(stream), SNR(user,user),baselineNoise);
       [codebookPub{user}{stream},codebookIndexPub{user}(stream)] = generateTransmitCodebook(dofSplitPub{user}(stream), SNR(user,user),baselineNoise);        
    end
end

%% Creation of messages to transmit

[U, sigma, V] = eigenchannel(H);  % Eigenchannel decomposition

covariancePri = cell(1,users);  % set up cells for covariance matrices
covariancePub = cell(1,users);

% Determine the direction matrices from the covariance matrices of the
% transmitted messages as per Karmakar, eq. 14 & 16 

directionPri = cell(users,users);                   
directionPub = cell(users,users);

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

privateMessage = cell(totalSymbols,users);                                  % Cells are being used for the time dimension so that the signal vectors don't need
publicMessage = cell(totalSymbols,users);                                   % manipulation.

privateSymbol = cell(totalSymbols,users);
publicSymbol = cell(totalSymbols,users);

privateCodeword = cell(totalSymbols,users);
publicCodeword = cell(totalSymbols,users);

transmittedMessage = cell(totalSymbols,users);
receivedMessage = cell(totalSymbols,users);


for symbol = 1:totalSymbols                                                 % Preallocate memory for codewords created
    
    for user = 1:users
        privateSymbol{symbol,user} = zeros(txAntennas(user),1);
        publicSymbol{symbol,user} = zeros(txAntennas(user),1);
    end
end
    
for symbol = 1:trainingSymbols

end

for symbol = (trainingSymbols+1):totalSymbols
    
    for user = 1:users
        for stream = 1:txAntennas(user)
            MPri = length(codebookPri{user}{stream});
            MPub = length(codebookPub{user}{stream});
            if (MPri > 1)
                privateCodeword{symbol,user}(stream) = (randi(MPri));
                privateSymbol{symbol,user}(stream) = codebookPri{user}{stream}(privateCodeword{symbol,user}(stream));
            end
            if (MPub > 1)
                publicCodeword{symbol,user}(stream) = (randi(MPub));
                publicSymbol{user}(stream) = codebookPub{user}{stream}(publicCodeword{symbol,user}(stream));
            end
        end
    end

    privateMessage{symbol,1} = V{1,2} * sqrt(directionPri{1,2}) * privateSymbol{symbol,1};
    privateMessage{symbol,2} = V{2,1} * sqrt(directionPri{2,1}) * privateSymbol{symbol,2};

    publicMessage{symbol,1} = V{1,2} * sqrt(directionPub{1,2}) * publicSymbol{symbol,1};
    publicMessage{symbol,2} = V{2,1} * sqrt(directionPub{2,1}) * publicSymbol{symbol,2};

    transmittedMessage{symbol,1} = privateMessage{symbol,1} + publicMessage{symbol,1};
    transmittedMessage{symbol,2} = privateMessage{symbol,2} + publicMessage{symbol,2};

end

%% Passing through the channel

for symbol = 1:totalSymbols
    for rxUser=1:users
        receivedMessage{symbol,rxUser} = zeros(rxAntennas(rxUser),1);
    end

    for rxUser = 1:users
        for txUser = 1:users
            receivedMessage{symbol,rxUser} = receivedMessage{symbol,rxUser} + (sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * transmittedMessage{symbol,txUser} );
        end
    end
    
    for rxUser = 1:users
        receivedMessage{symbol,rxUser} = receivedMessage{symbol,rxUser} + circSymAWGN(rxAntennas(rxUser),1,1);
    end
end



%% Equalisation

for symbol = 1:trainingSymbols                                              % Perform adaptive equalisation in this section
    
end

for symbol = (trainingSymbols+1):totalSymbols                               % Use calculated Equaliser to perform equalisation in this section

end


%% Detection & Estimation

publicInterference = cell(users,1);
crossInterference = cell(users,1);

decodedPub = cell(users,1);
decodedPri = cell(users,1);
decodedInt = cell(users,1);

for rxUser = 1:users

    publicInterference{rxUser} = zeros(rxAntennas(rxUser),1);
    crossInterference{rxUser} = zeros(rxAntennas(rxUser),1);
    
    commonSubspace = cell(users,1);
    
    for txUser = 1:users
        if (rxUser ~= txUser)

            if (SNR(rxUser,rxUser) >= SNR(txUser,rxUser))

                % if the SNR of the direct stream is greater than the INR
                % from user B, decode the common message from user A first

                commonSubspace{rxUser} = receivedMessage{rxUser};
                
                if (max(dofSplitPub{txUser}) > 0) % if public stream has been sent from the interfering user B
                               
                    if (rxAntennas(rxUser) > txAntennas(txUser))

                        % If user A has more antennas then strip out the extra stream(s) so that
                        % user B's interference can be subtracted
                        
                        numberOfDimensions = rxAntennas(rxUser);
                        hiddenDimensions = rxAntennas(rxUser) - txAntennas(txUser);
                        commonSubspace{rxUser} = removeHiddenStreams(commonSubspace{rxUser},(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
                    end

                    % Equalise user B's interference using SVD of the cross channel
                    
                    equalisedInt{rxUser} = pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * commonSubspace{rxUser};

                    equalisedInt{rxUser} = 1/sqrt(SNR(txUser,rxUser)) * equalisedInt{rxUser};       % Attenuate to bring constellation back to alphabet


                    % decode the public message from the unwanted user B as interference

                    for i = 1:txAntennas(txUser)
                    [dataInt{rxUser}(i), decodedInt{rxUser}(i)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(i),codebookIndexPub{txUser}(i));
                    end

                    % remove the effect of the decoded message from the received signal 

                    crossInterference{rxUser} = H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * decodedInt{rxUser}.';

                end
                
                commonSubspace{rxUser} = commonSubspace{rxUser} - crossInterference{rxUser};
                
                if (max(dofSplitPub{rxUser}) > 0)

                    equalisedPub{rxUser} = pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser};

                    equalisedPub{rxUser} = 1/sqrt(SNR(txUser,rxUser)) * equalisedPub{rxUser};       % Attenuate to bring constellation back to alphabet
                    
                    % decode the message and symbols sent by user A

                    for i = 1:length(equalisedPub{rxUser})
                        [dataPub{rxUser}(i), decodedPub{rxUser}(i)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(i),codebookIndexPub{rxUser}(i));
                    end
                    % remove the effect of the decoded message from the
                    % received signal

                    publicInterference{rxUser} = H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}.';                  
                    commonSubspace{rxUser} = commonSubspace{rxUser} - publicInterference{rxUser};

                end

                if (max(dofSplitPri{rxUser}) > 0)

                    isolatedPri{rxUser} = receivedMessage{rxUser} - publicInterference{rxUser} - crossInterference{rxUser};

                    equalisedPri{rxUser} = pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * isolatedPri{rxUser};

                    equalisedPri{rxUser} = 1/sqrt(SNR(txUser,rxUser)) * equalisedPri{rxUser};       % Attenuate to bring constellation back to alphabet
                    
                    % decode the private message from the desired user A

                    for i = 1:length(equalisedPri{rxUser})
                        [dataPri{rxUser}(i), decodedPri{rxUser}(i)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(i),codebookIndexPri{rxUser}(i));
                    end

                end
            end         % end of SNR/INR testing loop
        end
    end             % end of txUser loop  
    
end                 % end of rxUser loop


%% Plotting Results


