%% blockalignment.m - performs alignment on a block of data, length determined by the input parameters.


%% Initialisation Phase

addpath('Adaptive Modulation');
addpath('DoF Calculation');
addpath('Equalisation');
addpath('General functions');

dataSymbols = 2048;                              % number of data symbols to transmit
trainingSymbols = 0;                           % number of training symbols tranmitted before the data.

totalSymbols = dataSymbols + trainingSymbols;

users = 2;                                      % K = number of users in network

txAntennas = [3 2];                             % M = number of transmit antennas
rxAntennas = [3 2];                             % N = number of receive antennas

power = [1 1];                                  % transmit power levels (per user)

alpha = [1 0.6 ; 0.6 1];

baselinePower = 50000;
baselineNoise = 1;

SNR = (baselinePower/baselineNoise) .^ alpha;   % work out SNR value for given alpha and baseline power levels

DFEon = false;
noiseOn = false;

%% Creation of Channel matrices

H = generateChannel(users, txAntennas, rxAntennas, 'kronecker');

H{2,1} = H{1,2}';                                                           % Channel responses are reciprocal

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
dofSplitPub{1} = [0.2 0.2 0.0];                                               % PARAMETER TO CHANGE 
dofSplitPri{1} = [0.2 0.2 1.0];                                               %%%%%%%%%%%%%%%%%%%%%%
dofSplitPub{2} = [0.4 0.4];
dofSplitPri{2} = [0.4 0.4];


%% Creation of Source Alphabets

codebookPri = cell(users,1);
codebookPub = cell(users,1);

for user = 1:users
    codebookPri{user} = cell(txAntennas(user),1);
    codebookPub{user} = cell(txAntennas(user),1);
end

for user = 1:users
    for stream = 1:txAntennas(user)
        display(['Generating Private codebook for User ' num2str(user) ' stream ' num2str(stream) ':']);
       [codebookPri{user}{stream},codebookIndexPri{user}(stream)] = generateTransmitCodebook(dofSplitPri{user}(stream), SNR(user,user),6);
       display(['Generating Public codebook for User ' num2str(user) ' stream ' num2str(stream) ':']);
       [codebookPub{user}{stream},codebookIndexPub{user}(stream)] = generateTransmitCodebook(dofSplitPub{user}(stream), SNR(user,user),6);      
       display('');
    end
end

%% Creation of messages to transmit

[U, sigma, V] = eigenchannel(H);  % Eigenchannel decomposition

for rxUser = 1:users
    for txUser = 1:users
        condition(rxUser,txUser) = max(nonzeros(max(sigma{rxUser,txUser})))/min(nonzeros(max(sigma{rxUser,txUser})));
    end
end

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
            directionPri{i,j} = pinv(V{i,j}) * covariancePri{i} * pinv(V{i,j}');
            directionPub{i,j} = eye(txAntennas(i))/txAntennas(i) - directionPri{i,j};
            directionPri{i,j} = real(diag(diag(directionPri{i,j})));
            directionPub{i,j} = real(diag(diag(directionPub{i,j})));
        end
    end
end


%%
% Create the transmitted symbols U, W, X for each transmitter

privateCodeword = cell(users,1);                                % Codeword - the element transmitted that carries information.
publicCodeword = cell(users,1);

privateSymbol = cell(users,1);                                  % Symbol - the alphabet element corresponding to a codeword (complex).
publicSymbol = cell(users,1);

privateMessage = cell(users,1);                                 % Message - precoded and altered combination of symbols 
publicMessage = cell(users,1);

transmittedMessage = cell(users,1);                             % Messages at either end of the air interface
receivedMessage = cell(users,1);

equalisedPri = cell(users,1);                                   % Equalised message streams used for decoding
equalisedPub = cell(users,1);
equalisedInt = cell(users,1);

publicInterference = cell(users,1);                             % Equalised interference streams used in SIC
crossInterference = cell(users,1);

dataPub = cell(users,1);                                        % Decoded data and interference streams
dataPri = cell(users,1);
dataInt = cell(users,1);

decodedPub = cell(users,1);                                     % Symbols corresponding to decoded data (used for plotting)
decodedPri = cell(users,1);
decodedInt = cell(users,1);


for user = 1:users
    privateSymbol{user} = zeros(txAntennas(user),totalSymbols);             % dimension of cell: streams * symbols
    publicSymbol{user} = zeros(txAntennas(user),totalSymbols);
end
    
for user = 1:users                                                          % Generate training sequence for equaliser training
    for stream=1:txAntennas(user)
        MPri = length(codebookPri{user}{stream});
        MPub = length(codebookPub{user}{stream});
        publicCodeword{user}(stream,1:trainingSymbols) = ones(1,trainingSymbols);                       
        privateCodeword{user}(stream,1:trainingSymbols) = ones(1,trainingSymbols);                       
        publicCodeword{user}(stream,trainingSymbols+1:totalSymbols) = randi(MPub,1,totalSymbols-trainingSymbols);
        privateCodeword{user}(stream,trainingSymbols+1:totalSymbols) = randi(MPri,1,totalSymbols-trainingSymbols);
    end

end

if (trainingSymbols ~=0)
    publicCodeword{1}(1,1:trainingSymbols) = trainingSequence('onetwo',trainingSymbols);
    publicCodeword{1}(2,1:trainingSymbols) = trainingSequence('onetwo',trainingSymbols);
    publicCodeword{2}(1,1:trainingSymbols) = trainingSequence('onetwo',trainingSymbols);
    publicCodeword{2}(2,1:trainingSymbols) = trainingSequence('onetwo',trainingSymbols);
end

for user = 1:users
    for stream = 1:txAntennas(user)
        privateSymbol{user}(stream,:) = codebookPri{user}{stream}(privateCodeword{user}(stream,:));
        publicSymbol{user}(stream,:) = codebookPub{user}{stream}(publicCodeword{user}(stream,:));   
    end
end


for symbol = (trainingSymbols+1):totalSymbols
    
    for user = 1:users
        for stream = 1:txAntennas(user)
            MPri = length(codebookPri{user}{stream});
            MPub = length(codebookPub{user}{stream});
            if (MPri > 1)
                privateCodeword{user}(stream,symbol) = (randi(MPri));
                privateSymbol{user}(stream,symbol) = codebookPri{user}{stream}(privateCodeword{user}(stream,symbol));
            end
            if (MPub > 1)
                publicCodeword{user}(stream,symbol) = (randi(MPub));
                publicSymbol{user}(stream,symbol) = codebookPub{user}{stream}(publicCodeword{user}(stream,symbol));
            end
        end
    end

end
   
    privateMessage{1} = V{1,2} * sqrt(directionPri{1,2}) * privateSymbol{1};
    privateMessage{2} = V{2,1} * sqrt(directionPri{2,1}) * privateSymbol{2};

    publicMessage{1} = V{1,2} * sqrt(directionPub{1,2}) * publicSymbol{1};
    publicMessage{2} = V{2,1} * sqrt(directionPub{2,1}) * publicSymbol{2};

    transmittedMessage{1} = privateMessage{1} + publicMessage{1};
    transmittedMessage{2} = privateMessage{2} + publicMessage{2};
    
%% Calculate transmitter input power and normalise to <=1

for user = 1:users
    for symbol = 2:length(transmittedMessage{user})
        Q{symbol} = transmittedMessage{user}(:,symbol) * transmittedMessage{user}(:,symbol)';
        G{symbol} = zeros(length(Q{symbol}));
        for i = 1:length(Q{symbol})
            for j = 1:length(Q{symbol})
                G{symbol}(i,j) = G{symbol-1}(i,j) + Q{symbol}(i,j);
            end
        end
        T(symbol) = trace(G{symbol});
        T(symbol) = T(symbol)/symbol;
    end
    
    scaleFactor(user) = sqrt(T(2047));
    
    transmittedMessage{user} = transmittedMessage{user} / scaleFactor(user);
    
    end
    
%% Passing through the channel


for rxUser = 1:users
    receivedMessage{rxUser} = zeros(rxAntennas(rxUser),totalSymbols);
    for txUser = 1:users
        receivedMessage{rxUser} = receivedMessage{rxUser} + (sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * transmittedMessage{txUser});
    end
end

if noiseOn == true
    for rxUser = 1:users
        receivedMessage{rxUser} = receivedMessage{rxUser} + circSymAWGN(rxAntennas(rxUser),totalSymbols,1); % add AWGN to received signal
    end
end

for rxUser = 1:users
    receivedMessage{rxUser} = receivedMessage{rxUser} * scaleFactor(rxUser);
end

%% Equalisation & Detection

for symbol = 1:trainingSymbols
    for rxUser = 1:users

        commonSubspace = cell(users,1);

        for txUser = 1:users
            if (rxUser ~= txUser)

                if (SNR(rxUser,rxUser) >= SNR(txUser,rxUser));

                    % if the SNR of the direct stream is greater than the INR
                    % from user B, decode the common message from user A first

                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol);

                    if (min(dofSplitPub{txUser}) > 0) % if public stream has been sent from the interfering user B

                        if (rxAntennas(rxUser) > txAntennas(txUser))

                            % If user A has more antennas then strip out the extra stream(s) so that
                            % user B's interference can be subtracted

                            numberOfDimensions = rxAntennas(rxUser);
                            hiddenDimensions = rxAntennas(rxUser) - txAntennas(txUser);
                            commonSubspace{rxUser}(:,symbol) = removeHiddenStreams(commonSubspace{rxUser}(:,symbol),(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
                        end

                        % Equalise user B's interference using SVD of the cross channel

                        equalisedInt{rxUser}(:,symbol) = pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * commonSubspace{rxUser}(:,symbol);

%                         equalisedInt{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedInt{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet
                            
                        if (DFEon)
                            for stream = 1:txAntennas(txUser)
                                equalisedInt{rxUser}(stream,symbol) = equalisedInt{rxUser}(stream,symbol) * decisionFeedbackInt{user}(stream,symbol);
                                argErrInt{rxUser}(stream,symbol) = argd(publicSymbol{txUser}(stream,symbol)) - argd(equalisedInt{rxUser}(stream,symbol));
                                magErrInt{rxUser}(stream,symbol) = abs(publicSymbol{txUser}(stream,symbol)) / abs(equalisedInt{rxUser}(stream,symbol));
                                if magErrInt{rxUser}(stream,symbol) < 1
                                    stepInt{rxUser}(stream,symbol) = decisionFeedbackInt{user}(stream,symbol) * stepsize * -1;
                                else
                                    stepInt{rxUser}(stream,symbol) = decisionFeedbackInt{user}(stream,symbol) * stepsize;
                                end
                                decisionFeedbackInt{rxUser}(stream,symbol+1) = decisionFeedbackInt{rxUser}(stream,symbol) + stepInt{rxUser}(stream,symbol);
%                                 decisionFeedbackInt{rxUser}(stream,symbol+1) = decisionFeedbackInt{rxUser}(stream,symbol+1) / p2c(argErrInt{rxUser}(stream,symbol),1);
                            end
                        end
                        
                        % decode the public message from the unwanted user
                        % B as interference - not strictly necessary since
                        % data is not used for anything except training DFE
                        % above.

                        for stream = 1:txAntennas(txUser)
                        [dataInt{rxUser}(stream,symbol), decodedInt{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(stream),codebookPub{txUser}{stream});
                        end

                        % remove the effect of the decoded message from the received signal 

                        crossInterference{rxUser}(:,symbol) = H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * publicSymbol{rxUser}(1:txAntennas(txUser),symbol);
                    else
                        crossInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end

                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                    if (max(dofSplitPub{rxUser}) > 0) % If public messages sent by user A

                        equalisedPub{rxUser}(:,symbol) = pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);

%                         equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPub{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the message and symbols sent by user A
                        
                        if (DFEon)
                            for stream = 1:rxAntennas(rxUser)
                                equalisedPub{rxUser}(stream,symbol) = equalisedPub{rxUser}(stream,symbol) * decisionFeedbackPub{rxUser}(stream,symbol);    
                                argErrPub{rxUser}(stream,symbol) = argd(publicSymbol{rxUser}(stream,symbol)) - argd(equalisedPub{rxUser}(stream,symbol));
                                magErrPub{rxUser}(stream,symbol) = abs(publicSymbol{rxUser}(stream,symbol)) / abs(equalisedPub{rxUser}(stream,symbol));
                                if magErrPub{rxUser}(stream,symbol) < 1
                                    stepPub{rxUser}(stream,symbol) = decisionFeedbackPub{rxUser}(stream,symbol) * stepsize * -1;
                                else
                                    stepPub{rxUser}(stream,symbol) = decisionFeedbackPub{rxUser}(stream,symbol) * stepsize;
                                end
                                decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol) + stepPub{rxUser}(stream,symbol);
%                                 decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol+1) / p2c(argErrPub{rxUser}(stream,symbol),1);
                            end
                        end
                        
                        for stream = 1:length(equalisedPub{rxUser}(:,symbol))
                            [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});
                        end
                        % remove the effect of the decoded message from the
                        % received signal

                        publicInterference{rxUser}(:,symbol) = H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}(:,symbol);                  
                    else
                        publicInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end

                    if (max(dofSplitPri{rxUser}) > 0)   % If private streams from user A received

                        isolatedPri{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - publicInterference{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                        equalisedPri{rxUser}(:,symbol) = pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * isolatedPri{rxUser}(:,symbol);

%                         equalisedPri{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPri{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the private message from the desired user A
                        
                        if (DFEon)
                            for stream = 1:rxAntennas(rxUser)
                                equalisedPri{rxUser}(stream,symbol) = equalisedPri{rxUser}(stream,symbol) * decisionFeedbackPri{rxUser}(stream,symbol);    
                                argErrPri{rxUser}(stream,symbol) = argd(privateSymbol{rxUser}(stream,symbol)) - argd(equalisedPri{rxUser}(stream,symbol));
                                magErrPri{rxUser}(stream,symbol) = abs(privateSymbol{rxUser}(stream,symbol)) / abs(equalisedPri{rxUser}(stream,symbol));
                                if magErrPri{rxUser}(stream,symbol) < 1
                                    stepPri{rxUser}(stream,symbol) = decisionFeedbackPri{rxUser}(stream,symbol) * stepsize * -1;
                                else
                                    stepPri{rxUser}(stream,symbol) = decisionFeedbackPri{rxUser}(stream,symbol) * stepsize;
                                end
                                decisionFeedbackPri{rxUser}(stream,symbol+1) = decisionFeedbackPri{rxUser}(stream,symbol) + stepPri{rxUser}(stream,symbol);
%                                 decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol+1) / p2c(argErrPub{rxUser}(stream,symbol),1);
                            end
                        end
                        
                        for i = 1:length(equalisedPri{rxUser}(:,symbol))
                            [dataPri{rxUser}(i,symbol), decodedPri{rxUser}(i,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(i,symbol),codebookPri{rxUser}{i});
                        end

                    end
                end         % end of SNR/INR testing loop
            end
        end             % end of txUser loop  

    end                 % end of rxUser loop
end

for symbol = (trainingSymbols+1):totalSymbols

    for rxUser = 1:users

        commonSubspace = cell(users,1);

        for txUser = 1:users
            if (rxUser ~= txUser)

                if (SNR(rxUser,rxUser) >= SNR(txUser,rxUser));

                    % if the SNR of the direct stream is greater than the INR
                    % from user B, decode the common message from user A first

                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol);

                    if (min(dofSplitPub{txUser}) > 0) % if public stream has been sent from the interfering user B
                        display(['Decoding Interference from User ' num2str(txUser) '.']);
                        if (rxAntennas(rxUser) > txAntennas(txUser))

                            % If user A has more antennas then strip out the extra stream(s) so that
                            % user B's interference can be subtracted

                            numberOfDimensions = rxAntennas(rxUser);
                            hiddenDimensions = rxAntennas(rxUser) - txAntennas(txUser);
                            commonSubspace{rxUser}(:,symbol) = removeHiddenStreams(commonSubspace{rxUser}(:,symbol),(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
                        end

                        % Equalise user B's interference using SVD of the cross channel

                        equalisedInt{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) *  pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * commonSubspace{rxUser}(:,symbol);

                        [dataInt{rxUser}(stream,symbol), decodedInt{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(stream),codebookPub{txUser}{stream});
                     
                        % remove the effect of the decoded message from the received signal 

                        crossInterference{rxUser}(:,symbol) = H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * decodedInt{rxUser}(:,symbol);
                    else
                        crossInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end

                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                    % If public messages sent by user A
                    
                    if (max(dofSplitPub{rxUser}) > 0) 
                        display(['Decoding Public stream from User ' num2str(rxUser) '.']);
                        equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(rxUser,rxUser)) * pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);

%                         equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPub{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the message and symbols sent by user A
                        
                        if (DFEon)
                            for stream = 1:rxAntennas(rxUser)
                                equalisedPub{rxUser}(stream,symbol) = equalisedPub{rxUser}(stream,symbol) * decisionFeedbackPub{rxUser}(stream,symbol);
                                
                                [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});

                                argErrPub{rxUser}(stream,symbol) = argd(decodedPub{rxUser}(stream,symbol)) - argd(equalisedPub{rxUser}(stream,symbol));
                                magErrPub{rxUser}(stream,symbol) = abs(decodedPub{rxUser}(stream,symbol)) / abs(equalisedPub{rxUser}(stream,symbol));
                                if magErrPub{rxUser}(stream,symbol) < 1
                                    stepPub{rxUser}(stream,symbol) = decisionFeedbackPub{rxUser}(stream,symbol) * stepsize * -1;
                                else
                                    stepPub{rxUser}(stream,symbol) = decisionFeedbackPub{rxUser}(stream,symbol) * stepsize;
                                end
                                decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol) + stepPub{rxUser}(stream,symbol);
%                                 decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol+1) / p2c(argErrPub{rxUser}(stream,symbol),1);
                            end
                        else
                            for stream = 1:rxAntennas(rxUser)
                                [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});
                            end
                        end
                        
                        % remove the effect of the decoded message from the
                        % received signal
                        publicInterference{rxUser}(:,symbol) = H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}(:,symbol);
                    else
                        publicInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);

                    end

                    if (max(dofSplitPri{rxUser}) > 0)   % If private streams from user A received
                        display(['Decoding private stream from user ' num2str(rxUser) '.']);
                        isolatedPri{rxUser}(:,symbol) = (1/sqrt(SNR(rxUser,rxUser)) * receivedMessage{rxUser}(:,symbol)) - publicInterference{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                        equalisedPri{rxUser}(:,symbol) =  pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * isolatedPri{rxUser}(:,symbol);

%                         equalisedPri{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPri{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the private message from the desired user A
                        
                        if (DFEon)
                            for stream = 1:rxAntennas(rxUser)
                                equalisedPri{rxUser}(stream,symbol) = equalisedPri{rxUser}(stream,symbol) * decisionFeedbackPri{rxUser}(stream,symbol);
                                
                                [dataPri{rxUser}(stream,symbol), decodedPri{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(stream,symbol),codebookPri{rxUser}{stream});

                                argErrPri{rxUser}(stream,symbol) = argd(decodedPri{rxUser}(stream,symbol)) - argd(equalisedPri{rxUser}(stream,symbol));
                                magErrPri{rxUser}(stream,symbol) = abs(decodedPri{rxUser}(stream,symbol)) / abs(equalisedPri{rxUser}(stream,symbol));
                                if magErrPri{rxUser}(stream,symbol) < 1
                                    stepPri{rxUser}(stream,symbol) = decisionFeedbackPri{rxUser}(stream,symbol) * stepsize * -1;
                                else
                                    stepPri{rxUser}(stream,symbol) = decisionFeedbackPri{rxUser}(stream,symbol) * stepsize;
                                end
                                decisionFeedbackPri{rxUser}(stream,symbol+1) = decisionFeedbackPri{rxUser}(stream,symbol) + stepPri{rxUser}(stream,symbol);
%                                 decisionFeedbackPub{rxUser}(stream,symbol+1) = decisionFeedbackPub{rxUser}(stream,symbol+1) / p2c(argErrPub{rxUser}(stream,symbol),1);
                            end
%                         else
%                             for stream = 1:rxAntennas(rxUser)
%                                 [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});
%                             end
                        end

%                         for stream = 1:rxAntennas(rxUser)
%                             equalisedPri{rxUser}(stream,symbol) = scaleFactor{rxUser}(stream) * equalisedPri{rxUser}(stream,symbol);
%                         end
                        
                        for i = 1:length(equalisedPri{rxUser}(:,symbol))
                            [dataPri{rxUser}(i,symbol), decodedPri{rxUser}(i,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(i,symbol),codebookPri{rxUser}{i});
                        end

                    end
                end         % end of SNR/INR testing loop
            end
        end             % end of txUser loop  

    end                 % end of rxUser loop
end                     % end of symbol loop

%% Plotting Results

rows = users * 2;
cols = max(rxAntennas);

figure;
for user = 1:users
    for stream = 1:rxAntennas(user)
        position = (user-1)*cols + stream;
        if ~isempty(equalisedPub{user})
            subplot(rows,cols,position);
            scatter(real(equalisedPub{user}(stream,:)),imag(equalisedPub{user}(stream,:)),'.');
            hold on;
            scatter(real(codebookPub{user}{stream}),imag(codebookPub{user}{stream}),'o','r');
            hold off;
            title(['User ' num2str(user) ' Public Stream ' num2str(stream) ' (M=' num2str(2^codebookIndexPub{user}(stream)) ')']);
        end
        if ~isempty(equalisedPri{user})        
            position = position + (users *cols);
            subplot(rows,cols,position);
            scatter(real(equalisedPri{user}(stream,:)),imag(equalisedPri{user}(stream,:)),'.');
            hold on;
            scatter(real(codebookPri{user}{stream}),imag(codebookPri{user}{stream}),'o','r');
            hold off;
            title(['User ' num2str(user) ' Private Stream ' num2str(stream) ' (M=' num2str(2^codebookIndexPri{user}(stream)) ')']);        
        end
    end
end

    
%% Calculate BER for each stream

publicSER = zeros(users,1);
privateSER = zeros(users,1);

deltaPri = cell(users,1);
deltaPub = cell(users,1);

for user = 1:users
    streams = txAntennas(user)-spareDimensions(user);
    if ~isempty(dataPri{user})
        deltaPri{user} = dataPri{user}(:,trainingSymbols+1:end)-privateCodeword{user}(:,trainingSymbols+1:end);
        privateSER(user) = length(nonzeros(deltaPri{user}))/(streams*totalSymbols);
    end
    if ~isempty(dataPub{user})
        deltaPub{user} = dataPub{user}(1:streams,trainingSymbols+1:end)-publicCodeword{user}(1:streams,trainingSymbols+1:end);
        publicSER(user) = length(nonzeros(deltaPub{user}))/(streams*totalSymbols);
    end
end
     
condition
publicSER
privateSER