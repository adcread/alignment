%% blockalignment.m - performs alignment on a block of data, length determined by the input parameters.


%% Initialisation Phase

addpath('Adaptive Modulation');
addpath('DoF Calculation');
addpath('Equalisation');

dataSymbols = 1024;                              % number of data symbols to transmit
trainingSymbols = 0;                           % number of training symbols tranmitted before the data.

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
dofSplitPub{1} = [0.4 0.4 0];                                               % PARAMETER TO CHANGE 
dofSplitPri{1} = [0.0 0.0 0];                                               %%%%%%%%%%%%%%%%%%%%%%
dofSplitPub{2} = [0.4 0.4];
dofSplitPri{2} = [0.0 0.0];


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
    
for symbol = 1:trainingSymbols                                              % Generate codewords and symbols for training sequence
    for user = 1:users
        for stream=1:txAntennas(user)
            privateCodeword{user}(stream,symbol) = 1;                       % DUMMY VALUES TO TEST STRUCTURES
            privateSymbol{user}(stream,symbol) = 0;
        end
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
    
for user = 1:users
    scaleFactor{user} = zeros(1,txAntennas(user));
    for stream = 1:txAntennas(user)
        scaleFactor{user}(stream) = scaleFactor{user}(stream) + mean(real(transmittedMessage{user}(stream,:)).^2 +imag(transmittedMessage{user}(stream,:)).^2);
        if scaleFactor{user}(stream)~=0
        transmittedMessage{user}(stream,:) = transmittedMessage{user}(stream,:)/scaleFactor{user}(stream);
        end
    end
end

%% Passing through the channel


for rxUser = 1:users
    receivedMessage{rxUser} = zeros(rxAntennas(rxUser),totalSymbols);
    for txUser = 1:users
        receivedMessage{rxUser} = receivedMessage{rxUser} + (sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * transmittedMessage{txUser});
    end
end

for rxUser = 1:users
    receivedMessage{rxUser} = receivedMessage{rxUser} + circSymAWGN(rxAntennas(rxUser),totalSymbols,1); % add AWGN to received signal
end


%% Equalisation - not currently used

for symbol = 1:trainingSymbols                                              % Perform adaptive equalisation in this section
    
end

for symbol = (trainingSymbols+1):totalSymbols                               % Use calculated Equaliser to perform equalisation in this section

end


%% Detection & Estimation


for symbol = 1:totalSymbols

    for rxUser = 1:users

        publicInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
        crossInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);

        commonSubspace = cell(users,1);

        for txUser = 1:users
            if (rxUser ~= txUser)

                if (SNR(rxUser,rxUser) >= SNR(txUser,rxUser))

                    % if the SNR of the direct stream is greater than the INR
                    % from user B, decode the common message from user A first

                    commonSubspace{rxUser} = receivedMessage{rxUser}(:,symbol);

                    if (max(dofSplitPub{txUser}) > 0) % if public stream has been sent from the interfering user B

                        if (rxAntennas(rxUser) > txAntennas(txUser))

                            % If user A has more antennas then strip out the extra stream(s) so that
                            % user B's interference can be subtracted

                            numberOfDimensions = rxAntennas(rxUser);
                            hiddenDimensions = rxAntennas(rxUser) - txAntennas(txUser);
                            commonSubspace{rxUser} = removeHiddenStreams(commonSubspace{rxUser},(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
                        end

                        % Equalise user B's interference using SVD of the cross channel

                        equalisedInt{rxUser}(:,symbol) = pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * commonSubspace{rxUser};

                        equalisedInt{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedInt{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet


                        % decode the public message from the unwanted user B as interference

                        for stream = 1:txAntennas(txUser)
                        [dataInt{rxUser}(stream,symbol), decodedInt{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(stream),codebookIndexPub{txUser}(stream));
                        end

                        % remove the effect of the decoded message from the received signal 

                        crossInterference{rxUser}(:,symbol) = H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * decodedInt{rxUser}(:,symbol);

                    end

                    commonSubspace{rxUser} = receivedMessage{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                    if (max(dofSplitPub{rxUser}) > 0) % If public messages sent by user A

                        equalisedPub{rxUser}(:,symbol) = pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser};

                        equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPub{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the message and symbols sent by user A

                        for stream = 1:length(equalisedPub{rxUser}(:,symbol))
                            [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookIndexPub{rxUser}(stream));
                        end
                        % remove the effect of the decoded message from the
                        % received signal

                        publicInterference{rxUser}(:,symbol) = H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}(:,symbol);                  

                    end

                    if (max(dofSplitPri{rxUser}) > 0)   % If private streams from user A received

                        isolatedPri{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - publicInterference{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);

                        equalisedPri{rxUser}(:,symbol) = pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * isolatedPri{rxUser}(:,symbol);

                        equalisedPri{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) * equalisedPri{rxUser}(:,symbol);       % Attenuate to bring constellation back to alphabet

                        % decode the private message from the desired user A

                        for i = 1:length(equalisedPri{rxUser}(:,symbol))
                            [dataPri{rxUser}(i,symbol), decodedPri{rxUser}(i,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(i,symbol),codebookIndexPri{rxUser}(i));
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
        if ~isempty(equalisedPub{user})
            position = (user-1)*cols + stream;
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

publicBER = zeros(users,1);
privateBER = zeros(users,1);

for user = 1:users
    streams = txAntennas(user)-spareDimensions(user);
    if ~isempty(dataPri{user})
        delta = dataPri{user}(1:streams,:)-privateCodeword{user}(1:streams,:);
        privateBER(user) = 1 - sum(delta(:)==0)/(streams*totalSymbols)
    end
    if ~isempty(dataPub{user})
        delta = dataPub{user}(1:streams,:)-publicCodeword{user}(1:streams,:);
        publicBER(user) = 1 - sum(delta(:)==0)/(streams*totalSymbols)
    end
end
     
publicBER
privateBER