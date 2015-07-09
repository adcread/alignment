%% blockalignment.m - performs alignment on a block of data, length determined by the input parameters.


%% Initialisation Phase

addpath('Adaptive Modulation');
addpath('DoF Calculation');
addpath('Equalisation');
addpath('General functions');
addpath('Channel');

dataSymbols = 1024;                              % number of data symbols to transmit
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
noiseOn = true;
generateNewChannel = true;
displayConstellations = true;

%% Creation of Channel matrices

if generateNewChannel
    
    H = generateChannel(users, txAntennas, rxAntennas, 'kronecker');
%     H = generateChannel(users, txAntennas, rxAntennas, 'gaussian');

    H{2,1} = H{1,2}'; % Channel responses are reciprocal
end

[U, sigma, V] = eigenchannel(H);  % Eigenchannel decomposition

R0 = cell(users,1);
R1 = cell(users,1);

eigenspaceCorrelation = zeros(users,1);

for rxUser = 1:users

    R0{rxUser} = H{rxUser,rxUser} * H{rxUser,rxUser}';
    R1{rxUser} = zeros(rxAntennas(rxUser));
        
    for txUser = 1:users
        condition(rxUser,txUser) = max(nonzeros(max(sigma{rxUser,txUser})))/min(nonzeros(max(sigma{rxUser,txUser})));
        if (rxUser ~= txUser)
            R1{rxUser} = R1{rxUser} + H{txUser,rxUser}*H{txUser,rxUser}';
        end
    end
    
    eigenspaceCorrelation(rxUser) = J(R0{rxUser},R1{rxUser},1);
    
end

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
dofSplitPri{1} = [0.0 0.0 0.7];                                               %%%%%%%%%%%%%%%%%%%%%%
dofSplitPub{2} = [0.2 0.2];
dofSplitPri{2} = [0.4 0.4];

%% Creation of Source Alphabets

codebookPri = cell(users,1);
codebookPub = cell(users,1);

for user = 1:users
    codebookPri{user} = cell(txAntennas(user),1);
    codebookPub{user} = cell(txAntennas(user),1);
    
    publicStreams(user) = length(nonzeros(dofSplitPub{user}));
    privateStreams(user) = length(nonzeros(dofSplitPri{user}));
    
end

for user = 1:users
    for stream = 1:txAntennas(user)
       consoleOutput(displayConstellations,['Generating Private codebook for User ' num2str(user) ' stream ' num2str(stream) ':']);
       [codebookPri{user}{stream},codebookIndexPri{user}(stream)] = generateTransmitCodebook(dofSplitPri{user}(stream), SNR(user,user),6);
       consoleOutput(displayConstellations,['Generating Public codebook for User ' num2str(user) ' stream ' num2str(stream) ':']);
       [codebookPub{user}{stream},codebookIndexPub{user}(stream)] = generateTransmitCodebook(dofSplitPub{user}(stream), SNR(user,user),6);      
       consoleOutput(displayConstellations,[]);
    end
end

%% Creation of messages to transmit

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
            directionPub{i,j} = inv(V{i,j}) * covariancePub{i} * inv(V{i,j}');
            if rxAntennas(j) < txAntennas(i)
                b = txAntennas(i) - spareDimensions(i)+1;
                directionPub{i,j}(b,b) = 0;
            end    
            directionPri{i,j} = real(diag(diag(directionPri{i,j})));
            directionPub{i,j} = real(diag(diag(directionPub{i,j})));
        end
    end
end

% for i=1:users
%     for j = 1:users
%         if (i ~= j)
%             directionPri{i,j} = zeros(txAntennas(i));
%             directionPub{i,j} = zeros(txAntennas(i));
%             
%             Lam = diag(sigma{i,j}.^2);
%             
%             for k = 1:txAntennas(i)
%                 if k <= min(txAntennas(i),rxAntennas(j))
%                     directionPri{i,j}(k,k)= (1/txAntennas(i)) * 1/(1 + SNR(i,j) * Lam(k)^2);
%                     directionPub{i,j}(k,k) = (1/txAntennas(i)) - directionPri{i,j}(k,k);
%                 else
%                     directionPri{i,j}(k,k) = (1/txAntennas(i));
%                 end
%             end
%         end
%     end
% end


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
        publicCodeword{user}(stream,:) = randi(MPub,1,totalSymbols);
        privateCodeword{user}(stream,:) = randi(MPri,1,totalSymbols);
    end

end

for user = 1:users
    for stream = 1:txAntennas(user)
        privateSymbol{user}(stream,:) = codebookPri{user}{stream}(privateCodeword{user}(stream,:));
        publicSymbol{user}(stream,:) = codebookPub{user}{stream}(publicCodeword{user}(stream,:));   
    end
end
  
    privateMessage{1} = V{1,2} * sqrt(directionPri{1,2}) * privateSymbol{1};
    privateMessage{2} = V{2,1} * sqrt(directionPri{2,1}) * privateSymbol{2};

    publicMessage{1} = V{1,2} * sqrt(directionPub{1,2}) * publicSymbol{1};
    publicMessage{2} = V{2,1} * sqrt(directionPub{2,1}) * publicSymbol{2};

    transmittedMessage{1} = privateMessage{1} + publicMessage{1};
    transmittedMessage{2} = privateMessage{2} + publicMessage{2};
    
%% Calculate transmitter input power and normalise to <=1

% Q = cell(totalSymbols,1);
% G = cell(totalSymbols,1);
% 
% for user = 1:users
%     G{1} = zeros(txAntennas(user));
%     for symbol = 2:length(transmittedMessage{user})
%         Q{symbol} = transmittedMessage{user}(:,symbol) * transmittedMessage{user}(:,symbol)';
%         G{symbol} = zeros(txAntennas(user));
%         for i = 1:length(Q{symbol})
%             for j = 1:length(Q{symbol})
%                 G{symbol}(i,j) = G{symbol-1}(i,j) + Q{symbol}(i,j);
%             end
%         end
%         T(symbol) = trace(G{symbol});
%         T(symbol) = T(symbol)/symbol;
%     end
%     
%     scaleFactor(user) = sqrt(T(totalSymbols-1));
%     
%     if (scaleFactor(user)>0)
%         transmittedMessage{user} = transmittedMessage{user} / scaleFactor(user);
%     end
%     
%     end
    
%% Passing through the channel


for rxUser = 1:users
    receivedMessage{rxUser} = zeros(rxAntennas(rxUser),totalSymbols);
    for txUser = 1:users
        
        % received message is superposition of messages received from each user
        
        receivedMessage{rxUser} = receivedMessage{rxUser} + (sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * transmittedMessage{txUser});
    end
end

if noiseOn == true
    for rxUser = 1:users
        
        % add AWGN to received signal
        
        receivedMessage{rxUser} = receivedMessage{rxUser} + circSymAWGN(rxAntennas(rxUser),totalSymbols,1); 
    end
end

% remove scaling factor from input power constraint

% for rxUser = 1:users
%     receivedMessage{rxUser} = receivedMessage{rxUser} * scaleFactor(rxUser);
% end

%% Equalisation & Detection

for symbol = (trainingSymbols+1):totalSymbols
    
    commonSubspace = cell(users,1);
    
    for rxUser = 1:users      

        for txUser = 1:users
            
            % only start to decode when identified interfering users
            
            if (rxUser ~= txUser)
                
                % if desired public message is interference-free (i.e.
                % public DoF 'fits' between alpha(1,1) and alpha(2,1))
                % decoding order is then A_pub, B_pub (if needed) , A_pri
                
                commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol);
                
                if alpha(rxUser,rxUser) - max(dofSplitPub{rxUser}) >= alpha(txUser,rxUser)
                   
                    % if a desired public stream has been sent
                        
                    if (publicStreams(rxUser) > 0)
                
                        consoleOutput(displayConstellations,['Decoding Public message from User ' num2str(rxUser) '.']);
                                           
                        equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(rxUser,rxUser)) * pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);
                        
                        % decode the public stream
                        
                        for stream = 1:length(equalisedPub{rxUser}(:,symbol))
                            [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});
                        end
                        
                        % create an estimate of the public stream component
                        % of the received signal for cancellation
                        
                        publicInterference{rxUser}(:,symbol) = sqrt(SNR(rxUser,rxUser)) * H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}(:,symbol);
                   
                    % otherwise there is no public interference to cancel out -
                    % set the cancellation vector to zeros
                        
                    else
                        publicInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end

                    % cancel the effect of the desired public stream from
                    % the received signal
                    
                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - publicInterference{rxUser}(:,symbol);
                    
                    % if the interfering user has sent a public stream
                
                    if (publicStreams(txUser) > 0) 

                        consoleOutput(displayConstellations,['Decoding Interference from User ' num2str(txUser) '.']);
                                           
                        % if desired private message is higher-dimension
                        % than undesired public message
                   
%                         if (txAntennas(rxUser) > txAntennas(txUser))
                                                       
%                             hiddenDimensions = txAntennas(rxUser) - txAntennas(txUser);
                            
                            % remove the hidden streams
                            
%                             isolatedInt{rxUser}(:,symbol) = removeHiddenStreams(commonSubspace{rxUser}(:,symbol),(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
%                         else
                            isolatedInt{rxUser}(:,symbol) = commonSubspace{rxUser}(:,symbol);
%                         end
                        
                        
                        % Equalise interfering user's interference using SVD 
                        % of the cross channel
                        
                        equalisedInt{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) *  pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * isolatedInt{rxUser}(:,symbol);
                    
                        % decode the interfering public stream
                    
                        for stream = 1:length(equalisedInt{rxUser}(:,symbol))
                            [dataInt{rxUser}(stream,symbol), decodedInt{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(stream,symbol),codebookPub{txUser}{stream});
                        end
                    
                        % create an estimate of the public stream component
                        % of the received signal for cancellation
                    
                        crossInterference{rxUser}(:,symbol) = sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * decodedInt{rxUser}(:,symbol);
                    
                    % otherwise there is no public interference to cancel out -
                    % set the cancellation vector to zeros  
                    
                    else
                        crossInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end
                    
                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - publicInterference{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);
                
                    % if a desired private stream has been sent
                
                    if (privateStreams(rxUser) > 0)
                    
                        consoleOutput(displayConstellations,['Decoding private stream from User ' num2str(rxUser) '.']);
                    
                        % equalise the private stream using channel inversion
                    
                        equalisedPri{rxUser}(:,symbol) =  1/sqrt(SNR(rxUser,rxUser)) * pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);
                    
                        % decode private stream
                    
                        for stream = 1:length(equalisedPri{rxUser}(:,symbol))
                            [dataPri{rxUser}(stream,symbol), decodedPri{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(stream,symbol),codebookPri{rxUser}{stream});
                        end                
                    
                    end
                
                    % otherwise the public interfernce from undesired user
                    % mustbe removed first - decoding order becomes B_pub,
                    % A_pub, A_Pri
                
                elseif alpha(rxUser,rxUser) - dofSplitPub{rxUser}(1) < alpha(txUser,rxUser)
                    
                    % if a desired private stream has been sent
                    
                    if (publicStreams(txUser) > 0) 

                        consoleOutput(displayConstellations,['Decoding Interference from User ' num2str(txUser) '.']);
                                           
                        % it may be necessary to eliminate the private
                        % stream sent at full power from the desired user
                        
                        if (privateStreams(rxUser)>0) && (txAntennas(rxUser) > txAntennas(txUser))
                            
                            % hidden dimensions = number of tx antennas
                            % that user A has more than B
                            
                            hiddenDimensions = txAntennas(rxUser) - txAntennas(txUser);
                            
                            % remove the hidden dimension
                            
                            commonSubspace{rxUser}(:,symbol) = removeHiddenStreams(commonSubspace{rxUser}(:,symbol),(H{rxUser,rxUser}*V{rxUser,txUser}),hiddenDimensions);
                        end
                        
                        % Equalise interfering user's interference using SVD 
                        % of the cross channel
                   
                        equalisedInt{rxUser}(:,symbol) = 1/sqrt(SNR(txUser,rxUser)) *  pinv(sqrt(directionPub{txUser,rxUser})) * pinv(sigma{txUser,rxUser}) * U{txUser,rxUser}' * commonSubspace{rxUser}(:,symbol);
                    
                        % decode the interfering public stream
                    
                        for stream = 1:length(equalisedInt{rxUser}(:,symbol))
                            [dataInt{rxUser}(stream,symbol), decodedInt{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedInt{rxUser}(stream,symbol),codebookPub{txUser}{stream});
                        end
                    
                        % create an estimate of the public stream component
                        % of the received signal for cancellation
                    
                        crossInterference{rxUser}(:,symbol) = sqrt(SNR(txUser,rxUser)) * H{txUser,rxUser} * V{txUser,rxUser} * sqrt(directionPub{txUser,rxUser}) * decodedInt{rxUser}(:,symbol);
                    
                    % otherwise there is no public interference to cancel out -
                    % set the cancellation vector to zeros  
                    
                    else
                        crossInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end
                    
                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);
                    
                    % if a desired public stream has been sent
                    
                    if (publicStreams(rxUser) > 0)
                
                        consoleOutput(displayConstellations,['Decoding Public message from User ' num2str(rxUser) '.']);
                                           
                        equalisedPub{rxUser}(:,symbol) = 1/sqrt(SNR(rxUser,rxUser)) * pinv(sqrt(directionPub{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);
                        
                        % decode the public stream
                        
                        for stream = 1:length(equalisedPub{rxUser}(:,symbol))
                            [dataPub{rxUser}(stream,symbol), decodedPub{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPub{rxUser}(stream,symbol),codebookPub{rxUser}{stream});
                        end
                        
                        % create an estimate of the public stream component
                        % of the received signal for cancellation
                        
                        publicInterference{rxUser}(:,symbol) = sqrt(SNR(rxUser,rxUser)) * H{rxUser,rxUser} * V{rxUser,txUser} * sqrt(directionPub{rxUser,txUser}) * decodedPub{rxUser}(:,symbol);
                   
                    % otherwise there is no public interference to cancel out -
                    % set the cancellation vector to zeros
                        
                    else
                        publicInterference{rxUser}(:,symbol) = zeros(rxAntennas(rxUser),1);
                    end                    
                    
                    commonSubspace{rxUser}(:,symbol) = receivedMessage{rxUser}(:,symbol) - publicInterference{rxUser}(:,symbol) - crossInterference{rxUser}(:,symbol);
                    
                    if (privateStreams(rxUser) > 0)
                    
                        consoleOutput(displayConstellations,['Decoding private stream from User ' num2str(rxUser) '.']);
                    
                        % equalise the private stream using channel inversion
                    
                        equalisedPri{rxUser}(:,symbol) =  1/sqrt(SNR(rxUser,rxUser)) * pinv(sqrt(directionPri{rxUser,txUser})) * V{rxUser,txUser}' * pinv(H{rxUser,rxUser}) * commonSubspace{rxUser}(:,symbol);
                    
                        % decode private stream
                    
                        for stream = 1:length(equalisedPri{rxUser}(:,symbol))
                            [dataPri{rxUser}(stream,symbol), decodedPri{rxUser}(stream,symbol)] = maximumLikelihoodQAMDecoder(equalisedPri{rxUser}(stream,symbol),codebookPri{rxUser}{stream});
                        end                
                    
                    end
                end
            end
        end
    end
end


%% Plotting Results
if (displayConstellations)
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
    streams = txAntennas(user)-spareDimensions(user);
    if ~isempty(dataPub{user})
        deltaPub{user} = dataPub{user}(1:streams,trainingSymbols+1:end)-publicCodeword{user}(1:streams,trainingSymbols+1:end);
        publicSER(user) = length(nonzeros(deltaPub{user}))/(streams*totalSymbols);
    end
end
     
if (displayConstellations)
    condition
    publicSER
    privateSER
end
