%% Kronecker channel model decomposition

addpath('C:\PhD\alignment\Channel');
addpath('C:\PhD\alignment\Equalisation');
addpath('C:\PhD\alignment\General functions');

% For each antenna array create the correlation matrix for a uniform 
% linear array of antennas txDistance wavelengths apart

txAntennas = 4;
txDistance = 0.1;
txCorrelation = arrayCorrelation(txAntennas,txDistance);

rxAntennas = 3;
rxDistance = 0.3;
rxCorrelation = arrayCorrelation(rxAntennas,rxDistance);

% Structure of the training sequence:
% sequenceLength = length of random training signal
% noRepetitions = number of times random signal repeated (creates time corr)
% blockLength = sequenceLength * noRepetitions
% noBlocks = number of blocks sent (allows expectation over ensemble)

sequenceLength = 64;

noRepetitions = 1;

blockLength = sequenceLength * noRepetitions;

noBlocks = 256;

% Use the Kronecker channel model to generate a channel given the
% correlation matrices

H = KroneckerChannel(txAntennas,rxAntennas,txCorrelation,rxCorrelation);

channelPower = 1 ./ svd(H);

% Create array to hold entirety of signal

transmittedSequence = cell(noBlocks,1);
transformedSequence = cell(noBlocks,1);
receivedSequence = cell(noBlocks,1);

transmittedSequenceCorrelation = cell(noBlocks,1);
transformedSequenceCorrelation = cell(noBlocks,1);

y = cell(noBlocks,1);
z = cell(noBlocks,1);


%% Create the training sequence to be repeated in each block
% (dimensions are L x M)

for block = 1:noBlocks
    
    transmittedSequence{block} = zeros((blockLength), txAntennas);
    y{block} = zeros(blockLength*rxAntennas,1);
    z{block} = zeros(blockLength*rxAntennas);
  
    transmittedSequenceCorrelation{block} = zeros((blockLength),1);
    transformedSequenceCorrelation{block} = zeros((blockLength),1);
    transmittedSequenceCorrelation_temp = zeros((blockLength*2)-1,1);
    transformedSequenceCorrelation_temp = zeros((blockLength*2)-1,1);
    
%      sequence = circSymAWGN(sequenceLength,txAntennas(1),1);% * diag(channelPower);
%     for stream = 1:txAntennas
%         sequence(:,stream) = lteZadoffChuSeq(1,sequenceLength+1);
%     end
    % Repeat the sequence to create temporal correlations.

    % Calculate mean of Tx correlation matrix
      
    upperTriangle = chol(inv(txCorrelation));
    
    sequence = (upperTriangle' * (randn(txAntennas,sequenceLength) + 1i*randn(txAntennas,sequenceLength))/sqrt(2))';
    
%     sequence = circSymAWGN(sequenceLength,txAntennas(1),1);% * diag(channelPower);

    repeatedSequence = repmat(sequence,noRepetitions,1);

    % add the repeated sequence to the transmitted sequence array.
    
    transmittedSequence{block} = repeatedSequence;
    transformedSequence{block} = transmittedSequence{block} * sqrtm(txCorrelation);  
   
    % pass the signal through the channel

    receivedSequence{block} = H*transmittedSequence{block}.' + circSymAWGN(rxAntennas,blockLength,1);
    
    %vectorise Y to y

    y{block} = vec(receivedSequence{block});

    % Compute the covariance matrix of the received signal

    z{block} = z{block} + y{block} * y{block}';

    % Compute the autocorrelation over the ensemble of blocks

    for stream = 1:txAntennas
        transmittedSequenceCorrelation_temp = transmittedSequenceCorrelation_temp + xcorr(transmittedSequence{block}(:,stream));
        transformedSequenceCorrelation_temp = transformedSequenceCorrelation_temp + xcorr(transformedSequence{block}(:,stream));
    end

    % remove the additional autocorrelation entries and save

    transmittedSequenceCorrelation{block} = transmittedSequenceCorrelation_temp(blockLength:end);
    transformedSequenceCorrelation{block} = transformedSequenceCorrelation_temp(blockLength:end);

end

%% Work out the correlations over the ensemble

transmittedSequenceCorrelation_ave = ensembleAve(transmittedSequenceCorrelation);
transformedSequenceCorrelation_ave = ensembleAve(transformedSequenceCorrelation);

%% Calculate Q, the original sequence autocorrelation matrix

Q = zeros(blockLength);

Q(1,:) = transmittedSequenceCorrelation_ave.';

b = cell(noBlocks,1);

for i = 2:(blockLength)
   Q(i,:) = [circshift(transmittedSequenceCorrelation_ave.',i-1,2)];
end

% Calculate Q from the actual signal

[m,n] = size(rxCorrelation);

for block = 1:noBlocks
    
    b{block} = zeros(sequenceLength);

    % Method for calculating b comes from [79]

    for i = 1:(blockLength)
        for j = 1:(blockLength)
            b{block}(i,j) = trace(z{block}(((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)).'*rxCorrelation) / trace(rxCorrelation.'*rxCorrelation);
        end
    end
    
end

Q_calc = ensembleAve(b);


%% Plot the autocorrelations of the two sequences

figure;
hold on;
plot([0:(blockLength)-1],abs(transmittedSequenceCorrelation_ave),'r');
plot([0:(blockLength)-1],abs(transformedSequenceCorrelation_ave),'b');


figure();

subplot(3,1,1)
plot(abs(transmittedSequenceCorrelation{1}));
title('Training Sequence Autocorrelation');

subplot(3,1,2)
plot(abs(transformedSequenceCorrelation{1}));
title('Transformed Sequence Autocorrelation');

subplot(3,1,3)

figure();
imagesc(abs(Q));
title('Calculated Q Matrix');

figure();
imagesc(abs(Q_calc));
title('Computed b Matrix');