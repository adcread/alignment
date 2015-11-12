%% Kronecker channel model decomposition

addpath('C:\PhD\alignment\Channel');
addpath('C:\PhD\alignment\Equalisation');
addpath('C:\PhD\alignment\General functions');

% Structure of the training sequence:
% sequenceLength = length of random training signal
% noRepetitions = number of times random signal repeated (creates time corr)
% noBlocks = number of blocks sent (allows expectation over ensemble)

sequenceLength = 32;

noRepetitions = 4;

noBlocks = 128;

% For each antenna array create the correlation matrix for a uniform 
% linear array of antennas txDistance wavelengths apart

txAntennas = 3;
txDistance = 0.1;
txCorrelation = arrayCorrelation(txAntennas,txDistance);


rxAntennas = 5;
rxDistance = 0.3;
rxCorrelation = arrayCorrelation(rxAntennas,rxDistance);

% Use the Kronecker channel model to generate a channel given the
% correlation matrices

H = KroneckerChannel(txAntennas,rxAntennas,txCorrelation,rxCorrelation);

channelPower = 1 ./ svd(H);

% Create array to hold entirety of signal

transmittedSequence = zeros((sequenceLength*noRepetitions*noBlocks), txAntennas);

%% Creates a AWGN signal to be used as non-optimal training sequence
% (dimensions are L x M)

for block = 1:noBlocks
    
    sequence = circSymAWGN(sequenceLength,txAntennas(1),1);% * diag(channelPower);

    % Repeat the sequence to create temporal correlations.

    repeatedSequence = sequence;

    for i = 2:noRepetitions
        repeatedSequence = [repeatedSequence; sequence];
    end

    % add the repeated sequence to the transmitted sequence array.
    
    transmittedSequence(sequenceLength*noRepetitions*(block-1)+1:sequenceLength*noRepetitions*block,:) = repeatedSequence;
    
end
    
% multiply by the transmit correlation matrix to determine the transformed
% signal symbol matrix (needed to calculate Q).

transformedSequence = transmittedSequence * sqrtm(txCorrelation);

%% pass the signal through the channel

receivedSequence = H*transmittedSequence.';

%% Repartition received signal back into blocks for processing

y = zeros(sequenceLength*noRepetitions,1);
z = zeros(sequenceLength*noRepetitions*rxAntennas);

for block = 1:noBlocks
    
    startIndex = ((block-1) * sequenceLength * noRepetitions) + 1;
    endIndex = block * sequenceLength * noRepetitions;
    
    %vectorise Y to y

    y = vec(receivedSequence(:,startIndex:endIndex));

% Compute the covariance matrix of the received signal

    z = z + y * y';
end

%% Compute the autocorrelation over the ensemble of blocks

transmittedSequenceCorrelation = zeros((sequenceLength*noRepetitions*2)-1,1);
transformedSequenceCorrelation = zeros((sequenceLength*noRepetitions*2)-1,1);

for block = 1:noBlocks

    for stream = 1:txAntennas
        transmittedSequenceCorrelation = transmittedSequenceCorrelation + xcorr(transmittedSequence((sequenceLength*noRepetitions*(block-1)+1:sequenceLength*noRepetitions*block),stream),'coeff');
        transformedSequenceCorrelation = transformedSequenceCorrelation + xcorr(transformedSequence((sequenceLength*noRepetitions*(block-1)+1:sequenceLength*noRepetitions*block),stream),'coeff');
    end
end

% remove the additional autocorrelation entries and normalise

transmittedSequenceCorrelation = transmittedSequenceCorrelation(sequenceLength*noRepetitions:end);
transformedSequenceCorrelation = transformedSequenceCorrelation(sequenceLength*noRepetitions:end);

%% Plot the autocorrelations of the two sequences

figure;
hold on;
plot([0:(sequenceLength*noRepetitions)-1],abs(transmittedSequenceCorrelation),'r');
plot([0:(sequenceLength*noRepetitions)-1],abs(transformedSequenceCorrelation),'b');

%% Calculate Q, the original sequence autocorrelation matrix

Q = zeros(sequenceLength*noRepetitions,sequenceLength*noRepetitions);

Q(1,:) = transmittedSequenceCorrelation.';

for i = 2:(sequenceLength*noRepetitions)
   Q(i,:) = [circshift(transmittedSequenceCorrelation.',i-1,2)];
end

% Calculate Q from the actual signal

[m,n] = size(rxCorrelation);

b = zeros(sequenceLength);

% Method for calculating b comes from [79]

for i = 1:(sequenceLength*noRepetitions)
    for j = 1:(sequenceLength*noRepetitions)
        b(i,j) = trace(z(((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)).'*rxCorrelation) / trace(rxCorrelation.'*rxCorrelation);
    end
end
        
figure();

subplot(3,1,1)
plot(abs(transmittedSequenceCorrelation));
title('Training Sequence Autocorrelation');

subplot(3,1,2)
plot(abs(transformedSequenceCorrelation));
title('Transformed Sequence Autocorrelation');

subplot(3,1,3)


figure();
imagesc(abs(Q));

figure();
imagesc(abs(b));
