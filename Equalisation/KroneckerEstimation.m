%% Kronecker channel model decomposition

if ispc
    addpath('C:\PhD\alignment\Channel');
    addpath('C:\PhD\alignment\Equalisation');
    addpath('C:\PhD\alignment\General functions');

elseif ismac
    addpath('/Users/chris/PhD/alignment/Channel');
    addpath('/Users/chris/PhD/alignment/Equalisation');
    addpath('/Users/chris/PhD/alignment/General functions');   
end

% For each antenna array create the correlation matrix for a uniform 
% linear array of antennas txDistance wavelengths apart

% txAntennas = 6;
txDistance = 0.4;
txCorrelation = arrayCorrelation(txAntennas,txDistance);
% txCorrelation = eye(txAntennas);

% rxAntennas = 3;
rxDistance = 0.4;
rxCorrelation = arrayCorrelation(rxAntennas,rxDistance);
% rxCorrelation = eye(rxAntennas);

% Structure of the training sequence:
% sequenceLength = length of random training signal
% noRepetitions = number of times random signal repeated (creates time corr)
% blockLength = sequenceLength * noRepetitions
% noBlocks = number of blocks sent (allows expectation over ensemble)

% sequenceLength = 63;

noRepetitions = 1;

blockLength = sequenceLength * noRepetitions;

% Use the Kronecker channel model to generate a channel given the
% correlation matrices

H = [];

H = generateChannel(txAntennas,rxAntennas,'kronecker',txCorrelation,rxCorrelation);

% H = generateChannel(txAntennas,rxAntennas,'weichselberger',couplingMatrix);

% Create arrays to hold entirety of signal

sequence = [];
transmittedSequence = cell(noBlocks,1);
transformedSequence = cell(noBlocks,1);
receivedSequence = cell(noBlocks,1);

transmittedSequenceAutocorrelation = cell(noBlocks,1);
transformedSequenceAutocorrelation = cell(noBlocks,1);

transmittedSequenceCrosscorrelation = cell(noBlocks,1);
transformedSequenceCrosscorrelation = cell(noBlocks,1);

y = cell(noBlocks,1);
z = cell(noBlocks,1);

Q = cell(noBlocks,1);

%% Create the training sequence to be repeated in each block
% (dimensions are L x M)

for block = 1:noBlocks
    
    transmittedSequence{block} = zeros((blockLength), txAntennas);
    receivedSequence{block} = zeros(rxAntennas,blockLength);
    
    y{block} = zeros(blockLength*rxAntennas,1);
    z{block} = zeros(blockLength*rxAntennas);

    
    %% Create training sequences from Zadoff-Chu sequence (sequenceLength must be odd)
    
%     zadoffChuRoot = 49;
% 
%     zadoffChuSequence = lteZadoffChuSeq(zadoffChuRoot,sequenceLength);
%     
%     sequenceCrosscorrelation = abs(xcorr(zadoffChuSequence,'coeff'));
%     
%     [sortedCrosscorrelations, streamLags] = sort(sequenceCrosscorrelation(sequenceLength:end),'ascend');
%     
%     % Find the shifts that yield the lowest cross-correlation between
%     % sequences. This is as close to orthogal as ZC sequences will get!
% 
%     rootSequence = lteZadoffChuSeq(49,sequenceLength);
%     
%     for stream = 1:txAntennas
%         sequence(:,stream) = circshift(rootSequence,streamLags(stream));
%     end
    
    %% Create the training sequence with Walsh-Hadamard sequence (sequenceLength must be power of 2)
    
%     rootSequence = generateHadamardMatrix(2^(ceil(log2(sequenceLength))));
%     
%     for stream = 1:txAntennas
%         sequence(:,stream) = rootSequence(:,stream);
%     end
        
    %% Create the training sequence with AWGN, Choleksy decomposition of desired covariance matrix
    
    % Calculate mean of Tx correlation matrix
      
%     upperTriangle = chol(inv(txCorrelation));
%     
%     upperTriangle = eye(txAntennas);
%     
%     sequence = (upperTriangle' * (randn(txAntennas,sequenceLength) + 1i*randn(txAntennas,sequenceLength))/sqrt(2))';
%     
%     sequence = circSymAWGN(sequenceLength,txAntennas(1),1);% * diag(channelPower);    
    
    %% Generate a set of Gold codes and use them for training sequence
    
%    rootSequence = generateGoldCodes(round(log2(sequenceLength)));
    

    sequenceSelection = randsample(4:(sequenceLength+2),txAntennas);
    sequence = zeros(sequenceLength,txAntennas);
    for i = 1:txAntennas
        sequence(:,i) = rootSequence(:,sequenceSelection(i));
    end
    
    %% Repeat the sequence to create temporal correlations.

    repeatedSequence = repmat(sequence,noRepetitions,1);

    % add the repeated sequence to the transmitted sequence array.
    
    transmittedSequence{block} = repeatedSequence;% * pinv(sqrtm(txCorrelation));
    transformedSequence{block} = transmittedSequence{block} * sqrtm(txCorrelation);  
   
    %% Pass the signal through the channel

    receivedSequence{block} = SNR/txAntennas * H * transmittedSequence{block}.' + circSymAWGN(rxAntennas,blockLength,1);
    
    %vectorise Y to y

    y{block} = vec(receivedSequence{block});

    % Compute the covariance matrix of the received signal

    z{block} = y{block} * y{block}';
    
    Q{block} = transmittedSequence{block} * txCorrelation * transmittedSequence{block}';

end

%% Calculate Q, the original sequence autocorrelation matrix

% for block = 1:noBlocks
% 
%     Q{block} = zeros(blockLength);
% 
%     q = zeros(1,blockLength);
% 
%     for stream = 1:txAntennas
%         q = q + transmittedSequenceAutocorrelation{block}(stream,blockLength:end);
%     end
% 
%     Q{block}(1,:) = q;
% 
%     for i = 2:(blockLength)
%        Q{block}(i,:) = circshift(q,i-1,2);
%     end
%     
% end

% Calculate Q from the actual signal

[m,n] = size(rxCorrelation);

Z = ensembleAve(z);

b = zeros(blockLength);

% Method for calculating b comes from [79] - seperable least squares
% framework (theorem 3)
    

if (statisticsFlag)
    for i = 1:(blockLength)
        for j = 1:(blockLength)
            b(i,j) = trace(Z(((i-1)*m+1):(i*m),((j-1)*n+1):(j*n))'* rxCorrelation) / trace(rxCorrelation'*rxCorrelation);
        end
    end    
else
    for i = 1:(blockLength)
        b(i,i) = trace(Z(((i-1)*m+1):(i*m),((i-1)*n+1):(i*n))'* rxCorrelation) / trace(rxCorrelation'*rxCorrelation);
    end    
end

