function [ H ] = generateChannel( txAntennas, rxAntennas, channelModel, modelArg1, modelArg2, modelArg3 )
%GENERATECHANNEL Summary of this function goes here
%   Detailed explanation goes here

if strcmp(channelModel,'scalar')

elseif strcmp(channelModel, 'czink')

elseif strcmp(channelModel,'complex')

elseif strcmp(channelModel,'gaussian')

    for i = 1:users
        for j = 1:users
            H = normrnd(1,1,txAntennas,rxAntennas);
        end
    end

elseif strcmp(channelModel,'kronecker')

    % take in the Tx/Rx correlation matrices from the arguments
    
    txCorrelation = modelArg1;
    rxCorrelation = modelArg2;
    
    R = kron(txCorrelation,rxCorrelation); % Apply Kronecker product to resolve the full matrix
    
    %Generate a long vector form of NTx by NRx channels with N samples (note
    %the book only uses one sample

    h = (randn(rxAntennas * txAntennas,1) + 1j * randn(rxAntennas * txAntennas,1))/sqrt(2);

    % Either the sqrtm or Choelsky factorisation can be used here.
     W = sqrtm(R); 
    %W = chol(R);

    % Create Kronecker correlated channel H

    w = W * h;
    H = reshape(w, rxAntennas, txAntennas);
    
elseif strcmp(channelModel,'weichselberger')
    
    % take in the Coupling matrix as argument
    
    txCorrelation = modelArg1;
    rxCorrelation = modelArg2;
    couplingMatrix = modelArg3;
    
    % Create sets of eigenbases for the transmit and receive antennas
    
    [txEigenbase, ~, ~] = svd(txCorrelation);
    [rxEigenbase, ~, ~] = svd(rxCorrelation);
    
    H = rxEigenbase * (couplingMatrix .* randn(rxAntennas,txAntennas) + 1j * randn(rxAntennas,txAntennas)) * txEigenbase';
    
else
    % no operation, return empty cell of channels
    
end

    power = trace(H * H');
    H = H * sqrt(txAntennas * rxAntennas) / sqrt(power);


end

