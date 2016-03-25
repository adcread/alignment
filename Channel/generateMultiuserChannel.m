function [ H ] = generateChannel( users, txAntennas, rxAntennas, channelModel, modelArg1, ModelArg2 )
%GENERATECHANNEL Summary of this function goes here
%   Detailed explanation goes here

H = cell(users,users);

if strcmp(channelModel,'scalar')

elseif strcmp(channelModel, 'czink')

elseif strcmp(channelModel,'complex')

elseif strcmp(channelModel,'gaussian')

    for i = 1:users
        for j = 1:users
            H{i,j} = normrnd(1,1,txAntennas(j),rxAntennas(i));
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
    % W = sqrtm(R); 
    W = chol(R);

    % Create Kronecker correlated channel H

    w = W * h;
    H = reshape(w, N, M);
    
elseif strcmp(channelModel,'weichselberger')
    
    % take in the Coupling matrix as argument
    
    couplingMatrix = modelArg1;
    
    % Create sets of eigenbases for the transmit and receive antennas
    
    [txEigenbase, ~] = eigs(randn(rxAntennas, txAntennas) + 1j * randn(rxAntennas,txAntennas));
    [rxEigenbase, ~] = eigs(randn(rxAntennas,txAntennas) + 1j * randn(rxAntennas,txAntennas));
    
    H = txEigenbase * (couplingMatrix .* randn(rxAntennas,txAntennas) + 1j * randn(rxAntennas,txAntennas)) * rxEigenbase';
    
else
    % no operation, return empty cell of channels
    
end

for i = 1:users
    for j = 1:users
        power = trace(H{i,j} * H{i,j}');
        H{i,j} = H{i,j} * sqrt(txAntennas(i) * rxAntennas(j)) / sqrt(power);
    end
end

end

