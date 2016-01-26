function [ H ] = KroneckerChannel( M, N, R_tx, R_rx )
% Frist example of the Kronecker model as used in section 4.2.2 where it uses the 
% full Kronecker matrix.
% Tim Brown, January 2012

R = kron(R_tx,R_rx); % Apply Kronecker product to resolve the full matrix
    
%Generate a long vector form of NTx by NRx channels with N samples (note
%the book only uses one sample


h = (randn(N * M,1) + 1j * randn(N * M,1))/sqrt(2);

% Either the sqrtm or Choelsky factorisation can be used here.
% W = sqrtm(R); 
W = chol(R);

% Create Kronecker correlated channel H

w = W * h;
H = reshape(w, N, M);

end