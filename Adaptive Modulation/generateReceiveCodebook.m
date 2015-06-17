function [ codebook ] = generateReceiveCodebook( M )
%GENERATERECEIVECODEBOOK Summary of this function goes here
%   Detailed explanation goes here

alphabet = 0:2^M-1;

codebook = qammod(alphabet, 2^M);

energy = mean((real(codebook).^2)+(imag(codebook).^2));
codebook = codebook * (1/sqrt(energy));

end

