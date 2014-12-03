function [ U, sigma, V ] = eigenchannels( H, txAntennas, rxAntennas )
%EIGENCHANNELS Summary of this function goes here
%   Detailed explanation goes here

    [K, N, M] = cellDimensions(H);

    U = cell(K,K);
    sigma = cell(K,K);
    V = cell(K,K);

    for i = 1:K
        for j = 1:K
            u_int = zeros(txAntennas(i));
            sigma_int = zeros(rxAntennas(i),txAntennas(j));
            v_int = zeros(rxAntennas(j));
            [u_int, sigma_int, v_int] = eigenchannel(H{i,j});
            U{i,j} = u_int;
            sigma{i,j} = sigma_int;
            V{i,j} = v_int;
        end
    end

end

