% 
% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   1.0         19/11/14            Initial version


function [ U, H_eig, V ] = eigenchannel( H )
%EIGENCHANNEL Uses Single Value Decomposition (SVD) to create an eigenchannel from a channel matrix H 
%   

    [K, N, M] = cellDimensions(H);

    U = cell(K,K);
    V = cell(K,K);
    H_eig = cell(K,K);

    for i = 1:K
        for j = 1:K
            [U{i,j}, H_eig{i,j}, V{i,j}] = svd(H{i,j});
        end
    end
end

