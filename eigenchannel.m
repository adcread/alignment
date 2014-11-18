% 
% Version 0.1 Chris Waters

% Change History
%
%   Version     Date                Comments
%   1.0         18/11/14            Initial version


function [ U, H_eig, V ] = eigenchannel( H )
%EIGENCHANNEL Uses Single Value Decomposition (SVD) to create an eigenchannel from a channel matrix H 
%   

[U, H_eig, V] = svd(H);

end

