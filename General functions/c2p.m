function [ theta, rho ] = c2p( cart )
%C2P Returns the polar form of a complex number
%   Detailed explanation goes here

[theta, rho] = cart2pol(real(cart),imag(cart));

end

