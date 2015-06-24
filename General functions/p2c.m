function [ cart ] = p2c( theta, rho )
%P2C Summary of this function goes here
%   Detailed explanation goes here

cart = rho * cos((theta/360) * 2*pi) + 1i * rho* sin((theta/360) * 2*pi);

end

