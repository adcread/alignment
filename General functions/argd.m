function [ arg ] = argd( cart )
%ARGD Summary of this function goes here
%   Detailed explanation goes here

argr = cart2pol(real(cart),imag(cart));
arg = argr/pi * 180;

if (arg<0)
    arg = arg + 360;
end

end

