function [ arg ] = argr( cart )
%ARG Returns the argument of a complex number
%   Detailed explanation goes here

arg = cart2pol(real(cart),imag(cart));

if (arg<0)
    arg = arg + (pi/2);
end

end

