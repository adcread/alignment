function [ rho ] = arrayCorrelation( noAntennas, distance )
%ARRAYCORRELATION Creates correlation matrix of Uniform Linear Array (ULA)
%of antennas given the distance between them
%   Detailed explanation goes here

for i = 1:noAntennas
    for j = 1:noAntennas
        rho(i,j)=besselj(0,2*pi*distance*abs(j-i));
    end
end

end

