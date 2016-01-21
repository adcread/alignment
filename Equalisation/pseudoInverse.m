function [ piA ] = pseudoInverse( A )
%PSEUDOINVERSE Summary of this function goes here
%   Detailed explanation goes here

piA = inv(A'*A)*A';

xcorr(A);

end

