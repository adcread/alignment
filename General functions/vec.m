function [ vec ] = vec( matrix )
%VEC Summary of this function goes here
%   Detailed explanation goes here

[row, col] = size(matrix);

vec = matrix(:,1);

for i = 2:col
    vec = [vec ; matrix(:,i)];
end



end

