function [ K, N, M ] = matrixDimensions( matrix )
%MATRIXDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here

    K = size(matrix,1);                                                     % Get the number of users from the size of the cell matrix H                                             
    N = size(matrix,1);
    M = size(matrix,2);
    
end

