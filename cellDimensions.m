function [ K, N, M ] = cellDimensions( celldim )
%MATRIXDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here

    K = size(celldim,1);                                                     % Get the number of users from the size of the cell matrix H
    matrix = celldim{1,1};                                                
    N = size(matrix,1);
    M = size(matrix,2);
    
end

