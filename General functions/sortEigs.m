function [ sortedEigenvectors, sortedEigenvalues ] = sortEigs( inputArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[unsortedEigenvectors, unsortedEigenvalues] = eigs(inputArray);

eigenvalues = diag(unsortedEigenvalues);

rank = length(eigenvalues);

i = 1;

while i <= rank
    
    [value, location] = max(eigenvalues);
    
    if length(location) < 2
        
        sortedEigenvalues(i,i) = value;
        sortedEigenvectors(:,i) = unsortedEigenvectors(:,location);
        eigenvalues(location) = 0;
                
    else

        % pray this doesn't happen
        error('input matrix has two equal eigenchannels');
        
    end
    
    i = i+1;

end

end

