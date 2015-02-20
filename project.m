function [ projection, orthogonalComponent ] = project( vector, subspace )
%PROJECT Projects the vector given onto the subspace described by the
%subspace matrix
%   Detailed explanation goes here

    projection = (cpxIP(vector,subspace) / cpxIP(subspace,subspace)) * subspace;
    orthogonalComponent = vector - projection;
    
end

