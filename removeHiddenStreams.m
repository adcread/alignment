function [ commonSubspace ] = removeHiddenStreams( signal, signalSpace, numberOfStreams )
%REMOVEHIDDENSTREAMS removes the 'hidden' streams from a signal subspace.
%   signal = the signal to be modified
%   signalSpace = the subspace spanned by the received signal
%   numberOfDimensions = the number of stream to be removed from the signal.
   
    signalSubspace = cell(1,numberOfStreams);
    signalSpaceBasis = cell(1,numberOfStreams);

    for stream = 1:numberOfStreams
        signalSubspace{stream} = signalSpace;                                  % signal space is an orthogonal base already due to unitary precoding matrix
        signalSubspace{stream} = circshift(signalSubspace{stream},numberOfStreams,2);     % select the basis function to be used in G-S Orthogonalisation
        signalSpaceBasis{stream} = orthBasis(signalSubspace{stream});
    end
    
    commonSubspace = signal;
    
    for stream = 1:numberOfStreams                                         % successively remove each basis function identified as 'hidden'.
        commonSubspace = convertBases(commonSubspace,signalSpaceBasis{stream});
        commonSubspace(1) = 0;
        commonSubspace = signalSpaceBasis{stream} * commonSubspace;       
    end

end

