% create the opposing channel from one defined such that it meets the
% criteria for 'spatial alignment'

% channel = the MIMO channel to be used as a basis
% targetCorr = target level of correlation (i.e. J_target)

targetCorr = 0;

R_0 = channel * channel';

[M,N] = size(channel);

[u, Lambda] = eigs(R_0);

Gamma = Lambda;

noise = 1;

[Z, dummy1, dummy2] = svd(circSymAWGN(M,N,1));

normCorr = J(R_0,(Z*Gamma*Z'),noise);

