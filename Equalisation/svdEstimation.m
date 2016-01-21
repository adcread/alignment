% Perform SVD Estimation for rank-excessive channel

addpath('/Users/chris/PhD/alignment/Channel');
addpath('/Users/chris/PhD/alignment/');

H = generateChannel(2,[3 2],[3 2],'kronecker');

[U, sigma, V] = eigenchannel(H);
    
S = circSymAWGN(3,64,1);

Y = H{1,2}*S;

Z = Y * V{1,2}'

