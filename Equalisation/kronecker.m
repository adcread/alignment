% Create a training signal, pass it through a Kronecker channel and see
% what it looks like.

addpath('/Users/chris/PhD/alignment/Channel');
addpath('/Users/chris/PhD/alignment/General functions');
addpath('/Users/chris/PhD/alignment/Equalisation');

sequenceLength = 64;
repetitions = 16;

s = circSymAWGN(sequenceLength,4,1);

S = s;

for rep = 2:repetitions
    S = [S ; s];
end

R_r = [1 0.3 0.09 ; 0.3 1 0.3 ; 0.09 0.3 1];
R_t = [1 0.4 0.16 0.064 ; 0.4 1.0 0.4 0.16 ; 0.16 0.4 1 .4 ; 0.064 0.14 .4 1];

H = KroneckerChannel(4,3,R_t, R_r);

Y = H * S.';

y = vec(Y);

G = y * y';

S_hat = S * sqrtm(R_t);

R = zeros(1,(sequenceLength * repetitions));
R_hat = zeros(1,(sequenceLength * repetitions));

for k = 1:rank(R_r)
    for i = 1:(sequenceLength * repetitions)
        R_hat(i) = R_hat(i) + mean(S_hat(:,k).*circshift(S_hat(:,k),i));
    end
end

for k = 1:rank(R_r)
    for i = 1:(sequenceLength * repetitions)
        R(i) = R(i) + mean(S(:,k).*circshift(S(:,k),i));
    end
end

for i = 1:(sequenceLength * repetitions)
    Q(i,:) = circshift(R,i-1,2);
end

% Partition Q into N * N chunks and reverse the kronecker product

Qrank = size(Q,1);
Grank = size(G,1);
Rrank = size(G,1)/size(Q,1);

for i = 1:Rrank
    for j = 1:Rrank
        R_calc(i,j) = G(1+(i-1)*Qrank:i*Qrank,1+(j-1)*Qrank:j*Qrank) ./ Q;
    end
end


    


