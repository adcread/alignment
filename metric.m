addpath('C:\PhD\alignment\Channel');
addpath('C:\PhD\alignment\General functions');

H = generateChannel(2,[3 2],[3 2],'kronecker');

[U, sigma, V] = eigenchannel(H);

R_d = H{1,1}*H{1,1}';

R_i = H{2,1}*H{2,1}';

[Q, Lambda] = sortEigs(R_d);
[P, Sigma] = sortEigs(R_i);

for i = 1:3
    for j = 1:3
        direct_dp(i,j) = dot(Q(:,i),Q(:,j));
        int_dp(i,j) = dot(P(:,i),P(:,j));
        cross_dp(i,j) = dot(Q(:,i),P(:,j))
    end
end

