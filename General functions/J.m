function [ J_norm ] = J( R0, R1, sigma )
%J Summary of this function goes here
%   Detailed explanation goes here

d = length(R0);

[u, lambda] = eigs(R0);
[v, gamma] = eigs(R1);

J_abs = log2(det(eye(d) + R0 * pinv(R1 + sigma * eye(d))));

J_min = log2(det(eye(d) + R0 * pinv((u*gamma*u') + sigma * eye(d))));

J_max = log2(det(eye(d) + R0 * pinv((fliplr(u) * gamma * fliplr(u)') + sigma * eye(d))));

J_norm = (J_abs - J_min) / (J_max - J_min);

end

