
addpath('C:\PhD\alignment\Channel');

users = 2;
txAntennas = [3 2];
rxAntennas = [3 2];

H = generateChannel(users, txAntennas, rxAntennas, 'gaussian');

    R0 = H{1,1}*H{1,1}';
    R1 = H{2,1}*H{2,1}';
    
    [u, lambda] = eigs(R0);
    [v, gamma] = eigs(R1);

    
correlation = J(R0,R1,1)

minCorrelation = J(R0,(u * gamma * u'), 1);
maxCorrelation = J(R0,(fliplr(u)*gamma*fliplr(u)'),1);

normCorrelation = (correlation - minCorrelation)/(maxCorrelation - minCorrelation)

% J = log2 * det(eye(3) * lamda{i} * u{i}' * v{i} * pinv(eye(3) + gamma{i}) * v{i}' * u{i});
