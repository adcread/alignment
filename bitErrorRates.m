% Monte Carlo simulation of n blocks with new channel matrices for each -
% measure BER for each stream per iteration.

N = 1000;

displayConstellations = false;

publicSERs = zeros(2,N);
privateSERs = zeros(2,N);

for iteration = 1:N
    
    blockalignment;
    
    publicSERs(:,iteration) = publicSER;
    privateSERs(:,iteration) = privateSER;
    cond(iteration) = max(max(condition));
    eigenspace(iteration,:) = eigenspaceCorrelation;
end

for i = 1:2
    publicAveSER(i) = mean(publicSERs(i,:));
    privateAveSER(i) = mean(privateSERs(i,:));
end
