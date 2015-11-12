%% Derive a plot of the distribution of channel conditions for channel generating functions

iterations = 2;    % Monte Carlo iterations

txAntennas = 2;
rxAntennas = 2;

for i = 1:iterations
    
    H = generateChannel(1, txAntennas, rxAntennas, 'gaussian');
    
    [U, Lambda, V] = eigenchannel(H);
    
    condition(i) = max(nonzeros(Lambda{1}))/min(nonzeros(Lambda{1}));
    SERpub(:,i) = publicSER;
    SERpri(:,i) = privateSER;
end

% edges = [min(condition):1:50 max(condition)];
% 
% histogram(condition,edges,'Normalization','pdf');
% xlabel('Condition Number');
% ylabel('Number of channels');
