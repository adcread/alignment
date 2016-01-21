% Perform iterations of the block alignment, capturing results of each
% simulation.

numberOfIterations = 2048;

errorRate = cell(2,1);
for i = 1:2
    errorRate{i} = zeros(2,numberOfIterations);
    conditionStat = zeros(3,numberOfIterations);
end

for iteration = 1:numberOfIterations
    display(['Iteration: ' num2str(iteration)]);
    blockalignment;
    conditionStat(1,iteration) = condition(1,1);
    conditionStat(2,iteration) = condition(2,2);
    conditionStat(3,iteration) = condition(2,1);
    errorRate{1}(:,iteration) = [publicSER(1) privateSER(1)].';
    errorRate{2}(:,iteration) = [publicSER(2) privateSER(2)].';
end