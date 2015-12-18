%% Investigate the effect of parameter 'sequenceLength' on rank Estimate.

sequenceLengthParameter = 2.^(1:8);

txAntennas = 5;

rxAntennas = 5;

noIter = 100;

estimate = [];
est_mean = [];
temp = [];
maxCrosscorr = [];
B = cell(length(sequenceLengthParameter),noIter);

for iteration = 1:length(sequenceLengthParameter)
    
    sequenceLength = sequenceLengthParameter(iteration);
    

    
    for iter = 1:noIter
        display([num2str(sequenceLengthParameter(iteration)) ' ' num2str(iter)]);
        KroneckerEstimation;
        B{iteration,iter} = b;
        estimate(iteration,iter) = mean(abs(diag(b)));
    end
    
end
