%% Investigate the effect of parameter 'sequenceLength' on rank Estimate.

addpath('/Users/chris/PhD/alignment/Equalisation');

blockNumberParameter = 2.^(1:10);

sequenceLength = 64;

txAntennas = 5;

rxAntennas = 5;

noIter = 100;

est_mean = [];
temp = [];
maxCrosscorr = [];
B = cell(length(blockNumberParameter),noIter);

estimate = zeros(length(blockNumberParameter),noIter);

for iteration = 1:length(blockNumberParameter)
    
    noBlocks = blockNumberParameter(iteration);
       
    for iter = 1:noIter
        display([num2str(blockNumberParameter(iteration)) ' ' num2str(iter)]);
        KroneckerEstimation;
        B{iteration,iter} = b;
        estimate(iteration,iter) = mean(abs(diag(b)));
    end
    
end

figure;

for imageIndex = 1:10
    subplot(4,4,imageIndex);
    histogram(estimate(imageIndex,:));
    title([num2str(2^imageIndex) ' Sequences']);
end
