%% Investigate the effect of parameter 'sequenceLength' on rank Estimate.

%addpath('/Users/chris/PhD/alignment/Equalisation');


noIter = 1000;

txAntennas = 5;
rxAntennas = 5;
sequenceLength = 127;
SNR = 10;
noBlocks = 64;

rootSequence = generateGoldCodes(round(log2(sequenceLength)));

est_mean = [];
temp = [];
maxCrosscorr = [];

B = cell(1,noIter);
Q_dist = cell(1,noIter);
     
bb = [];
qq = [];

    for iterDistribution = 1:noIter
        display(['Iteration: ' num2str(iterDistribution)]);
        KroneckerEstimation;
        B{iterDistribution} = b;
        Q_dist{iterDistribution} =ensembleAve(Q);
        bb = [bb; diag(B{iterDistribution})];
        qq = [qq; diag(Q_dist{iterDistribution})];
        end
    
