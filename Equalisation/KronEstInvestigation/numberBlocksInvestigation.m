%% Investigate the effect of parameter 'sequenceLength' on rank Estimate.

%addpath('/Users/chris/PhD/alignment/Equalisation');

noBlockParameter = 2.^(1:9);

sequenceLength = 127;
rootSequence = generateGoldCodes(round(log2(sequenceLength)));

txAntennas = 5;

rxAntennas = 5;

SNR = 100;

noIter = 100;

est_mean = [];
temp = [];
maxCrosscorr = [];
B = cell(length(noBlockParameter),noIter);

estimate = zeros(length(noBlockParameter),noIter);

for iterNoBlocks = 1:length(noBlockParameter)
    
    noBlocks = noBlockParameter(iterNoBlocks);
       
    for iterDistribution = 1:noIter
        display(['Number of blocks: ' num2str(noBlockParameter(iterNoBlocks)) ' Iteration: ' num2str(iterDistribution)]);
        KroneckerEstimation;
        B{iterNoBlocks,iterDistribution} = b;
        estimate(iterNoBlocks,iterDistribution) = mean(abs(diag(b)));
    end
    
end
%%
figure;

for imageIndex = 1:length(noBlockParameter)
    means(imageIndex) = mean(estimate(imageIndex,:)/SNR);
    vars(imageIndex) = var(estimate(imageIndex,:)/SNR);
    subplot(3,3,imageIndex);
    %histogram(estimate(imageIndex,:)/SNR);
    histfit(estimate(imageIndex,:)/SNR);
    xlim([means(imageIndex)-6*sqrt(vars(imageIndex)) means(imageIndex)+6*sqrt(vars(imageIndex))]);
    ylim([0 150]);
    hold on;
    vardist = var(estimate(imageIndex,:)/SNR);
    text(means(imageIndex),140,['\mu = ' num2str(means(imageIndex)) ', \sigma^2 = ' num2str(vars(imageIndex))],'HorizontalAlignment','center');
    %plot(txAntennas-3*(sqrt(vardist)):0.05:txAntennas+3*(sqrt(vardist)),normpdf(txAntennas-3*(sqrt(vardist)):0.05:txAntennas+3*(sqrt(vardist)),mean(estimate(imageIndex,:)/SNR),vardist),'r');
    title([num2str(2^imageIndex) ' Sequences']);
end
