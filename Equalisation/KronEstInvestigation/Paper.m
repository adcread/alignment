
statisticsFlag = true;

txAntennas = 5;
rxAntennas = 5;

noBlocksParameter = [256 512 1024];
sequenceLengthParameter = [31 63];
snrParameter = [1 2 5 10 15 20 25 30];

noIterations = 1;

B = cell(length(noBlocksParameter),length(sequenceLengthParameter),length(snrParameter));
estimate = cell(length(noBlocksParameter),length(sequenceLengthParameter),length(snrParameter));
MSE = zeros(length(noBlocksParameter),length(sequenceLengthParameter),length(snrParameter));

for noBlocksIndex = 1:length(noBlocksParameter)

    noBlocks = noBlocksParameter(noBlocksIndex);

    for sequenceLengthIndex = 1:length(sequenceLengthParameter)

        sequenceLength = sequenceLengthParameter(sequenceLengthIndex);
        rootSequence = generateGoldCodes(round(log2(sequenceLength)));
            
        for snrIndex = 1:length(snrParameter)
    
            display(['Number of Blocks: ' num2str(noBlocksParameter(noBlocksIndex)) ' Sequence Length: ' num2str(sequenceLengthParameter(sequenceLengthIndex)) ' SNR = ' num2str(snrParameter(snrIndex)) 'dB.']);

            SNR = 10^(snrParameter(snrIndex)/20);
            
            estimate{noBlocksIndex,sequenceLengthIndex,snrIndex} = zeros(1,noIterations);
            
            for iterationIndex = 1:noIterations
                            
                KroneckerEstimation;
                
                B{noBlocksIndex,sequenceLengthIndex,snrIndex,iterationIndex} = b;
                
            end
            
            
        end
    end
end

%%
for noBlocksIndex = 1:length(noBlocksParameter)
    for sequenceLengthIndex = 1:length(sequenceLengthParameter)
        for snrIndex = 1:length(snrParameter)
            for iterationIndex = 1:noIterations
                %estimate{noBlocksIndex,sequenceLengthIndex,snrIndex}(iterationIndex) = (txAntennas/(10^(0.2*snrParameter(snrIndex)))) * mean(abs(diag(B{noBlocksIndex,sequenceLengthIndex,snrIndex,iterationIndex})));
                estimate{noBlocksIndex,sequenceLengthIndex,snrIndex}(iterationIndex) = mean(abs(diag(B{noBlocksIndex,sequenceLengthIndex,snrIndex,iterationIndex})));
            end
            MSE(noBlocksIndex,sequenceLengthIndex,snrIndex) = mean((txAntennas - ((txAntennas^2/(10^(snrParameter(snrIndex)/20))) * (estimate{noBlocksIndex,sequenceLengthIndex,snrIndex} -1))).^2);
        end
    end
end

%%

figure;
hold on;

legendArray = [];

for noBlocksIndex = 1:length(noBlocksParameter)
    for sequenceLengthIndex = 1:length(sequenceLengthParameter)
        
        plot(snrParameter,squeeze(MSE(noBlocksIndex,sequenceLengthIndex,:)));
        legendArray = [legendArray strcat(num2str(noBlocksParameter(noBlocksIndex)),'Blocks, Sequence Length = ',num2str(sequenceLengthParameter(sequenceLengthIndex)))];
 
    end
end

% legend(legendArray);
        