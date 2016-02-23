SNRdB = 0:30;

txAntennas = 4;

rxAntennas = 4;

sequenceLength = 127;

noIter = 1000;

estimateSNR = zeros(length(SNRdB),noIter);

for step = 1:length(SNRdB)
    
    for iter = 1:noIter
        disp(['SNR= ' num2str(SNRdB(step)) ' dB - Iteration ' num2str(iter)]);
        
        SNR = 10^(SNRdB(step)/10);

        KroneckerEstimation;

        estimateSNR(step,iter) = mean(abs(diag(b)));

        errors(step,iter) = txAntennas - estimateSNR(step,iter) / power(10,(SNRdB(step)/10));
    end
    
end

figure;
plot(SNRdB,abs(mean(errors,2)));
    