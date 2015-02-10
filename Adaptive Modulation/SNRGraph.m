SNR = 0:0.1:100;
scalarSNR = 10.^(SNR/10);

awgnChannelCapacity = log2(1 + 10.^(SNR/10));

squareConstellationSizes = 2:2:100;
crossConstellationSizes = 3:2:101;

squareConstellationSNRs = 10*log10((((2.^squareConstellationSizes)-1)*23.4423)/3);
crossConstellationSNRs = 10*log10((((2.^crossConstellationSizes)-1)*23.4423)/3);

squareQAMChannelCapacity = zeros(1,length(SNR));

for snr = 1:length(SNR)
    squareQAMChannelCapacity(snr) = findLargestConstellation(SNR(snr),squareConstellationSNRs,squareConstellationSizes);
    crossQAMChannelCapacity(snr) = findLargestConstellation(SNR(snr),crossConstellationSNRs,crossConstellationSizes);
end

plot(SNR,awgnChannelCapacity);
hold on;
plot(SNR,squareQAMChannelCapacity,'r');
plot(SNR,crossQAMChannelCapacity,'g');
ylabel('Channel capacity (bits/sec/Hz)');
xlabel('SNR (dB)');
title('Channel Capacity');
legend('Gaussian Codebook','Adaptive QAM (Square Constellation)','Adaptive QAM (Cross Constellation)')

hold off
figure;
plot(SNR,(awgnChannelCapacity-(max(squareQAMChannelCapacity,crossQAMChannelCapacity))));
ylabel('QAM Rate Penalty (bits/sec/Hz)');
xlabel('SNR (dB)');
title('Rate penalty from using Adaptive QAM (\Gamma)')
