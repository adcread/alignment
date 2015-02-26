SNR = 0:0.1:50;
scalarSNR = 10.^(SNR/10);

awgnChannelCapacity = log2(1 + 10.^(SNR/10));



squareQAMChannelCapacity = floor(log2((1 + (3*scalarSNR)/23.4423)));
crossQAMChannelCapacity = ((log2(((32/31)* (3*scalarSNR + 1)/23.4423))));

plot(SNR,awgnChannelCapacity);
hold on;
plot(SNR,squareQAMChannelCapacity,'r');
plot(SNR,crossQAMChannelCapacity,'g');
ylabel('Channel capacity (bits/sec/Hz)');
xlabel('SNR (dB)');
title('Channel Capacity');
legend('Gaussian Codebook','Adaptive QAM (Square Constellation)','Adaptive QAM (Cross Constellation)')
