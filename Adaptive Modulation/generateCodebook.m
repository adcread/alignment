function [ codebook ] = generateCodebook( DoF, signalPower, noisePower )
%GENERATECODEBOOK Creates a complex Gaussian codebook of length N
%   Calculate the SNR of the link to be coded for and the desired DoF to be
%   achieved, then create a codebook of that size and return it.

SNR = signalPower / noisePower;
SNRdB = DoF * (10*log10(SNR));

% Create the table to perform SNR lookups from in dB
SNRtable = 0:0.1:100;

% Calculate the maximum constellation sizes for square and cross
% constellations (diamond constellations can be added as a separate line at
% a later date if required)

squareConstellationSizes = 2:2:100;
crossConstellationSizes = 1:2:99;

squareConstellationSNRs = 10*log10((((2.^squareConstellationSizes)-1)*23.4423)/3);
crossConstellationSNRs = 10*log10((((2.^crossConstellationSizes)-1)*23.4423)/3);

squareQAMChannelCapacity = zeros(1,length(SNR));

for snr = 1:length(SNRtable)
    squareQAMChannelCapacity(snr) = findLargestConstellation(SNRtable(snr),squareConstellationSNRs,squareConstellationSizes);
    crossQAMChannelCapacity(snr) = findLargestConstellation(SNRtable(snr),crossConstellationSNRs,crossConstellationSizes);
end

% Combine the channel capacity curves to find supremum
adaptiveQAMChannelCapacity = max(squareQAMChannelCapacity,crossQAMChannelCapacity);

for n = 1:length(SNRtable)
    if ((SNRtable(n) < SNRdB && SNRtable(n+1) > SNRdB) || (SNRtable(n) == SNRdB))
        display(['Creating QAM constellation for SNR = ' num2str(SNRtable(n)) ' dB.'])
        break
    end
end

M = 2^adaptiveQAMChannelCapacity(n);
if M >=4
    alphabet = 0:(M-1);

    codebook = qammod(alphabet,M);

    display(['Created QAM constellation M = ' num2str(M) '.']);

%    scatterplot(codebook);
else
    display(['Unable to create codebook: SNR too low']);
    codebook = 0;
end

