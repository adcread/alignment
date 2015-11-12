% calculates the waterfilling power division given a channel.

% eigVal = svd(channel);

singularVal = [3.2 1.2 0.44 0.24];

rank = length(singularVal);

noisePower = 1;

waterBins = (eigVal)/noisePower;

totalPower = 1;



