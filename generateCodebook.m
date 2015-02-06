function [ codebook ] = generateCodebook( DoF, signalPower, noisePower )
%GENERATECODEBOOK Creates a complex Gaussian codebook of length N
%   Detailed explanation goes here

numberOfCodewords = floor(DoF * (1 + (signalPower/noisePower)));

codebook = zeros(numberOfCodewords, 1);

originalConstellation = gaussianCodeword(200*numberOfCodewords, signalPower);

noiseSphere = sqrt(2*noisePower);

decodeable = false;

while decodeable == false
    decodeable = true;
    originalConstellation = originalConstellation(abs(originalConstellation)<=sqrt(2*signalPower));
    codebook = randsample(originalConstellation,numberOfCodewords);
    for i =1:numberOfCodewords
        for j = 1:numberOfCodewords
            if i ~=j
                if abs(codebook(i)-codebook(j)) <= 2*noiseSphere
                    decodeable = false;
                end
            end
        end
    end     
end



end

