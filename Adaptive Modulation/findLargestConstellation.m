function [ largestConstellation ] = findLargestConstellation( SNR, tableOfSNRs, tableOfRates )
%FINDLARGESTCONSTELLATION Summary of this function goes here
%   Detailed explanation goes here

    largestConstellation = tableOfRates(find(tableOfSNRs<=SNR,1,'last'));
    if isempty(largestConstellation)
        largestConstellation = 0;
    end

end

