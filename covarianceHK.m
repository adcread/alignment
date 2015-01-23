function [ covariancePri, covariancePub ] = covarianceHK( SNR, H )
%COVARIANCEHK Calculates the covariance matrices of the public and private messages for a simple Han-Kobayashi scheme. 
%   Detailed explanation goes here

    [N, M] = size(H);

    covariancePri = (1/M) * inv(eye(M) + SNR * (ctranspose(H) * H));
    covariancePub = (eye(M)/M) - covariancePri;

    
end

