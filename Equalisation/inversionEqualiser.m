
% Performs the channel inversion equalisation from [66] Weikert

addpath('/Users/chris/PhD/alignment/Channel/');
addpath('/Users/chris/PhD/alignment/General functions/');

H = generateChannel(1,3,2,'kronecker');
H = H{1,1}

SNR = 500;

sequenceLength = 64;

sequence = circSymAWGN(3,sequenceLength,1);

for i = 1:2
    
    receivedSignal(i,:) =  sqrt(SNR) * H(i,:) * sequence;% + circSymAWGN(1,sequenceLength,1);
        
end

for i = 1:2

    H_est(i,:) =  1/sqrt(SNR) * receivedSignal(i,:) * (pinv(sequence(1:2,:)'*sequence(1:2,:)) * sequence(1:2,:)');
   
end

H_est

err = H - H_est