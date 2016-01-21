% Measure efficacy of Kronecker Estimation for varying numbers of Tx/Rx
% antennas

noTx = 4;
noRx = 4;
sequenceLength = 256;

diagonal = zeros(noTx,noRx,sequenceLength);
estimates = zeros(noTx,noRx);

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        
        disp(['Performing estimation Tx = ' num2str(txAntennas) ', Rx = ' num2str(rxAntennas) '.']);
        
        % Perform estimation
        KroneckerEstimation;
        
        % Store values of Q's diagonal in matrix for later analysis
        diagonal(txAntennas,rxAntennas,:) = diag(b);
        
        % Calculate estimate of txAntennas from this diagonal
        estimates(txAntennas,rxAntennas) = mean(abs(diag(b)));
        
    end
    
end

% Present the results in two figures
figure();

imagesc(estimates);
colorbar;

figure();

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        
        subplot(noTx,noRx,(noTx * (txAntennas-1) + rxAntennas));
        
        histogram(abs(squeeze(diagonal(txAntennas,rxAntennas,:))),20);
        title(['Tx = ' num2str(txAntennas) ', Rx = ' num2str(rxAntennas) ', Est = ' num2str(estimates(txAntennas,rxAntennas))]);
    
    end

end

    
        