% Measure efficacy of Kronecker Estimation for varying numbers of Tx/Rx
% antennas

noTx = 4;
noRx = 4;
sequenceLength = 127;

diagonal = zeros(noTx,noRx,sequenceLength);
estimates = zeros(noTx,noRx);

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        tic;
        disp(['Performing estimation Tx = ' num2str(txAntennas) ', Rx = ' num2str(rxAntennas) '.']);
        for iter = 1:100
            disp(['Performing Iteration ' num2str(iter) '.']);
%         couplingMatrix = zeros(rxAntennas,txAntennas);
%         couplingMatrix(1,:) = randn(1,txAntennas);
%         couplingMatrix(:,1) = randn(1,rxAntennas);
       
        
        % Perform estimation
        KroneckerEstimation;
        
        % Store values of Q's diagonal in matrix for later analysis
        B(txAntennas,rxAntennas,iter,:,:) = b;
        
        % Calculate estimate of txAntennas from this diagonal
        %estimates(txAntennas,rxAntennas) = mean(abs(diag(b)));
        end
        toc;
        
    end
    
end

%% Present the results in two figures
figure();

imagesc(estimates);
colorbar;

figure();

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        
        subplot(noTx,noRx,(noTx * (txAntennas-1) + rxAntennas));
        
        histogram(abs(diag(squeeze(B(txAntennas,rxAntennas,:,:)))),20);
        title(['Tx = ' num2str(txAntennas) ', Rx = ' num2str(rxAntennas) ', Est = ' num2str(estimates(txAntennas,rxAntennas))]);
    
    end

end

    
        