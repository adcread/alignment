% Measure efficacy of Kronecker Estimation for varying numbers of Tx/Rx
% antennas

noTx = 9;
noRx = 9;
sequenceLength = 127;
noBlocks = 64;
noIter = 500;

diagonal = zeros(noTx,noRx,sequenceLength);
estimates = zeros(noTx,noRx);

AntennasB = cell(noTx,noRx,noIter);

rootSequence = generateGoldCodes(round(log2(sequenceLength)));

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        tic;
        disp(['Performing estimation Tx = ' num2str(txAntennas) ', Rx = ' num2str(rxAntennas) '.']);
        for iter = 1:noIter
            disp(['Performing Iteration ' num2str(iter) '.']);
       
            % Perform estimation
            KroneckerEstimation;
        
            % Store values of Q's diagonal in matrix for later analysis
            AntennasB{txAntennas,rxAntennas,iter} = b;
        end
        
        % Calculate estimate of txAntennas from this diagonal
        a = [];
        a = cell(1,noIter);
        for iter = 1:noIter
            a{iter} = AntennasB{txAntennas,rxAntennas,iter};
        end
        
        estimates(txAntennas,rxAntennas) = mean(abs(diag(ensembleAve(a))));
        
        toc;
    end
    
end

%% Present the results in two figures
figure();

imagesc(estimates);
colorbar;

figure();

errors = zeros(noTx, noRx, noIter);

for txAntennas = 1:noTx
    
    for rxAntennas = 1:noRx
        
        for iter = 1:noIter
            errors(txAntennas,rxAntennas,iter) = (mean(diag(AntennasB{txAntennas,rxAntennas,iter}))/100000 - txAntennas)^2;
        end
        
    MSE(txAntennas,rxAntennas) = mean(errors(txAntennas,rxAntennas,:));
    
    end

end

for rxAntennas = 1:noRx
    semilogy(real(MSE(:,rxAntennas)));
    hold on;
end
    
        