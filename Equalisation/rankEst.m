% Estimate the rank using the magnitude of Q's diagonal for M = 1:10, N =
% 1:10

for txAntennas = 1:6
    for rxAntennas = 1:6
        for iter = 1:16
            string = ['Tx = ' num2str(txAntennas) ' Rx = ' num2str(rxAntennas) ' Iteration = ',num2str(iter)];
            display(string);
            KroneckerEstimation;
            temp = ensembleAve(transmittedSequenceCrosscorrelation);
            if txAntennas >1;
                for i = 1:txAntennas
                    temp(i,i,:) = 0;
                end
            end
            maxCrosscorr(txAntennas,rxAntennas,iter) = max(abs(temp(blockLength+1:end)))/txAntennas;
            delta(txAntennas,rxAntennas,iter) = norm((Z - kron(b,rxCorrelation)),'fro');
            estimate(txAntennas,rxAntennas,iter) = mean(abs(diag(b)));
        end
    end
end

for i = 1:txAntennas
    for j = 1:rxAntennas
        est_mean(i,j) = mean(squeeze(estimate(i,j,:)));
    end
end

noColours = 2 *max([rxAntennas txAntennas]);

figure;
subplot(2,1,1);
imagesc(est_mean);
title('Rank estimate');
xlabel('RX antennas');
ylabel('TX antennas');
colormap(jet(noColours));
colorbar;
subplot(2,1,2);
imagesc(mean(maxCrosscorr,3));
title('Maximum Crosscorrelation exhibited (normalised)');
xlabel('RX antennas');
ylabel('TX antennas');
colormap(jet(noColours));
colorbar;

