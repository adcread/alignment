% Estimate the rank using the magnitude of Q's diagonal for M = 1:10, N =
% 1:10

noTx = 5;
noRx = 5;
noIter = 50;

B = cell(noTx,noRx,noIter);

estimate = [];
est_mean = [];
temp = [];
maxCrosscorr = [];

for txAntennas = 1:noTx
    for rxAntennas = 1:noRx
        for iter = 1:noIter
            string = ['Tx = ' num2str(txAntennas) ' Rx = ' num2str(rxAntennas) ' Iteration = ',num2str(iter)];
            display(string);
            KroneckerEstimation;
            temp = ensembleAve(transmittedSequenceCrosscorrelation);
            if txAntennas >1;
                for i = 1:txAntennas
                    temp(i,i,:) = 0;
                end
            end
            B{txAntennas,rxAntennas,iter} = b;
            maxCrosscorr(txAntennas,rxAntennas,iter) = max(abs(temp(blockLength+1:end)))/txAntennas;
            %delta(txAntennas,rxAntennas,iter) = norm((Z - kron(b,rxCorrelation)),'fro');
            estimate(txAntennas,rxAntennas,iter) = mean(abs(diag(b)));ense
        end
    end
end

for i = 1:txAntennas
    for j = 1:rxAntennas
        est_mean(i,j) = mean(squeeze(estimate(i,j,:)));
    end
end

for t = 1:noTx
    for r = 1:noRx
        temp = [];
        for i = 1:noIter
            temp = [temp; abs(diag(B{t,r,iter}))];
        end
        d{t,r} = temp;
    end
end

figure();

for drawRow = 1:noTx
    for drawCol = 1:noRx
        plotLocation = (drawRow-1) * noTx + drawCol;
        subplot(noTx,noRx,plotLocation);
        hist(d{drawRow,drawCol});
        title(['Tx ' num2str(drawRow) ', Rx' num2str(drawCol)]);
    end
end


noColours = 2 *max([rxAntennas txAntennas]);

figure;
subplot(2,1,1);
imagesc(est_mean);
title('Rank estimate');
xlabel('RX antennas');
ylabel('TX antennas');
colormap(parula(noColours));
colorbar;
subplot(2,1,2);
imagesc(mean(maxCrosscorr,3));
title('Maximum Crosscorrelation exhibited (normalised)');
xlabel('RX antennas');
ylabel('TX antennas');
colormap(parula(noColours));
colorbar;

