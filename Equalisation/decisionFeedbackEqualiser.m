% create a vector 'signal' in preparation for running this script
% containing the equalised signal to be decoded, and a vector 'sequence'
% that contains the training sequence to be used to control the filter.
% Also pass through a vector 'codebook' containin the codebook used.

N = min([length(signal) length(sequence)]);

delta = 0.9;

for symbol = 1:N
    expectedSignal(symbol) = codebook(sequence(symbol));
end

% DFE solutions

% DFECoefficient(1) = 1;
% recoveredSignal(1) = signal(1)/DFECoefficient(1);
% recoveredSignal(2) = signal(1)/DFECoefficient(1);
% 
% for symbol = 3:trainingSymbols
%     DFECoefficient(symbol) = (signal(symbol) / expectedSignal(symbol)) - (delta(1) * DFECoefficient(symbol-1) + delta(2) * DFECoefficient(symbol-2));
%     recoveredSignal(symbol) = signal(symbol) / DFECoefficient(symbol);
% end
% 
% for symbol = trainingSymbols+1:N
%     recoveredSignal(symbol) = signal(symbol) / DFECoefficient(trainingSymbols);
% end

% LMS Filter solution

DFECoefficient(1) = 1;
windowLength = 8;

windows = trainingSymbols/windowLength;

for symbol = 1:N
    recoveredSignal(symbol) = signal(symbol) * DFECoefficient(symbol);
    recoveredSignal(symbol)
    err(symbol) = expectedSignal(symbol) / recoveredSignal(symbol);
    err(symbol)
    step(symbol) = delta * err(symbol);
    step(symbol)
    DFECoefficient(symbol+1) = DFECoefficient(symbol) * step(symbol);
    DFECoefficient(symbol+1)
end

% for symbol = trainingSymbols:N
%     recoveredSignal(symbol) = signal(symbol) * DFECoefficient(trainingSymbols);
% end

scatter(real(recoveredSignal),imag(recoveredSignal),'.');
hold on;
scatter(real(codebook),imag(codebook),'o','g');
scatter(real(expectedSignal),imag(expectedSignal),'o','r');
hold off;
title('Recovered signal through DFE')