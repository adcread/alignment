% takes an Nx1 signal 'x' and calculates the mean over multiple symbols.

for symbol = 2:length(signal)
    Q{symbol} = signal(:,symbol) * signal(:,symbol)';
    G{symbol} = zeros(length(Q{symbol}));
    for i = 1:length(Q{symbol})
        for j = 1:length(Q{symbol})
            G{symbol}(i,j) = G{symbol-1}(i,j) + Q{symbol}(i,j);
        end
    end
    T(symbol) = trace(G{symbol});
    T(symbol) = T(symbol)/symbol;
end
plot(T);

