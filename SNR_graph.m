
% Plots graphs for each 

P_e = 10.^(0:-1:-9);
M = 2.^(2:2:10);

for i = 1:length(P_e)

    for j = 1:length(M)
        d = sqrt(12/(M(j)-1));
        SNR(i,j) = qfuncinv(P_e(i)/(4*(1-1/sqrt(M(j)))))*d*(M(j)-1)/6;
    end
end

figure;
hold on;
