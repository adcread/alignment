
DoF = 0.5;
signalPower = 40;
noisePower = 1;

constellation = generateCodebook(DoF, signalPower, noisePower);

scatterplot(constellation);
hold on;

plotCircle(0,0,sqrt(2*signalPower));

for i = 1:length(constellation)
    plotCircle(real(constellation(i)),imag(constellation(i)),sqrt(2*noisePower));
end

for i = 1:7
    for j = 1:7
    distance(i,j) = abs(constellation(i)-constellation(j));
    end
end