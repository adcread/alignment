% Attempts to perform distributed optimisation between 2 users using Goal
% Programming to achieve optimal DoF in (M1, M2, N1, N2) MIMO IC.

initialisationConstants;

SNR = zeros(users);

for i = 1:users
    SNR(i,:) = sqrt(transmitPower(i)) * fading(i,:);
end

if (symmetricChannel == true)
    rho = directChannelFade;
else
    rho = 2;                                                                % Determine a baseline 'rho' to be used as the base for alpha calculations.
end

alpha = log2(SNR)/log2(rho);
beta = zeros(users);

for i=1:users
    for j = 1:users
        beta(i,j) = max([0 (alpha(i,i)-alpha(i,j))]);
    end
end

d1p = 0;                                                                    % User 1 private DoF
d2p = 0;

d1c1 = 1;                                                                   % User 1 common DoF calculated by User 1
d1c2 = NaN;                                                                   % User 1 common DoF calculated by User 2
d2c1 = NaN;
d2c2 = min([rxAntennas(1), txAntennas(1), rxAntennas(2)]);
d1c = 0;
d2c = 0;

% User 1 determines DoF using the absolute maximum of d1c as side
% information

iteration = 0;

while ((d1c1 >= 0) && ((d1c1 ~=d1c2) || (d2c1 ~=d2c2)))
    [d1p, d1c1, d2c1] = degreesOfFreedom2UserMIMOIC_LGP(1, maxIter, txAntennas, rxAntennas, alpha, beta, d2c2, 2);

    [d2p, d2c2, d1c2] = degreesOfFreedom2UserMIMOIC_LGP(2, maxIter, txAntennas, rxAntennas, alpha, beta, d1c1, 2);
    iteration = iteration + 1;
    d1c = min([d1c1 d1c2]);
    d2c = min([d2c1 d2c2]);
    d1c1 = d1c;
    d2c2 = d2c;
    sumDoF = [d1p d1c d2p d2c];
    string = ['Iteration ',num2str(iteration),': ',num2str(sumDoF)];
    disp(string);
end