% calculates the correlation matrix for a Monte Carlo simulation of the
% chennel generating functions used in blockalignment.m

txAntennas = 3;

rxAntennas = 2;

a = 0.0;        % correlation coefficient between Tx/Rx antennas

for iteration = 1:10000

    R_t = eye(txAntennas);
    R_r = eye(rxAntennas);
    for g = 1:txAntennas
        R_t = R_t + diag(((a^g)*ones(txAntennas-g,1)),g) + diag(((a^g)*ones(txAntennas-g,1)),-g); % generates diagionals in the antennna covariance matrix 
    end                                                                                                 % - this should be updated with the Bessel function
     for g = 1:rxAntennas                                                                          % approach described in http://www.ece.nus.edu.sg/stfpage/elehht/Teaching/EE6832/Lecture%20Notes%5CMultiple%20Antennas%20for%20MIMO%20Communications%20-%20Channel%20Correlation.pdf
        R_r = R_r + diag(((a^g)*ones(rxAntennas-g,1)),g) + diag(((a^g)*ones(rxAntennas-g,1)),-g);
    end           
    H{iteration} = KroneckerChannel(txAntennas,rxAntennas,R_t,R_r);
end

H_mean = zeros(rxAntennas);

for iteration = 1:10000
    H_mean = H_mean + (H{iteration}*H{iteration}');
end

H_mean = H_mean / 10000


