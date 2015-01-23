function [ H ] = generateChannel( users, txAntennas, rxAntennas, channelModel )
%GENERATECHANNEL Summary of this function goes here
%   Detailed explanation goes here

H = cell(users,users);

if strcmp(channelModel,'scalar')
    

elseif strcmp(channelModel,'complex')


elseif strcmp(channelModel,'kronecker')

    a = 0.5;
    
    for i = 1:users             % transmitting node
        for j = 1:users            % receiving node
            R_t = eye(txAntennas(i));
            R_r = eye(rxAntennas(j));
            for g = 1:txAntennas(i)
                R_t = R_t + diag(((a^g)*ones(txAntennas(i)-g,1)),g) + diag(((a^g)*ones(txAntennas(i)-g,1)),-g); % generates diagionals in the antennna covariance matrix 
            end                                                                                                 % - this should be updated with the Bessel function
             for g = 1:rxAntennas(j)                                                                            % approach described in http://www.ece.nus.edu.sg/stfpage/elehht/Teaching/EE6832/Lecture%20Notes%5CMultiple%20Antennas%20for%20MIMO%20Communications%20-%20Channel%20Correlation.pdf
                R_r = R_r + diag(((a^g)*ones(rxAntennas(j)-g,1)),g) + diag(((a^g)*ones(rxAntennas(j)-g,1)),-g);
            end           
            H{i,j} = KroneckerChannel(txAntennas(i),rxAntennas(j),R_t,R_r);
            power = trace(H{i,j} * H{i,j}');
            H{i,j} = H{i,j} / sqrt(power);
        end
    end
    
else
    % no operation, return empty cell of channels
    
end

end

