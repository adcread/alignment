function [ H ] = generateChannel( users, txAntennas, rxAntennas, fading, channelModel )
%GENERATECHANNEL Summary of this function goes here
%   Detailed explanation goes here

H = cell(users,users);

if strcmp(channelModel,'scalar')
    

elseif strcmp(channelModel,'complex')


elseif strcmp(channelModel,'kronecker')

    a = 0.5;
    
    for i = 1:users
        for j = 1:users
            R_t = eye(txAntennas(j));
            R_r = eye(rxAntennas(i));
            for g = 1:txAntennas(j)
                R_t = R_t + diag(((a^g)*ones(txAntennas(j)-g,1)),g) + diag(((a^g)*ones(txAntennas(j)-g,1)),-g);
            end
             for g = 1:rxAntennas(i)
                R_r = R_r + diag(((a^g)*ones(rxAntennas(i)-g,1)),g) + diag(((a^g)*ones(rxAntennas(i)-g,1)),-g);
            end           
            H{i,j} = KroneckerChannel(txAntennas(j),rxAntennas(i),R_t,R_r);
        end
    end
    
else
    % no operation, return empty cell of channels
    
end

end

