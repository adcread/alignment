% 'Uses a combination of SVD and EVD to find the precoding matrices and
% ranks of each user in a 2-user MIMO interference channel with unequal
% numbers of Tx and Rx antennas.'

if ispc
    addpath('C:\PhD\alignment\Channel\');
    addpath('C:\PhD\alignment\General functions\');
elseif ismac
    addpath('/Users/chris/PhD/alignment/Channel/');
    addpath('/Users/chris/PhD/alignment/General functions/');
end

txAntennas = [3 2];
rxAntennas = [3 2];

H = generateChannel(2,rxAntennas,txAntennas,'kronecker');     % create a 2-user MIMO network with M=[3 2], N=[3 2].

%celldisp(H);

for i = 1:2
    for j = 1:2
        
        H_r{i,j} = H{j,i}';
        
    end
end

for i = 1:2
    for j = 1:2
        
        [Ux{i,j}, sigmax{i,j},Vx{i,j}] = svd(H{i,j});
        [Ux_r{i,j}, sigmax_r{i,j}, Vx_r{i,j}] = svd(H_r{i,j});
        
    end
end

SNR = [500 100 ; 100 500];  % slight attenuation on the cross links - doesn't really make a difference.

sequenceLength = 10240;

for i = 1:2
    
    sequence{i} = circSymAWGN(txAntennas(i),sequenceLength,1);
    
end

% message 1 - Transmitter 1 sends rank-3 training sequence

y{1} = sqrt(SNR(1,1)) * H{1,1} * eye(3) ;% * sequence{1};

y{2} = sqrt(SNR(1,2)) * H{1,2} * eye(3) ;%* sequence{1};

[U{1,1}, lambda{1,1}] = sortEigs(y{1}*y{1}');

lambda{1,1} = 1/(SNR(1,1) * sequenceLength) * lambda{1,1};

[U{1,2}, lambda{1,2}] = sortEigs(y{2}*y{2}');

lambda{1,2,1} = 1/(SNR(1,2) * sequenceLength) * lambda{1,2,1};

V_r{2,1} = U{1,2};

% message 2 - Receiver 2 sends rank-2 training sequence

y{1} = sqrt(SNR(2,1)) * H_r{2,1} * V_r{2,1}' ;%* sequence{2};

y{2} = sqrt(SNR(2,2)) * H_r{2,2} * V_r{2,1}' ;%* sequence{2};

[U_r{2,1}, lambda_r{2,1}] = sortEigs(y{1}*y{1}');

lambda_r{2,1} = 1/(SNR(2,1) * sequenceLength) * lambda_r{2,1};

[U_r{2,2}, lambda_r{2,2}] = sortEigs(y{2}*y{2}');

lambda_r{2,2} = 1/(SNR(2,2) * sequenceLength) * lambda_r{2,2};

sigma_r{2,1} = sqrt(lambda_r{2,1});

s = sigma_r{2,1}(1:rxAntennas(1),1:txAntennas(2));

V_r{2,1} = sqrt(pinv(s)*U_r{2,1}'*receivedSignal{1}*pinv(sequence{2}));

V_r{2,1} = V_r{2,1}';