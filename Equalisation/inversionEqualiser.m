
% Performs the channel inversion equalisation from [66] Weikert

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
        
        [U{i,j}, sigma{i,j},V{i,j}] = svd(H{i,j});
        [U_r{i,j}, sigma_r{i,j}, V_r{i,j}] = svd(H_r{i,j});
        
    end
end

SNR = 500;

sequenceLength = 10240;

for i = 1:2
    
    sequence{i} = circSymAWGN(txAntennas(i),sequenceLength,1);
    
end

receivedSignal{1} = sqrt(SNR) * eye(3)' * H{2,1} * V{2,1} * sequence{2};

receivedSignal{2} = sqrt(SNR) * eye(2)' * H{1,2} * V{1,2} * sequence{1};

[U_est{2,1}, sigma_est{2,1}] = sortEigs(receivedSignal{1}*receivedSignal{1}');

[U_est{1,2}, sigma_est{1,2}] = sortEigs(receivedSignal{2}*receivedSignal{2}');