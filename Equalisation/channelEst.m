addpath('C:\PhD\alignment');
addpath('C:\PhD\alignment\Channel');


users = 2;

txAntennas = [3 2];
rxAntennas = [3 2];

sequenceLength = 128;

H = generateChannel(users,txAntennas, rxAntennas,'kronecker');

for j = 1:users
    for i = 1:users
        d(i,j) = rank(H{i,j});
    end
end

for i = 1:users
    x{i} = circSymAWGN(txAntennas(i),sequenceLength,1);
end

y = cell(users,1);

for j = 1:users
    y{j} = zeros(rxAntennas(j),sequenceLength);
end

% for j = 1:users
%     for i = 1:users
%         y{j} = y{j} + H{i,j}*x{i};
%     end
% end

y{1} = H{1,1}*x{1};
y{2} = H{1,2}*x{1};

for j = 1:users
    P{j} = (y{j}*x{j}')/sequenceLength;
    V{j} = zeros(txAntennas(j),max(d(:,j)),sequenceLength);
    U{j} = zeros(rxAntennas(j),max(d(:,j)),sequenceLength);
end

% focus on user 1 - equalisation

W1{1}=zeros(rxAntennas(1),max(d(:,1)),sequenceLength);
W2{1}=zeros(rxAntennas(1),max(d(:,1)),sequenceLength);
W1{1}(:,:,1) = eye(txAntennas(1));
W1{1}(:,:,2) = eye(txAntennas(1));
W1{2}(:,:,1) = eye(txAntennas(1));
W1{2}(:,:,2) = eye(txAntennas(1));
user = 1;

for stream = 1:max(d(:,user))
       
    [V{user}(:,:,1), dummy] = eigs(sqrt(2)*(randn(rxAntennas(1))+ 1i*randn(rxAntennas(1))));
    [U{user}(:,:,1), dummy] = eigs(sqrt(2)*(randn(rxAntennas(1))+ 1i*randn(rxAntennas(1))));
    
    W1{1}(:,:,1) = eye(txAntennas(1)) * U{1}(:,:,1);
    W2{1}(:,:,1) = V{1}(:,:,1) * eye(rxAntennas(1));
    
    for iter = 2:sequenceLength

        for stream = 1:max(d(:,user))
            
            if stream == 1
                field = [2:max(d(:,user))];
            else
                field = [1:(stream-1) (stream+1):max(d(:,user))];
            end
            
            a = W1{user}(:,field,iter-1);
            b = W2{user}(:,field,iter-1);
            
            W1{user}(:,stream,iter) = (eye(rxAntennas(user)) - a * inv(a'*a) * a') * P{user}*V{user}(:,stream,(iter-1));
            U{user}(:,stream,iter) = W1{user}(:,stream,iter)*sqrt(inv((W1{user}(:,stream,iter)'*W1{user}(:,stream,iter))));
            
            W2{user}(:,stream,iter) = (eye(rxAntennas(user)) - b * inv(b'*b) * b') * P{user}'*U{user}(:,stream,iter);
            V{user}(:,stream,iter) = W2{user}(:,stream,iter)*sqrt(inv((W2{user}(:,stream,iter)'*W2{user}(:,stream,iter))));
            
            sigma(:,stream,iter) = sqrt(inv((W2{user}(:,stream,iter)'*W2{user}(:,stream,iter))));
            
            H_est{1}(:,:,iter) = U{1}(:,:,iter)*diag(sigma(:,:,iter))*V{1}(:,:,iter)';
            
            err(iter) = norm(H{1,1} - H_est{1}(:,:,iter),'fro');
        end
    end
end
