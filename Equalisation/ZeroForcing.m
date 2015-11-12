
y = receivedMessage{1} - publicInterference{1};

H_bar = H{1,1}*V{1,2};

G = orthBasis(circshift(H_bar,1,2));

G(:,1) = 0;

G = circshift(G,2,2);

y_1 = G'*y;

y_2 = U{2,1}'*y_1;

z = 1/sqrt(SNR(2,1)) * pinv(directionPub{2,1}) * pinv(sigma{2,1}) * y_2;