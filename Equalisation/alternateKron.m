[m,n] = size(Z);

m1 = m / rxAntennas;
m2 = rxAntennas;

n1 = n / rxAntennas;
n2 = rxAntennas;

C = eye(m2,n2);
B = eye(m1,n1);

gamma = trace(C.'*C);

for i = 1:m1
    for j = 1:n1
        B(i,j) = trace(C.'*Z(((i-1)*m2+1):(i*m2),((j-1)*n2+1):(j*n2)))/gamma;
    end
end
beta = trace(B.'*B);
for i = 1:m2
    for j = 1:n2
        C(i,j) = trace(B.'*Z(((i-1)*m1+1):(i*m1),((j-1)*n1+1):(j*n1)))/beta;
    end
end