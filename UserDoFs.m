% varies the side information DoF for Users 1 and 2 and plots their optimal
% responses.

M = [3 2];
N = [3 2];
alpha = [1 0.6 ; 0.6 1];
beta = [0 0.4 ; 0.4 0];

x = 0:0.1:3;

for n = 1:length(x)

    [d1p(n),d1c(n),d2c1(n)] = degreesOfFreedom2UserMIMOIC_LGP(1,100,M,N,alpha,beta,x(n),4);
    [d2p(n),d2c(n),d1c2(n)] = degreesOfFreedom2UserMIMOIC_LGP(2,100,M,N,alpha,beta,x(n),4);

end

sum1 = d1p + d1c;
sum2 = d2p + d2c;

figure;
subplot(2,1,1);
plot(x,d1p);
hold on;
plot(x,d1c);
plot(x,d2c1);
plot(x,sum1);
legend('d1p','d1c','d2c1','sum');
xlabel('d2c2');
ylabel('DoF');

subplot(2,1,2);
plot(x,d2p);
hold on;
plot(x,d2c);
plot(x,d1c2);
plot(x,sum2);
legend('d2p','d2c','d1c2','sum');
xlabel('d1c1');
ylabel('DoF');