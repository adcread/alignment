function [ userPrivate, userCommon, crossCommon ] = degreesOfFreedom2UserMIMOIC_Memory( user, maxIter, M, N, alpha, beta, crossCommonIn, crossWeight, previousUserCommon, previousWeight )
%degreesOfFreedom2UserMIMOIC_Memory Calculates DoF for user in 2-user MIMO IC

%   Detailed explanation goes here

users = size(alpha,1);

for i = 1:users
    for j = 1:users
        m(i,j) = min([M(i) N(j)]);
    end
end

if (user==1)
    cross = 2;
else
    cross = 1;
end

g1 = m(user,user);
g2 = crossCommonIn;
g3 = f(N(user),beta(user,cross),m(user,cross),alpha(user,user),pos((M(user)-N(cross))))*alpha(user,user); % optimise private DoF
g4 = previousUserCommon;

c1 = f(N(user),beta(user,cross),m(user,cross),alpha(user,user),pos((M(user)-N(cross))))*alpha(user,user);
c2 = min([N(user) M(user) N(cross)])*alpha(user,user);
c3 = (min([N(user) M(cross)]))*alpha(cross,user);
c4 = min([M(user) N(user)]) * alpha(user,user);
c5 = g(N(user),alpha(cross,user),M(cross),beta(user,cross),m(user,cross),1,pos((M(user)-N(cross))));
c6 = f(N(user),alpha(cross,user),M(cross),alpha(user,user),m(user,cross));
c7 = f(N(user),alpha(cross,user),M(cross),alpha(user,user),M(user));

target = [0 0 0 0 1 0 crossWeight 0 2 0 previousWeight];                                              % weight determines the severity of reducing the other users's common rate
ineqConstraints = [c1 c2 c3 c4 c5 c6 c7];
% Coefficients = d1p d2c d2c y1p y1m y2p y2m y3p y3m
ineqCoefficients = [alpha(user,user) 0 0 0 0 0 0 0 0 0 0; 
    0 alpha(user,user) 0 0 0 0 0 0 0 0 0; 
    0 0 alpha(cross,cross) 0 0 0 0 0 0 0 0;
    alpha(user,user) alpha(user,user) 0 0 0 0 0 0 0 0 0 ; 
    alpha(user,user) 0 alpha(cross,cross) 0 0 0 0 0 0 0 0;
    0 alpha(user,user) alpha(cross,cross) 0 0 0 0 0 0 0 0; 
    alpha(user,user) alpha(user,user) alpha(cross,cross) 0 0 0 0 0 0 0 0];

eqConstraints = [g1 g2 g3 g4];
eqCoefficients = [1 1 0 -1 1 0 0 0 0 0 0;
    0 0 1 0 0 -1 1 0 0 0 0;
    1 0 0 0 0 0 0 -1 1 0 0;
    0 1 0 0 0 0 0 0 0 0 0];

options = optimset('Display', 'none', 'Diagnostics', 'off', 'Simplex', 'on', 'Algorithm','Simplex', 'MaxIter', maxIter);

lowerBounds = [0 0 0 0 0 0 0 0 0 0 0];
upperBounds = [c1 c2 min([c3 crossCommonIn]) Inf Inf Inf Inf Inf Inf Inf Inf];

[result] = linprog(target, ineqCoefficients,ineqConstraints,eqCoefficients,eqConstraints,lowerBounds,upperBounds,[],options);
userPrivate = result(1);
userCommon = result(2);
crossCommon = result(3);
end

