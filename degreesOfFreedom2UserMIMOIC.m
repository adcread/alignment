function [ userPrivate, userCommon, crossCommon ] = degreesOfFreedom2UserMIMOIC( user, M, N, alpha, beta, crossCommonIn )
%DEGREESOFFREEDOM Summary of this function goes here
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

c1 = f(N(user),beta(user,cross),m(user,cross),alpha(user,user),pos((M(user)-N(cross))))*alpha(user,user);
c2 = min([N(user) M(user) N(cross)])*alpha(user,user);
c3 = (min([N(user) M(cross)]))*alpha(cross,user);
c4 = min([M(user) N(user)]) * alpha(user,user);
c5 = g(N(user),alpha(cross,user),M(cross),beta(user,cross),m(user,cross),1,pos((M(user)-N(cross))));
c6 = f(N(user),alpha(cross,user),M(cross),alpha(user,user),m(user,cross));
c7 = f(N(user),alpha(cross,user),M(cross),alpha(user,user),M(user));
c8 = crossCommonIn;                                                         % crossCommonIn is the public DoF limit from the other user.

target = [-1 -1 -1];
constraints = [c1 c2 c3 c4 c5 c6 c7 c8];
coefficients = [alpha(user,user) 0 0 ; 
    0 alpha(user,user) 0 ; 
    0 0 alpha(cross,cross) ;
    alpha(user,user) alpha(user,user) 0 ; 
    alpha(user,user) 0 alpha(cross,cross);
    0 alpha(user,user) alpha(cross,cross) ; 
    alpha(user,user) alpha(user,user) alpha(cross,cross);
    0 0 1];

options = optimset('Display', 'iter', 'Diagnostics', 'off', 'Simplex', 'on', 'Algorithm','Simplex');

lowerBounds = [0 0 0];
upperBounds = [c1 c2 min([c3 c8])];

[result] = linprog(target, coefficients,constraints,[],[],lowerBounds,upperBounds,[],options);
userPrivate = result(1);
userCommon = result(2);
crossCommon = result(3);
end

