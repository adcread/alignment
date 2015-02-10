function [ userPrivate, userCommon, crossPrivate, crossCommon ] = universalDegreesOfFreedom2UserMIMOIC( M, N, alpha, beta )
%DEGREESOFFREEDOM Summary of this function goes here
%   Detailed explanation goes here

limit = max([M N]);
users = size(alpha,1);

for i = 1:users
    for j = 1:users
        m(i,j) = min([M(i) N(j)]);
    end
end

c1 = f(N(1),beta(1,2),m(1,2),alpha(1,1),pos((M(1)-N(2))))*alpha(1,1);
c2 = min([N(1) M(1) N(2)])*alpha(1,1);
c3 = (min([N(1) M(2)]))*alpha(2,1);
c4 = min([M(1) N(1)]) * alpha(1,1);
c5 = g(N(1),alpha(2,1),M(2),beta(1,2),m(1,2),1,pos((M(1)-N(2))));
c6 = f(N(1),alpha(2,1),M(2),alpha(1,1),m(1,2));
c7 = f(N(1),alpha(2,1),M(2),alpha(1,1),M(1));
c8 = f(N(2),beta(2,1),m(2,1),alpha(2,2),pos((M(2)-N(1))))*alpha(2,2);
c9 = min([N(2) M(2) N(1)])*alpha(2,2);
c10 = (min([N(2) M(1)]))*alpha(1,2);
c11 = min([M(2) N(2)]) * alpha(2,2);
c12 = g(N(2),alpha(1,2),M(1),beta(2,1),m(2,1),1,pos((M(2)-N(1))));
c13 = f(N(2),alpha(1,2),M(1),alpha(2,2),m(2,1));
c14 = f(N(2),alpha(1,2),M(1),alpha(2,2),M(2));
target = [-1 -1 -1 -1];
constraints = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14];
coefficients = [alpha(1,1) 0 0 0; 
    0 alpha(1,1) 0 0; 
    0 0 0 alpha(2,2) ;
    alpha(1,1) alpha(1,1) 0 0; 
    alpha(1,1) 0 0 alpha(2,2);
    0 alpha(1,1) 0 alpha(2,2) ; 
    alpha(1,1) alpha(1,1) 0 alpha(2,2);
    0 0 alpha(2,2) 0 ; 
    0 0 0 alpha(2,2); 
    0 alpha(1,1) 0 0;
    0 0 alpha(2,2) alpha(2,2); 
    0 alpha(1,1) alpha(2,2) 0;
    0 alpha(1,1) 0 alpha(2,2); 
    0 alpha(1,1) alpha(2,2) alpha(2,2)];

options = optimset('Display', 'iter', 'Diagnostics', 'on', 'Simplex', 'on', 'Algorithm','simplex');

lowerBounds = zeros(1,4);
upperBounds = [c1 c2 c8 c9];

[result, sumDoF, exitflag] = linprog(target, coefficients,constraints,[],[],lowerBounds,upperBounds,[],options);
userPrivate = result(1);
userCommon = result(2);
crossPrivate = result(3);
crossCommon = result(4);
end

