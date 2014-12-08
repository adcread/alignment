function [ result ] = g( u0, a1, u1, a2, u2, a3, u3 )
%G Summary of this function goes here
%   Detailed explanation goes here

a = [a1 a2 a3];
u = [u1 u2 u3];

[ai,sortVector] = sort(a,'descend');

for i =1:3
    ui(i) = u(sortVector(i));
end

result = pos(min([u0, ui(1)])*ai(1)) + pos(min([pos(u0-ui(1)) ui(2)])*ai(2)) + pos(min([pos(u0-ui(1)-ui(2)) ui(3)])*ai(3));
end

