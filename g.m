function [ result ] = g( u, a1, u1, a2, u2, a3, u3 )
%G Summary of this function goes here
%   Detailed explanation goes here

result = pos(min([u, u1])*a1) + pos(min([pos(u-u1) u2])*a2) + pos(min([pos(u-u1-u2) u3])*a3);
end

