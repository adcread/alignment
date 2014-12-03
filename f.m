function [ result ] = f( u, a1, u1, a2, u2 )
%F Summary of this function goes here
%   Detailed explanation goes here

    if (a1 >=a2)
        result = max([0 (min([u u1])*a1)]) + max([0 (min([(u-u1) u2])*a2)]);
    else
        result = max([0 (min([u u2])*a2)]) + max([0 (min([(u-u2) u1])*a1)]);
    end

end

