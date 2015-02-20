function [ b ] = orthBasis( x )
%ORTHBASIS Create a set of orthogonal bases spanning the space occupied by
%v
%   Detailed explanation goes here

dim = length(x);

b = zeros(dim,dim);

b(:,1) = x(:,1);

for i = 2:dim
    b(:,i) = x(:,i);
    sum = 0;
    for j = 1:i-1
        sum = sum + project(x(:,i),b(:,j));
    end
    b(:,i) = x(:,i) - sum;
end

for i = 1:dim
   b(:,i) = b(:,i)/sqrt(cpxIP(b(:,i),b(:,i)));
end

end

