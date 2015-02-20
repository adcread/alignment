function [ output ] = convertBases( input, newspace )
%CONVERTBASES Summary of this function goes here
%   Detailed explanation goes here

dim = length(input);

output = zeros(dim,1);

for i = 1:dim
   output(i) = cpxIP(input,newspace(:,i));
end

end

