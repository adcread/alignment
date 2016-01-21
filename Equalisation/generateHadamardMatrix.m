function [H]=generateHadamardMatrix(codeSize)
 
%[H]=generateHadamardMatrix(codeSize);
 
% Function to generate Walsh-Hadamard Matrix where "codeSize" is the code
 
% length of walsh code. The first matrix gives us two codes; 00, 01. The second
 
% matrix gives: 0000, 0101, 0011, 0110 and so on
 
% Author: Mathuranathan for http://www.gaussianwaves.com
 
% License: Creative Commons: Attribution-NonCommercial-ShareAlike 3.0
 
% Unported
 
 
%codeSize=64; %For testing only
 
N=2;
H=[0 0 ; 0 1];
if bitand(codeSize,codeSize-1)==0
while(N~=codeSize)
       N=N*2;
       H=repmat(H,[2,2]);
       [m,n]=size(H);
 
      %Invert the matrix located at the bottom right hand corner
      for i=m/2+1:m,
          for j=n/2+1:n,
                H(i,j)=~H(i,j);
         end
     end
end
else
disp('INVALID CODE SIZE:The code size must be a power of 2');
end