%Gold code sequence generator
%Run from editor Debug (F5)
%JC 4/25/09
%Gold code generation can sometimes be confusing and hopefully this m-file
%will be helpful. The m-file uses two preferred pairs of m-sequences (length 2^n-1)
%chips long where n=5 with length 31. The 33 gold codes
%(two original m-sequences plus 2^n-1 gold codes) produced in a 31x33 matrix
%here (with n=5) are probably not to useful for multiple access situations (CDMA)
%but is used here to provide a model for learning. The preferred pair used
%here is 52 and 5432 or 45 and 75 octal. Other preferred pairs can include values
%such as 53 5321, 5431, and 5421. The m-file will also check the cross-correlation
%values between any two codes in the matrix and should be three valued%(-9,-1,7). The 
%autocorrelation of a Gold code can also be checked. Some, not all, of the
%Gold codes in the matrix are balanced (16 1's and 15 -1's).

function [ goldCodeMatrix ] = generateGoldCodes( order )

    length = (2^order)-1;

    % Pick a known preferred pair from lookup table
    switch order
        case 3
            error('No such generator polynomials exist');
        case 4
            error('No such generator polynomials exist');
        case 5
            PNtaps1 = [0 1 0 0 1]; % z^5 + z^2 + 1
            PNtaps2 = [0 1 1 1 1]; % z^5 + z^4 + z^3 + z^2 + 1
        case 6
            PNtaps1 = [0 1 0 0 0 1];
            PNtaps2 = [1 1 0 0 1 1];
        case 7
            PNtaps1 = [0 0 1 0 0 0 1];
            PNtaps2 = [1 1 1 0 0 0 1];
        case 8
            PNtaps1 = [1 1 0 0 1 1 1 1];
            PNtaps2 = [1 0 0 0 0 1 1 1];
        case 9
            PNtaps1 = [0 0 0 1 0 0 0 0 1];
            PNtaps2 = [0 0 1 1 0 1 0 0 1];
        case 10
            PNtaps1 = [0 0 1 1 1 1 1 1 1 1];
            PNtaps2 = [1 0 0 1 0 1 1 0 1 1];
        case 11
            PNtaps2 = [0 1 0 0 0 0 0 0 0 0 1];
            PNtaps1 = [0 1 0 0 1 0 0 1 0 0 1];
        case 12
            PNtaps1 = [0 1 1 0 0 0 0 0 1 0 0 1];
            PNtaps2 = [1 1 1 0 0 1 1 1 0 0 1 1];
        case 13
            PNtaps1 = [1 0 1 1 0 0 0 0 0 0 0 0 1];
            PNtaps2 = [0 0 0 1 1 0 1 0 1 1 0 0 1];
    end              
    
    %Generate 1st m-sequence of length 2^(order)-1
    
    PNregister1 = ones(1,order);
    PNcode1=zeros(1,length);
    
    for i=1:length
       temp = mod(sum(and(PNtaps1,PNregister1)),2);
       PNcode1(i) = PNregister1(order);
       for j=order:-1:2
          PNregister1(j)=PNregister1(j-1);
       end
       PNregister1(1) = temp;
    end
    
    %Generate 2nd m-sequence
    
    PNregister2=ones(1,order);%initial fill
    PNcode2=zeros(1,length);
    
    for i=1:length
       temp = mod(sum(and(PNtaps2,PNregister2)),2);
       PNcode2(i) = PNregister2(order);
       for j=order:-1:2
          PNregister2(j)=PNregister2(j-1);
       end
       PNregister2(1) = temp;
    end

    %Generate a set of unique Gold codes in a matrix (31x33) which includes
    %the original 1st and 2nd m-sequences plus (2^n-1)=31 other Gold codes with n=5. 
    %These unique codes are with initial fills of [1 1 1 1 1] in register1
    %and register2. Other unique sets of Gold codes can be generated with 
    %different initial fill values.

    goldCodeMatrix = zeros(length,length+1);
    
    goldCodeMatrix(:,1) = PNcode1;
    goldCodeMatrix(:,2) = PNcode2;
    
    for lag=0:length
        shiftedCode = circshift(PNcode2,lag,2);
        goldCodeMatrix(:,3+lag) = xor(PNcode1,shiftedCode).';
    end

    %Change matrix codes from 1/0 to 1/-1 and show 33 codes in command window.
    
    goldCodeMatrix = 2*goldCodeMatrix-1;

end