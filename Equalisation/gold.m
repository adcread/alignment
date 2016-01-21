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


%Generate 1st m-sequence (52) of length 31 (1+x^2+x^5)
register1=[1 1 1 1 1];%initial fill
code1=zeros(1,31);
for i=1:31
   temp = mod(register1(2)+register1(5),2);
   code1(i) = 2*register1(5)-1;
   for j=5:-1:2
      register1(j)=register1(j-1);
   end
   register1(1) = temp;
end
%Generate 2nd m-sequence (5432) of length 31 (1+x^2+x^3+x^4+x^5)
register2=[1 1 1 1 1];%initial fill
code2=zeros(1,31);
for i=1:31
   temp = mod(register2(2)+register2(3)+register2(4)+register2(5),2);
   code2(i) = 2*register2(5)-1;
   for j=5:-1:2
      register2(j)=register2(j-1);
   end
   register2(1) = temp;
end

m_sequence_1=code1;%1/-1(bipolar sequence) output
m_sequence_2=code2;%1/-1(bipolar sequence) output


m_sequence_1=m_sequence_1';%transpose to a column
m_sequence_2=m_sequence_2';%transpose to a column


%Generate a set of unique Gold codes in a matrix (31x33) which includes
%the original 1st and 2nd m-sequences plus (2^n-1)=31 other Gold codes with n=5. 
%These unique codes are with initial fills of [1 1 1 1 1] in register1
%and register2. Other unique sets of Gold codes can be generated with 
%different initial fill values.


m_sequence_1 =m_sequence_1>0;%change 1/-1 to 1/0 
m_sequence_2 =m_sequence_2>0;%change 1/-1 to 1/0

Gold_code_matrix(:,1) = m_sequence_1;
Gold_code_matrix(:,2) = m_sequence_2;

for phase_shift=0:30
    shifted_code=circshift(m_sequence_2,phase_shift);
    Gold_code_matrix(:,3+phase_shift)=mod(m_sequence_1+shifted_code,2);
end

%Change matrix codes from 1/0 to 1/-1 and show 33 codes in command window.
Gold_code_matrix=2*Gold_code_matrix-1%change 1/0 to 1/-1

%Choose 2 codes from Gold code matrix and plot the cross-correlation
%values. If codes are from preferred pairs, the three values would be
%(-9,-1,+7).

codeA=Gold_code_matrix(:,9);
codeB=Gold_code_matrix(:,22);

%Determine cross-correlations

for shift=0:40
    shifted_code1 = circshift(codeA,shift);
    crosscorrelation(shift+1) = codeB'*shifted_code1;
end

subplot(3,2,1)
plot(crosscorrelation) 
grid on
xlabel('shifts');ylabel('value of correlations');
title('Cross-correlations of two 31 chip codes')

%Choose 1 code from the Gold code matrix and plot the autocorrelation
%values.

codeC=Gold_code_matrix(:,17);

for shift=0:40
    shifted_code_A= circshift(codeC,shift);
    autocorrelation_1(shift+1) = codeC'*shifted_code_A;     
end

%Show that all autocorrelation values of m-sequence 1 and 2
%(with non-zero shift) equals 31/-1.  

subplot(3,2,2)
plot(autocorrelation_1)
grid on
xlabel('shifts');ylabel('value of correlations');
title('31 chip Gold code autocorrelation')

%References

%[1] R. C. Dixon, Spread Spectrum Systems Applications, 2nd ed, Wiley (1984)

%[2] Matlab Helpdesk, Gold Sequence Generator,
%www.mathworks.com/access/helpdesk/help/toolbox/commblks/ref/goldsequencegenerator.html

%[3] Peterson's Table of Irreducible Polynomials over GF(2),
%www.cs.umbc.edu/~lomonaco/f97/442/Peterson_Table.html

%[4] John G. Prokis, Digital Communications, 4th ed, 13.2.5 Generation of
%PN Sequences and Gold codes,pages 776-771