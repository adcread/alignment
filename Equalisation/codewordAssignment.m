b = 5;

l = 2^b;

codebook = qammod(0:l-1,l);

scatterplot(codebook);

for i = 1:l
    x = real(codebook(i));
    y = imag(codebook(i))-0.15;
    text(x,y,num2str(i-1));
end
