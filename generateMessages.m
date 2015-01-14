privateMessageStream = cell(2,1);
publicMessageStream = cell(2,1);

for i = 1:users
    privateMessageStream{i} = randi([0 3],100,1);
    publicMessageStream{i} = randi([0 3],100,1);
end