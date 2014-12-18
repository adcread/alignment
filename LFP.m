
%% maximise Z = \frac{x_1 + 3x_2 + 2x_3}{2x_1 + x_2 + 4x_3 + 1}

% subject to x_1 + 3x_2 + 6x_3 \leq 8
% 2x_1 + x_2 + 4x_3 \leq 5

c = [1 3 2];
alpha = [];
d = [2 1 4];
beta = 1;

A = [1 3 6; 2 2 3];
b = [8 5];

