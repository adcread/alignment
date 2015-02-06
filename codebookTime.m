
steps = 1000;
iters = 50;

time = zeros(steps,iters);

for S = 1:steps
    for iter = 1:iters
        tic
        generateCodebook(1,S,1);
        time(S, iter) = toc;
    end
    timeAve(S) = mean(time(S,:));
end
