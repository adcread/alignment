function [ ] = plotSignals( signals )
%PLOTSIGNALS Summary of this function goes here
%   Detailed explanation goes here

vector = cell(1,length(signals));

for i=1:length(signals)
    vector{i} = [0 signals(i)];
end
for i=1:length(signals)
    plot(vector{i});
    hold all;
end

end

