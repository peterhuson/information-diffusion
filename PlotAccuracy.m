clear;

votes = readmatrix("digg2009/votes_714.csv");

xes = [1 8];
ts = [1 15];

x = xes(1):xes(2);
t = ts(1):ts(2);
[X, Y] = meshgrid(x,t);
Z = zeros(ts(2), xes(2));

time_0 = votes(1,1);
block = 1;
for v = votes'
    time_block = floor((v(1) - time_0) / (60 * 60)) + 1;
    if time_block == 1
        v(1) - 1245966166;
    end
    if v(4) > 0 && v(4) <= xes(2) && time_block <= ts(2) 
        if time_block > block
            Z(time_block,:) = Z(time_block-1,:);
            block = block + 1;
        end
        Z(time_block, v(4)) = Z(time_block, v(4)) + 1;
    end
end

% Normalize by the total population of the digg data
populations = [1 259 49014 1126905 2052095 2170024 2194024 2199247 2200477]; % Number of people in each distance
% [259 49014 1126905 2052095 2170024 2194024 2199247 2200477]
% Z = bsxfun(@rdivide,Z, populations);
Z(ts(1),:)
% Normalize by the final number of exposed people:
final_votes = Z(ts(2),:);
Z = bsxfun(@rdivide,Z, final_votes);
Z = Z .* 10000;

writematrix(Z, "Accuracy.txt")

m = mesh(X,Y,Z,'FaceAlpha','0.5','FaceColor','flat');
% set(gca,'ZScale','log');
xlabel("x Friend Distance"); ylabel("t Hour"); zlabel("z Global Density");
title("Story 714 Density");
