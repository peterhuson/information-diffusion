clear;

votes = readmatrix("digg2009/votes_714.csv");

xes = [0 8];
ts = [1 15];

x = xes(1):xes(2);
t = ts(1):ts(2);
[X, Y] = meshgrid(x,t);
Z = zeros(ts(2), xes(2)+1);

time_0 = votes(1,1);
block = 1;
for v = votes'
    
    time_block = floor((v(1) - time_0) / (60 * 60)) + 1;
%     if time_block > 0
%         v
%     end
    if v(4)+1 <= xes(2) && time_block <= ts(2) 
        if time_block > block
            Z(time_block,:) = Z(time_block-1,:);
            block = block + 1;
        end
        Z(time_block, v(4)+1) = Z(time_block, v(4)+1) + 1;
    end
end

populations = [1 259 49014 1126905 2052095 2170024 2194024 2199247 2200477]; % Number of people in each distance
% [259 49014 1126905 2052095 2170024 2194024 2199247 2200477]

Z = bsxfun(@rdivide,Z, populations);

m = mesh(X,Y,Z,'FaceAlpha','0.5','FaceColor','flat');
set(gca,'ZScale','log');
xlabel("x Friend Distance"); ylabel("t Hour"); zlabel("z Global Density");
title("Story 714 Density");
