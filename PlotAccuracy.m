clear;
file = "digg2009/votes_714.csv";
votes = readmatrix(file);

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
%     if time_block > 0
%         v
%     end
    if v(4) > 0 && v(4) <= xes(2) && time_block <= ts(2)
        if time_block > block
            Z(time_block,:) = Z(time_block-1,:);
            block = block + 1;
        end
        Z(time_block, v(4)) = Z(time_block, v(4)) + 1;
    end
end

% Normalize by the total population of the digg data
% populations = [259 49014 1126905 2052095 2170024 2194024 2199247 2200477]; % Number of people in each distance
% populations = [259 49014 1126905 2052095 2170024 2194024 2199247 2200477]; % Number of people in each distance
% populations = [359 99017 1214129 1885767 1990504 2014729 2020741 2022199]; 
% populations = [1020 273388 1571443 1919342 1997461 2016886 2021359 2022396]; % Number of people in each distance
populations = [1557 346744 1641109 1937290 2003317 2018062 2021615 2022430];
populations = [populations(1) diff(populations)];
% figure(1);
% plot(x, populations(1:xes(2)));
% [259 49014 1126905 2052095 2170024 2194024 2199247 2200477]
% Z = bsxfun(@rdivide,Z,populations(1:xes(2)));

% Normalize by the final number of exposed people:
final_votes = Z(ts(2),:);
Z = bsxfun(@rdivide,Z, final_votes);

writematrix(Z, "Accuracy.txt")
figure(2);
m = mesh(X,Y,Z,'FaceAlpha','0.5','FaceColor','flat');
% set(gca,'ZScale','log');
xlabel("x Friend Distance"); ylabel("t Hour"); zlabel("z Global Density");
title(file);
