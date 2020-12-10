function [accuracy] = LogisticalDiffusiveAccuracy(delta,r,K,plot_bool)
%UNTITLED delta, r, K
%   Detailed explanation goes here
% Logistical
% d = 0.002;
% r = 0.4;
% K = 25;
% txt = sprintf('%f | %f | %f', delta,r,K)
if delta <= 0
    accuracy = 1000;
    return
end
    
% Driving Parameters
numxpoints = 101;
xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
xx = linspace(xes(1),xes(2),numxpoints);
% Initial condition from paper i.e. \phi(x)
% y = [0 6 2 2 1 1 1 0.5 0];
y = [101 560 1085 343 76 18 2 0];
I_initial = spline(x,[0 y 0]);

Z = zeros(ts(2), numxpoints);
Z_model = zeros(ts(2), xes(2));

Ispline = I_initial;
for t_i = t
    Z(t_i,:) = ppval(Ispline,xx);
    Z_model(t_i,:) = ppval(Ispline,x);
%     I_t = delta * ppval(fnder(Ispline,2),x) + prod(t_i,:).*ppval(Ispline,x);
%     I_t = spline(x,[0 0 I_t(2:xes(2)) 0 0]);
%     t = fn2fm(fncmb(Ispline,1/K),'B-');
    vals = ppval(Ispline,xx);
    ik = spline(xx,((vals.*(1-(vals./K)))).*r);
    
    ixxr = fncmb(fnder(Ispline,2),delta);
    smooth = xx(1:4:numxpoints);
    ixx = spline(smooth,ppval(ixxr,smooth));
    I_t = fncmb(ixx,'+',...
        ik ... (1-I/K)
        );
%     figure(1);
%     fnplt(Ispline);
%     hold on;fnplt(ik);fnplt(ixxr);fnplt(ixx);fnplt(I_t);legend("I","Ik","Ixxr","Ixx","I_t");
%     hold off;
%     set(gca, 'YScale', 'log')
    
    Ispline = fncmb(I_t,'+',Ispline);
end


if plot_bool
    [X, Y] = meshgrid(xx,t);
    figure(3);
    mesh(X,Y,Z,'FaceAlpha','0.6','EdgeAlpha','0.5','FaceColor','interp');
    view(-80,20);
    xlabel("x Distance");
    ylabel("t Time");
    zlabel("z Votes");
    title("Logistical Diffusive Information Diffusion");
    hold on;
    xlim(xes)
    ylim(ts)
    zlim([0 7000])
end

file = "digg2009/votes_714.csv";
votes = readmatrix(file);

xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
Z_digg = zeros(ts(2), xes(2));

time_0 = votes(1,1);
block = 1;
for v = votes'

    time_block = floor((v(1) - time_0) / (60 * 60)) + 1;
%     if time_block > 0
%         v
%     end
    if v(4) > 0 && v(4) <= xes(2) && time_block <= ts(2)
        if time_block > block
            Z_digg(time_block,:) = Z_digg(time_block-1,:);
            block = block + 1;
        end
        Z_digg(time_block, v(4)) = Z_digg(time_block, v(4)) + 1;
    end
end

if plot_bool
    [X, Y] = meshgrid(x,t);
    
    % writematrix(Z, "Accuracy.txt")
    % m = mesh(X,Y,Z,[0 0 0],'LineStyle','--','FaceAlpha','0.2','EdgeAlpha','0','EdgeColor','white');%,'FaceColor','interp'
    X = reshape(X.',1,[]);
    Y = reshape(Y.',1,[]);
    Z_digg_plot = reshape(Z_digg.',1,[]);
    figure(3);
    scatter3(X,Y,Z_digg_plot,'black');%,'FaceColor','interp'
    hold off;
end

accuracy = CalcAccuracy(Z_model(:,2:4), Z_digg(:,2:4));
accuracy = relError(Z_model(:,2:4), Z_digg(:,2:4),'display',false);
accuracy = accuracy(2,1);

end