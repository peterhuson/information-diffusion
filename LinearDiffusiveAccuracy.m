function [accuracy] = LinearDiffusiveAccuracy(alpha,Beta,gamma,delta,rho,sigma,y3,plot_bool)
%LINEARDIFFUSIVEACCURACY (alpha,Beta,gamma,delta,rho,sigma,plot_bool)
%   Detailed explanation goes here
txt = sprintf('%f | %f | %f | %f | %f | %f | %f', alpha,Beta,gamma,delta,rho,sigma,y3)

% Driving Parameters
numxpoints = 101;
xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
% Initial condition from paper i.e. \phi(x)
% y = [0 6 2 2 1 1 1 0.5 0];
% y = [101 560 1085 343 76 18 2 0];
y = [101 560 y3 343 76 18 2 0];
xx = linspace(xes(1),xes(2),numxpoints);

% Linear Diffusive
%%% h(x)
%TODO: Use Mary's distance data to fit a more accurate h(x) function
h = -(x - rho).*(x - sigma);

asymptote = Beta / alpha;
r = asymptote - exp(-alpha*(t - 1))*(asymptote - gamma);

prod = r'*h;
% Show growth mesh

% The key to this implementation is that you can easily take the derivative
% of a spline using `fnder()`
I_initial = spline(x,[0 y 0]);

% Mesh stuff
Z = zeros(ts(2), numxpoints);
Z_model = zeros(ts(2), xes(2));

if plot_bool
    [X_p, Y_p] = meshgrid(x,t);
    [X, Y] = meshgrid(xx,t);
    
%     figure(1);
%     mesh(X_p,Y_p,prod,'FaceAlpha','0.5','FaceColor','flat');
%     xlabel("x Distance"); ylabel("t Time"); zlabel("z Growth");
%     title("Growth Function");
end

Ispline = I_initial;
for t_i = t
    Z(t_i,:) = ppval(Ispline,xx);
    Z_model(t_i,:) = ppval(Ispline,x);
%     I_t = delta * ppval(fnder(Ispline,2),x) + prod(t_i,:).*ppval(Ispline,x);
%     I_t = spline(x,[0 0 I_t(2:xes(2)) 0 0]);
    I_t = fncmb(fncmb(fnder(Ispline,2),delta),'+', ...
        spline(xx,ppval(spline(x, prod(t_i,:)),xx).*ppval(Ispline,xx)));
    Ispline = fncmb(I_t,'+',Ispline);
end

if plot_bool
    figure(3);
    view(-80,20)
    mesh(X,Y,Z,'FaceAlpha','0.6','EdgeAlpha','0.5','FaceColor','interp')%,'FaceColor','flat')
    xlabel("x Distance");
    ylabel("t Time");
    zlabel("z Votes");
    title("Linear Diffusive Information Diffusion");
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
    view(-80,20)
    figure(3);
    scatter3(X,Y,Z_digg_plot,'black');%,'FaceColor','interp'
    view(30,20)
    hold off;
end

accuracy = CalcAccuracy(Z_model(:,1:5), Z_digg(:,1:5));
accuracy = relError(Z_model(:,1:5), Z_digg(:,1:5),'display',false);
accuracy = accuracy(2,1);
end

