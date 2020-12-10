clear;
% Driving Parameters
numxpoints = 101;
xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
% Initial condition from paper i.e. \phi(x)
% y = [0 6 2 2 1 1 1 0.5 0];
y = [101 560 1085 343 76 18 2 0];
y3 = 1085;
xx = linspace(xes(1),xes(2),numxpoints);

% Linear Diffusive
%%% h(x)
%TODO: Use Mary's distance data to fit a more accurate h(x) function
rho = -0.9478;
sigma = 8.9149;
h = -(x - rho).*(x - sigma);
delta = 0.0001;

% %%% r(t)
% Beta = 0.00005;
% alpha = 0.3526;
bs = [0.00000:0.0002:0.002];
bs = [0.0001];
b_range = [1:numel(bs)];
as = [0.1:0.2:4.5];
as = [0.4428    ];
a_range = [1:numel(as)];
gs = [0.01:0.01:0.1]; % Fix a and iterate over g and B
gs = [0.0283    ];
g_range = [1:numel(gs)];
[X_a, Y_a] = meshgrid(gs,bs);
accuracy = zeros(numel(bs), numel(gs));
a_i = 1;
for b_i = b_range
    for g_i = g_range
        Beta = bs(b_i);
        alpha = as(a_i);
        gamma = gs(g_i);
        accuracy(b_i, g_i) = LinearDiffusiveAccuracy(alpha,Beta,gamma,delta,rho,sigma,y3,true)
    end
end
figure(3);
accuracy
% mesh(X_a,Y_a,accuracy,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x \gamma");
ylabel("y \Beta");

options = optimset('PlotFcns',@optimplotfval);
fun = @(x)LinearDiffusiveAccuracy(x(1),x(2),x(3),x(4),x(5),x(6),x(7),true);

x5 = [0.42,0.000086,0.041,0.00014,-0.5,7.39,1055.3]; % Final good run
x6 = [0.413,0.000034,0.0181,-0.000042,-1.88,9.876,1093.3]; % 2.66 relative error
x7 = [0.413049,0.000032,0.017005,-0.000052,-2.01397,10.1360,1094.79]; % 2.6642 relative error


x,f = fminsearch(fun,x7,options)
fun(f)
dim = [0.7 0.65 .15 .2];
str = sprintf('alpha=%.4f \n beta=%.5f \n gamma=%.3f \n d=%.5f \n rho=%.3f \n sigma=%.2f', f(1),f(2),f(3),f(4),f(5),f(6));
annotation('textbox',dim,'String',str)

% relError(Z_model, Z_digg)

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
% final_votes = Z_digg(ts(2),:);
% Z = bsxfun(@rdivide,Z, final_votes);

% writematrix(Z, "Accuracy.txt")
% m = mesh(X,Y,Z,[0 0 0],'LineStyle','--','FaceAlpha','0.2','EdgeAlpha','0','EdgeColor','white');%,'FaceColor','interp'
% X = reshape(X.',1,[]);
% Y = reshape(Y.',1,[]);
% Z_digg = reshape(Z_digg.',1,[]);
% figure(2);
% m = scatter3(X,Y,Z_digg);%,'FaceColor','interp'

hold off;
%% Logistical


options = optimset('PlotFcns',@optimplotfval);
fun = @(x)LogisticalDiffusiveAccuracy(x(1),x(2),x(3),true);
guess1 = [0.00001 0.062 7972];
accuracy,f = fminsearch(fun,guess1,options)
fun(f)
dim = [0.7 0.65 .15 .15];
str = sprintf('d=%.5f \n r=%4f \n K=%.2f', f(1),f(2),f(3));
annotation('textbox',dim,'String',str)


%%% Debugging Stuff: Uncomment if you want to see each individual spline
% prspline = spline(x, [prod(1,:)]);
% I2s = fncmb(fncmb(fnder(I_initial,2),delta),'+',spline(xx,ppval(prspline,xx).*ppval(I_initial,xx)));
% I2 = ppval(fnder(I_initial,2),x) * delta + prod(1,:).*ppval(I_initial,x);
% I2ss = spline(x,[0 0 I2(2:xes(2)) 0 0]);
% I2s = fncmb(I2s,'+',I_initial);
% 
% I2xx = fnder(I2s,2);
% I3 = delta * ppval(I2xx,x) + prod(2,:).*ppval(I2s,x);
% I3s = spline(x,[0 0 I3(2:xes(2)) 0 0]);
% I3s = fncmb(I3s,'+',I2s);
% 
% I3xx = fnder(I3s,2);
% I4 = delta * ppval(I3xx,x) + prod(3,:).*ppval(I3s,x);
% I4s = spline(x,[0 0 I4(2:xes(2)) 0 0]);
% I4s = fncmb(I4s,'+',I3s);
% 
% figure(4);
% plot(x,y,'o',xx,ppval(I_initial,xx),'-',xx,ppval(I2s,xx),'-',xx,ppval(I2ss,xx),'-',xx,ppval(I3s,xx),'-' ...
%     ,xx,ppval(I4s,xx),'-',xx,ppval(prspline,xx),'-');
% legend("points","I", "I2s","I2ss", "I3s", "I4s", "prspline");
