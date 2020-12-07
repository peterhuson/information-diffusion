% Driving Parameters
numxpoints = 101;
xes = [0 6];
ts = [1 6];

x = xes(1):xes(2);
t = ts(1):ts(2);

% Initial condition
%%% Susceptible starts out with All people people (density of 1)
%%% Infected starts out with 1 person (density of 1/N)
N = 1000; % Just an arbitrary population size.
y_sus = [0 N 0 0 0 0 0];
y_inf = [0 1 0 0 0 0 0];
xx = linspace(xes(1),xes(2),numxpoints);

%%% beta (transmission coefficient)
beta = 3; %just a guess.. need to know average links between people in the network

%%% diffusion constant? also just random
delta = 0.005;


% The key to this implementation is that you can easily take the derivative
% of a spline using `fnder()`
S_initial = spline(x, [0 y_sus 0]);
I_initial = spline(x, [0 y_inf 0]);

% Mesh stuff
[X, Y] = meshgrid(xx,t);
Z_sus = zeros(ts(2), numxpoints);
Z_inf = zeros(ts(2), numxpoints);

Sspline = S_initial;
Ispline = I_initial;
for t_i = t
    Sspline_orig = Sspline; % it gets overwritten so you need to save a copy of the original.
    
    Z_sus(t_i,:) = ppval(Sspline, xx);
    S = delta * ppval(fnder(Sspline,2),x) - beta * ppval(Ispline,x).*ppval(Sspline,x)/(ppval(Ispline,x) + ppval(Sspline,x));
    Sspline = spline(x,[0 0 S(2:6) 0 0]);
    
    Z_inf(t_i,:) = ppval(Ispline, xx);
    I = delta * ppval(fnder(Ispline,2),x) + beta * ppval(Ispline,x).*ppval(Sspline_orig,x) /(ppval(Ispline,x) + ppval(Sspline_orig,x));
    Ispline = spline(x,[0 0 I(2:6) 0 0]);
end

% Normalize density between 0 and 1 (plz work)
norm_Z_sus = Z_sus - min(Z_sus(:));
norm_Z_sus = norm_Z_sus./ max(norm_Z_sus(:));

norm_Z_inf = Z_inf - min(Z_inf(:));
norm_Z_inf = norm_Z_inf./ max(norm_Z_inf(:));

% Graph the susceptible and infected
mesh(X,Y,norm_Z_sus,'FaceAlpha','0.5','FaceColor','flat')
hold on
mesh(X,Y,norm_Z_inf,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");
title("Susceptible and Infected Populations");


%%% Debugging Stuff: Uncomment if you want to see each individual spline
% I2 = delta * ppval(fnder(I_initial,2),x) + prod(1,:).*ppval(I_initial,x);
% I2s = spline(x,[0 0 I2(2:6) 0 0]);
% 
% I2xx = fnder(I2s,2);
% I3 = delta * ppval(I2xx,x) + prod(2,:).*ppval(I2s,x);
% I3s = spline(x,[0 0 I3(2:6) 0 0]);
% 
% I3xx = fnder(I3s,2);
% I4 = delta * ppval(I3xx,x) + prod(3,:).*ppval(I3s,x);
% I4s = spline(x,[0 0 I4(2:6) 0 0]);
% 
% figure(2);
% plot(x,y,'o',xx,ppval(I_initial,xx),'-',xx,ppval(I2s,xx),'-',xx,ppval(I3s,xx),'-' ...
%     ,xx,ppval(I4s,xx),'-');
% legend("points","I", "I2s", "I3s", "I4s");

