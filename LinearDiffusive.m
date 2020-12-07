% Driving Parameters
numxpoints = 101;
xes = [0 6];
ts = [1 6];

x = xes(1):xes(2);
t = ts(1):ts(2);
% Initial condition from paper i.e. \phi(x)
y = [0 18 7 8 5 4 0];
xx = linspace(xes(1),xes(2),numxpoints);

%%% h(x)
%TODO: Use Mary's distance data to fit a more accurate h(x) function
rho = -0.95;
sigma = 8.9;
h = -(x - rho).*(x - sigma);

%%% r(t)
Beta = 0.0059;
alpha = 1.55;
gamma = 0.078;
delta = 0.002;
asymptote = Beta / alpha;
r = asymptote - exp(-alpha*(t - 1))*(asymptote - gamma);

prod = r'*h;



% The key to this implementation is that you can easily take the derivative
% of a spline using `fnder()`
I_initial = spline(x,[0 y 0]);

% Mesh stuff
[X, Y] = meshgrid(xx,t);
Z = zeros(ts(2), numxpoints);

Ispline = I_initial;
for t_i = t
    Z(t_i,:) = ppval(Ispline,xx);
    I = delta * ppval(fnder(Ispline,2),x) + prod(t_i,:).*ppval(Ispline,x);
    Ispline = spline(x,[0 0 I(2:6) 0 0]);
end

mesh(X,Y,Z,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");
title("Linear Diffusive Information Diffusion");


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

