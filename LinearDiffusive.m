% Driving Parameters
numxpoints = 101;
xes = [0 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
% Initial condition from paper i.e. \phi(x)
y = [0 6 2 2 1 1 1 0.5 0];
xx = linspace(xes(1),xes(2),numxpoints);

%%% h(x)
%TODO: Use Mary's distance data to fit a more accurate h(x) function
rho = -0.9478;
sigma = 8.9149;
h = -(x - rho).*(x - sigma);

%%% r(t)
Beta = 0.0019;
alpha = 1.5526;
gamma = 0.078;
delta = 0.002;
asymptote = Beta / alpha;
r = asymptote - exp(-alpha*(t - 1))*(asymptote - gamma);

prod = r'*h;
% Show growth mesh
figure(1);
[X_p, Y_p] = meshgrid(x,t);
mesh(X_p,Y_p,prod,'FaceAlpha','0.5','FaceColor','flat')


% The key to this implementation is that you can easily take the derivative
% of a spline using `fnder()`
I_initial = spline(x,[0 y 0]);

% Mesh stuff
[X, Y] = meshgrid(xx,t);
Z = zeros(ts(2), numxpoints);
xlabel("x Distance"); ylabel("t Time"); zlabel("z Growth");
title("Growth Function");


Ispline = I_initial;
for t_i = t
    Z(t_i,:) = ppval(Ispline,xx);
%     I_t = delta * ppval(fnder(Ispline,2),x) + prod(t_i,:).*ppval(Ispline,x);
%     I_t = spline(x,[0 0 I_t(2:xes(2)) 0 0]);
    I_t = fncmb(fncmb(fnder(Ispline,2),delta),'+', ...
        spline(xx,ppval(spline(x, prod(t_i,:)),xx).*ppval(Ispline,xx)));
    Ispline = fncmb(I_t,'+',Ispline);
end

figure(2);
mesh(X,Y,Z,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");
title("Linear Diffusive Information Diffusion");
zlim([0 100])


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
% figure(3);
% plot(x,y,'o',xx,ppval(I_initial,xx),'-',xx,ppval(I2s,xx),'-',xx,ppval(I2ss,xx),'-',xx,ppval(I3s,xx),'-' ...
%     ,xx,ppval(I4s,xx),'-',xx,ppval(prspline,xx),'-');
% legend("points","I", "I2s","I2ss", "I3s", "I4s", "prspline");

