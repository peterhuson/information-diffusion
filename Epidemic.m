
% Driving Parameters
numxpoints = 8;
xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);
xx = linspace(xes(1),xes(2),numxpoints);

% Initial condition
y_inf = [1.350662534142069 1.241994141180107 1.209190322169032 1.203465346130819 1.202576116353606 1.199396102035386 1.083287067674959 1.0];
y_sus = exp(1 - log(y_inf));

N = sum(y_sus);

%%% beta (transmission coefficient)
beta = 1.65;

%%% diffusion constant
delta_I = 0.02;
delta_S = 0.005;


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
    % Update and replace NaN
    update_S = beta * ppval(Ispline,x).*ppval(Sspline,x)./ N;
    update_S(isnan(update_S)) = 0;
    
    S_t = fncmb(fncmb(fnder(Sspline,2),delta_S),'-', ...
        spline(x, update_S));
    Sspline = fncmb(S_t,'+',Sspline);
    
    Z_inf(t_i,:) = ppval(Ispline, xx);
    % Update and replace NaN
    update_I = beta * ppval(Ispline,x).*ppval(Sspline_orig,x)./N;
    update_I(isnan(update_I)) = 0;
    
    I_t = fncmb(fncmb(fnder(Ispline,2),delta_I),'+', ...
        spline(x, update_I));
    Ispline = fncmb(I_t,'+',Ispline);
end

fnplt(Sspline)

% Normalize density between 0 and 1
norm_Z_sus = Z_sus - min(Z_sus(:));
norm_Z_sus = norm_Z_sus./ max(norm_Z_sus(:));

norm_Z_inf = Z_inf - min(Z_inf(:));
norm_Z_inf = norm_Z_inf./ max(norm_Z_inf(:));

% Graph the susceptible and infected
mesh(X,Y,Z_sus,'FaceAlpha','0.5','FaceColor','green')
hold on
mesh(X,Y,Z_inf,'FaceAlpha','0.5','FaceColor','red')
hold on

Z = readmatrix('Accuracy.txt');
X_data = readmatrix('X.txt');
Y_data = readmatrix('Y.txt');

%Graphing the data from story 714
mesh(X_data,Y_data,Z,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");

title('\fontsize{14} SI Model Information Diffusion');

lgd = legend('Susceptible', 'Infected', 'Plot of data for story 714');
lgd.FontSize = 14;

%calculating relative error between Digg data and our SI model 
relError(Z_inf, Z);