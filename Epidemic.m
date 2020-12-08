% Driving Parameters
numxpoints = 101;
xes = [1 8];
ts = [1 15];

x = xes(1):xes(2);
t = ts(1):ts(2);

% Initial condition
%%% Susceptible starts out with all people (density of 1)
%%% Infected starts out with 1 person (density of 1/N)
N = 1000; % Just an arbitrary population size.
%y_sus = [0 N-1 N-4 N-94 N-171 N-181 N-183 N-183 N-184]; % one hour into the spread; create some sort of smooth dist
%y_inf = [0 1 4 94 171 181 183 183 184];

y_sus = [ 1 4 94 171 181 183 183 183];
%y_sus = [0 1 5 99 270 451 634 817 1000];
y_inf = [ 1 1 1 1 1 1 1 1];
xx = linspace(xes(1),xes(2),numxpoints);


%%% beta (transmission coefficient)
beta = 5.14; %0.000008; %just a guess.. need to know average links between people in the network

%%% diffusion constant
delta_I = 0.01;
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

%fnplt(Sspline)
%fnplt(Ispline)
for t_i = t
    Sspline_orig = Sspline; % it gets overwritten so you need to save a copy of the original.
    
    %I = delta * ppval(fnder(Ispline,2),x) + prod(t_i,:).*ppval(Ispline,x);
    %Ispline = spline(x,[0 0 I(2:6) 0 0]);
    %        --> 
    %I_t = fncmb(fncmb(fnder(Ispline,2),delta),'+', ... spline(xx,ppval(spline(x, prod(t_i,:)),xx).*ppval(Ispline,xx)));
    %Ispline = fncmb(I_t,'+',Ispline);
    
    Z_sus(t_i,:) = ppval(Sspline, xx);
    % Update and replace NaN
    update_S = beta * ppval(Ispline,x).*ppval(Sspline,x)./N;
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

% Normalize density between 0 and 1 (plz work)
norm_Z_sus = Z_sus - min(Z_sus(:));
norm_Z_sus = norm_Z_sus./ max(norm_Z_sus(:));

norm_Z_inf = Z_inf - min(Z_inf(:));
norm_Z_inf = norm_Z_inf./ max(norm_Z_inf(:));

% Graph the susceptible and infected
mesh(X,Y,Z_sus,'FaceAlpha','0.5','FaceColor','green')
hold on
mesh(X,Y,Z_inf,'FaceAlpha','0.5','FaceColor','red')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");
title("Susceptible and Infected Populations");
