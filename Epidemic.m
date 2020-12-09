% import ('relError.m')

% populations = readmatrix("digg2009/votes_714.csv");


populations = [1557 346744 1641109 1937290 2003317 2018062 2021615 2022430];
populations = [populations(1) diff(populations)];
% Driving Parameters
numxpoints = 8;
xes = [1 8];
ts = [1 50];

x = xes(1):xes(2);
t = ts(1):ts(2);

% Initial condition
%%% Susceptible starts out with all people (density of 1)
%%% Infected starts out with 1 person (density of 1/N)
% N = 1000; % Just an arbitrary population size.
% N = 24099;
%y_sus = [0 N-1 N-4 N-94 N-171 N-181 N-183 N-183 N-184]; % one hour into the spread; create some sort of smooth dist
%y_inf = [0 1 4 94 171 181 183 183 184];

% y_sus = [ 1 4 94 171 181 183 183 183]; % old

% y_sus = [1557 346744 1641109 1937290 2003317 2018062 2021615 2022430]
% y_sus = [populations(1) diff(populations)];
% y_sus = [y_sus(1:end-1)];

y_sus = [0.93513166 0.99838498 0.99933886 0.99982295 0.99996206 0.99999108 0.99999901 1];
% y_sus = y_sus .* 10000
y_sus = exp(y_sus)
% y_inf = [648.6833654464 16.1502434072 6.6113829124 1.7705144816 0.3793708135 0.0891944846 0.0098930805 0.0]

y_inf = [1.350662534142069 1.241994141180107 1.209190322169032 1.203465346130819 1.202576116353606 1.199396102035386 1.083287067674959 1.0];

y_sus = exp(1 - log(y_inf));
% y_inf = [1.06701853 1.00161633 1.00066136 1.00017707 1.00003794 1.00000892 1.00000099 1];


N = sum(y_sus);


%y_sus = [0 1 5 99 270 451 634 817 1000];
% y_inf = [ 1 20  10 10 10 10 10 10]; %old 


xx = linspace(xes(1),xes(2),numxpoints);


%%% beta (transmission coefficient)
beta = 0.4; %0.000008; %just a guess.. need to know average links between people in the network

%%% diffusion constant
delta_I = 0.001;
delta_S = 0.001;


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

% Normalize density between 0 and 1
norm_Z_sus = Z_sus - min(Z_sus(:));
norm_Z_sus = norm_Z_sus./ max(norm_Z_sus(:));

norm_Z_inf = Z_inf - min(Z_inf(:));
norm_Z_inf = norm_Z_inf./ max(norm_Z_inf(:));

Z_sus
% Z_sus = Z_sus ./ populations

% Graph the susceptible and infected
mesh(X,Y,Z_sus,'FaceAlpha','0.5','FaceColor','green')
hold on
mesh(X,Y,Z_inf,'FaceAlpha','0.5','FaceColor','red')
hold on

Z = readmatrix('Accuracy.txt');
X_data = readmatrix('X.txt');
Y_data = readmatrix('Y.txt');

mesh(X_data,Y_data,Z,'FaceAlpha','0.5','FaceColor','flat')
xlabel("x Distance");
ylabel("t Time");
zlabel("z Density");

title('\fontsize{14}Susceptible and Infected Populations');

lgd = legend('Susceptible', 'Infected');
lgd.FontSize = 14;



Error = 0;

Z(1,2)
%|approx - exact| / exact 

% for t = 2:length(Z)
%     for d = 1:8
% %         abs(norm_Z_inf(t,d) - Z(t,d))
%         Error = Error + abs(norm_Z_inf(t,d) - Z(t,d))/Z(t,d)
%     end
% end
% 
% Error

relError(Z_inf, Z)
% 
% deltaSignal = abs(norm_Z_inf - Z);
% percentageDifference = norm_Z_inf ./ Z; % Percent by element.
% meanPctDiff = mean(percentageDifference); % Average percentage over all elements.