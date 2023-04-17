clear all
close all
%% Homework 3 Problem 2
% Author: Darrell Hamilton
% Date: 04/12/2023
% Description: 


% Initial Conditions
a = 1;
e = 0.5;
incl = 45/180*pi;
t = 1:0.1:100;
mu = 1;
n = sqrt(mu/a^3);

% Starting orbit at a "zero-ed out" position
M = 0;
AOP = 0;
LAN = 0;

% Defining initial conditions
l_0 = M;
g_0 = AOP;
h_0 = LAN;
L_0 = n*a^2;
G_0 = L_0 * sqrt(1 - e^2);
H_0 = G_0 * cos(incl);

% Using ode45 to integrate over 100 time units
x_0 = [L_0 G_0 H_0 l_0 g_0 h_0];
options = odeset('AbsTol',1e-13,'RelTol',1e-13);
[t,X] = ode45(@eom_Delaunay_non_canon, t, x_0, options, [a, e, incl]);

% Renaming output
M = X(:,4)';
AOP = X(:,5)';
LAN = X(:,6)';

% Solve for true anomaly by solving keplar's equation
for i = 1:length(M)
    if M(i) <= pi
        E_0 = M(i) + e/2;
    elseif M(i) > pi
        E_0 = M(i) - e/2;
    end
    
    %   Solve Keplar's equation to calculate True Anomaly (theta)
    E(1) = E_0;
    thres = 10^-8;
    FE = thres + 1;
    inc = 1;
    while(abs(FE) > thres)
        FE = (E(inc) - e*sin(E(inc))-M(i))/(1-e*cos(E(inc)));
        E(inc+1) = E(inc) - FE;
        inc = inc + 1;
    end
    E = E(end);
    f(i) = 2*atan2(sqrt(1+e)*tan(0.5*E),sqrt(1-e));
end


% Defining theta for rotation matrix
theta = AOP + f;
r = a*(1-e^2)./(1+e*cos(f));
incl(1:length(theta)) = incl;

% Calculating position
for i = 1:length(theta)
    r_vec(i,1:3) = r(i)*[(cos(theta(i))*cos(LAN(i)) - cos(incl(i))*sin(LAN(i))*sin(theta(i))), ...
        cos(theta(i))*sin(LAN(i)) + cos(incl(i))*cos(LAN(i))*sin(theta(i)), ...
        sin(incl(i))*sin(theta(i))];
end

% Making graphs
[x,y,z] = sphere(20);
figure
plot3(r_vec(:,1),r_vec(:,2),r_vec(:,3))
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
hold on
axis equal
surf(x/5, y/5, z/5);
hold off

figure
subplot(3,1,1)
plot(t,f)

subplot(3,1,2)
plot(t,X(:,5))

subplot(3,1,3)
plot(t,X(:,6))
