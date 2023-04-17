clear all
close all

%% Homework 3 Problem 3
% Author: Darrell Hamilton
% Date: 04/12/2023
% Description: This script is design to execute the work for Homework 3
% Question 3 of AE502. This script is to compare the analytic solutions of
% question 1 to those of question 2 for at least 20 different initial
% conditions with the same integration time of 100 time units and then make
% two different equinoctial elements plots h vs k and p vs q where. This is
% reproduced for w = {0.02, 0.1, 0.5}TU^-1

% Defining 20 different orbital initial conditions
a = 0.25:0.25:5;
e = 0.1:0.036:.8;
inc = 0:4:76;
ini_con = [a;e;inc];

t = 100;
mu = 1;
w = [0.01 0.02 0.1 0.5];
for mi = 1:length(w)

    for z = 1:size(ini_con,2)

        % Initial Conditions
        a = ini_con(1, z);
        e = ini_con(2, z);
        incl = ini_con(3, z)/180*pi;
        n = sqrt(mu/a^3);

        % Defining arbitrary initial conditions
        M = 0;
        Omega = 0;
        AOP = 0;
        l_0 = M;
        g_0 = AOP;
        h_0 = Omega;
        L_0 = n*a^2;
        G_0 = L_0 * sqrt(1 - e^2);
        H_0 = G_0 * cos(incl);

        %   Deluanay analytical expressions
        l_d = 1/L_0^3 * t + l_0;
        g_d = g_0;
        h_d = -w(mi)*t + h_0;
        L_d = L_0;
        G_d = G_0;
        H_d = H_0;

        % Canonical analytical expressions
        L_c = L_0 + w(mi)*H_0*L_0^3;
        l_c = 1/(L_c + w(mi)*H_0*L_c^3)^3 * t + l_0 - w(mi)*3*H_0*L_c^2*l_0;
        h_c = h_0 - w(mi).*L_c^3.*l_c;
        g_c = g_0;
        G_c = G_0;
        H_c = H_0;

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


        % Delaunay
        h_de(z) = e*sin(g_d + h_d);
        k_de(z) = e*cos(g_d + h_d);
        p_de(z) = tan(inc/2)*sin(h_d);
        q_de(z) = tan(inc/2)*cos(h_d);

        %  Canonical
        h_ce(z) = e*sin(h_c);
        k_ce(z) = e*cos(h_c);
        p_ce(z) = tan(inc/2)*sin(h_c);
        q_ce(z) = tan(inc/2)*cos(h_c);



    end
% Graphing equations
    figure('name',"w = " + string(w(mi)))
    subplot(2,1,1)
    scatter(k_ce, h_ce)
    hold on
    scatter(k_de, h_de)
    legend('Canonical','Delaunay');
    title('h vs k')

    hold off

    subplot(2,1,2)
    scatter(p_ce, q_ce)
    hold on
    scatter(p_de, q_de)
    legend('canonical','Delaunay');
    title('q vs p')
    hold off

end
