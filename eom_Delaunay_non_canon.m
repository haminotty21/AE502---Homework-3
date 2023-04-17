function [xout] = eom_Delaunay_non_canon(t,xin,var)
a = var(1);
e = var(2);
i = var(3);

w = 0.01;

L = xin(1);
G = xin(2);
H = xin(3);
l = xin(4);
g = xin(5);
h = xin(6);


xout = [0;0;0; ...
        1/(L^3); ...
        0; ...
        w];


end