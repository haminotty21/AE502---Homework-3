function [xout] = eom_Delaunay_canonical(t,xin,var)
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

% Canonical Variable transformations
l = l - 3*w*H*L^2*l;
h = h - w*L^3*l;

% L = L - w*H*L^3;

xout = [0;0;0; ...
        1/(L + w*H*L^3)^3 *(1 + 3*w*H*L^2); ...
        0; ...
        1/(L + w*H*L^3)^3 * (w*L^3)];


end