function [AA, Res] = MHBequation_shifted(X,n)
% MHB equation for Lure systems
%   Detailed explanation goes here

global W G gamma tau 

n = length(X);
I = eye((n-1)/3);

fi = X(1:(n-1)/3);
a = X(((n-1)/3+1):2*(n-1)/3);
b = X(2*(n-1)/3+1:n-1);
w = X(end);

%% Output

q =  a.*exp(1i*fi);
y0 = b;
y = y0 + q;

%% Describing functions

K0 = 49*b/50 - b.^2/3 - a.^2/2 + 49*a.^2./(100*b);
K1 = - a.^2/4 - b.^2 + 49*b/25;

%% Parameters 

c0 = (1/double(W(0)))*I;
c1 = (1/double(W(1i*w)))*I;
de = exp(-1i*w*tau);
M0 = gamma*G - gamma*I;
M1 = gamma*de*G - gamma*I;

%% MHB Equation

%  MHB = C0 - C1*y0 + (C1 - c1 + M)*y;

MHB1 = (c0 - M0 - diag(K0))*y0;
MHB2 = (c1 - M1 - diag(K1))*q;

AA = (c0 - M0 - K0)+(c1 - M1 - K1);

Res = MHB1+MHB2;
end

