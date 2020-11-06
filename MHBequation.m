function [AA, Res] = MHBequation(g,k)
% MHB equation for Lure systems
%   Detailed explanation goes here

global G W gamma tau

fi = g(1:(k-1)/2);
alf = g(((k-1)/2+1):k-1);
w = g(k);
q = alf.*exp(1i*fi);
I = eye((k-1)/2);

c1 = double(1/W(1i*w));
c2 = exp(-1i*w*tau);
k = (1/2)*abs(alf).^2;


AA = diag(k) + c1*I + gamma*(I - c2*G);
Res = AA*q;
end

