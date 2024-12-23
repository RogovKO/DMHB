function F = costfunctiondelay(X)
% Creates objective function that we need to solve optimization problem
%   Detailed explanation goes here

global Wz Wy G gamma tau x_eq

n = length(X);

fi = X(1:(n-1)/2);
alf = X(((n-1)/2+1):n-1);
f = X(n);
q = alf.*exp(1i*fi);

c = double(Wz(1i*f)/Wy(1i*f));
c1 = double((1/Wy(1i*f)));
c2 = exp(-1i*f*tau);

AA = c*diag(((3/4)*abs(alf).^2 - (147/200)*abs(alf))) + c1*eye((n-1)/2) + gamma*eye(length(G(n-1)/2)) - c2*gamma*G;

z = AA*q;
k = alf-abs(q);

F = z'*z + (k'*k);

end