function F = opti(X)
% Creates objective function that we need to solve optimization problem
%   Detailed explanation goes here

global Wz Wy G

n = length(X);

fi = X(1:(n-1)/2);
alf = X(((n-1)/2+1):n-1);
w = X(n);
q = alf.*exp(1i*fi);

c = double(Wz(1i*w)/Wy(1i*w));
c1 = double((1/Wy(1i*w)));

AA = c.*diag(3*abs(alf).^2/4) + c1*eye((n-1)/2) + G;

z = AA*q;
k = alf-abs(q);

F = z'*z + (k'*k);

end