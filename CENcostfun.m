function F = CENcostfun(X)
%   Creates objective function that we need to solve optimization problem
%   Detailed explanation goes here

global G gamma tau A B C Z

n = length(X);
I = eye((n-1)/3);

fi = X(1:(n-1)/3);
a = X(((n-1)/3+1):2*(n-1)/3);
b = X(2*(n-1)/3+1:n-1);
w = X(end);


Wz = Z*((1i*w*eye(length(A)) - A)^-1)*B;
Wz0 = Z*((0*eye(length(A)) - A)^-1)*B;

Wy = C*((1i*w*eye(length(A)) - A)^-1)*B;
Wy0 = C*((0*eye(length(A)) - A)^-1)*B;


%% Output

q =  a.*exp(1i*fi);
y0 = b;
y = y0 + q;

%% Describing functions

K0 = b.^2 + 3*a.^2/2;
K1 = 3*b.^2 + 3*a.^2/4;

%% Parameters 

c0 = (1/Wy0)*I;
c02 = (Wz0/Wy0)*I;
c1 = (1/Wy)*I;
c12 = (Wz/Wy)*I;
de = exp(-1i*w*tau);
M0 = gamma*G - gamma*I;
M1 = gamma*de*G - gamma*I;

%% MHB Equation

%  MHB = C0 - C1*y0 + (C1 - c1 + M)*y;

MHB1 = (c0 - M0 + c02*diag(K0))*y0;
MHB2 = (c1 - M1 + c12*diag(K1))*q;

%% Objective function

k1 = (a) - abs(q);

F = MHB1'*MHB1 + MHB2'*MHB2 +(k1'*k1);

end