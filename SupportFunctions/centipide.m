function dxdt = centipide(t,x0,H)
% Centipide paper

global A B C Z G gamma 

[m,m]=size(G);
q = x0;

y =  kron(eye(m),C)*q;
z =  kron(eye(m),Z)*q;
ylag = kron(eye(m),C)*(H(:,1));
nl = z.^3;
u = (-gamma*eye(m)*y + gamma*G*ylag);

dxdt = kron(eye(m),A)*q - kron(eye(m),B)*nl + kron(eye(m),B)*u;
end

