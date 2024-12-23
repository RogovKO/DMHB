function dxdt = FHM_zeroeq(t,x0,H)
% FHN chapter.

global Az B C G gamma x_eq

[m,m]=size(G);
q = x0;

y =  kron(eye(m),C)*q;
ylag = kron(eye(m),C)*(H(:,1));
nl = (-1/3)*y.^3 - x_eq(2)*y.^2;
u = (-gamma*eye(m)*y + gamma*G*ylag);

dxdt = kron(eye(m),Az)*q + kron(eye(m),B)*nl + kron(eye(m),B)*u;
end
