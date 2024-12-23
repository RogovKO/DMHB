function dxdt = ddex3de(t,x0,H)
% FHN chapter.

global A B C Z G P gamma

[m,m]=size(G);
x = x0;

y =  kron(eye(m),C)*x;
ylag = kron(eye(m),C)*H(:,1);
z =  kron(eye(m),Z)*x.^3;
u = (-gamma*eye(m)*y + G*ylag);

dxdt = kron(eye(m),A)*x + kron(ones(m,1),P) - kron(eye(m),B)*z + kron(eye(m),B)*u;

end

