function [t,z,z1,g_s,amp_s,omega_s,q_s,fi_s] = Simulation(x0, tend, A,B,C,Z,P,offset)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global G

n=length(offset);

OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

[t,y] = ode45(@FHN,[0 tend],x0,OPTIONS,A,B,C,Z,G,P);

z = kron(eye(n),Z)*y';
z1 = z + fliplr(offset)';

[omega_s, T_s, amp_s, fi_s] = osciprof(z,t,n);

q_s = amp_s.*exp(1i*fi_s);
g_s = [fi_s;amp_s;omega_s];


end

