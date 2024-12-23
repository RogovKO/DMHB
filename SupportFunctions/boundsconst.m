function [lbound,ubound] = boundsconst(x,fre,Beta)
%Lower and upper bounds constarint
%   Detailed explanation goes here
n = length(x);

lbound = [zeros((n-1)/3,1); zeros(((n-1)/3),1); zeros(((n-1)/3),1); 0.9*fre];
ubound = [0; 2*pi*ones(((n-1)/3)-1,1);  Beta*ones(((n-1)/3),1); Beta*ones(((n-1)/3),1); 1.1*fre];
end

