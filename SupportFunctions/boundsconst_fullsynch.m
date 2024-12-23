function [lbound,ubound] = boundsconst_fullsynch(x,fre,Beta)
%   Lower and upper bounds constarint
%   Detailed explanation goes here
n = length(x);

lbound = [zeros((n-1)/3,1); -10e-4*ones(((n-1)/3),1);-10e-4*ones(((n-1)/3),1); 0.8*fre];
ubound = [2*pi*ones(((n-1)/3),1); Beta*ones(((n-1)/3),1);Beta*ones(((n-1)/3),1); 1.5*fre];
end

