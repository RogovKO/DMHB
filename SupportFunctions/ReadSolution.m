function [g,q,omega,amp,fi] = ReadSolution(g,k)
% Read the solution to the problem
%   Detailed explanation goes here
fi = g(1:(k-1)/2);
amp = g(((k-1)/2+1):k-1);
omega = g(k);
q = amp.*exp(1i*fi);

T = 2*pi/omega;
fir= rad2deg(fi);


di1 = ['Frequency is equal to ',num2str(omega)];
di2 = ['Period is equal to ',num2str(T)];
disp(di1)
disp(di2)

disp('Amplitudes are equal to ');
disp(amp)

disp('The phase shift is');
disp(fir)
 

end

