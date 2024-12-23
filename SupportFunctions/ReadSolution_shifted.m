function [g,q,omega,amp,fi,B] = ReadSolution_shifted(g,n)
% Read the solution to the problem
%   Detailed explanation goes here

% fi = g(1:(k-2)/2);
% amp = g(((k-2)/2+1):k-2);
% B = g(k-1);
% omega = g(end);


fi = g(1:(n-1)/3);
amp = g(((n-1)/3+1):2*(n-1)/3);
B = g(2*(n-1)/3+1:n-1);
omega = g(end);

q = amp.*exp(1i*fi);
T = 2*pi/omega;
fir= rad2deg(fi);


di1 = ['Frequency is equal to ',num2str(omega)];
di2 = ['Period is equal to ',num2str(T)];
% di3 = ['DC is ',num2str(B)];

disp(di1)
disp(di2)
% disp(di3)

disp('DC is ');
disp(B)

disp('Amplitudes are equal to ');
disp(amp)

disp('The phase shift is (RAD)');
disp(fi)
 
disp('The phase shift is (DEG)');
disp(fir)
 
end

