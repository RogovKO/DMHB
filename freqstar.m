function [fre] = freqstar(L)
% Looks for bifurcation point of the network
%   Detailed explanation goes here

ref = -10;
K = 0;
w = 0.00001;
count = 0;

fre=[];

while ref < 10e-6
    
   Li = double(L(w))
    
   ref = real(Li)
   
   f=w;
   
   fre=[fre;f];           

   w = w + 0.000001;
   count = count + 1;
 
end

end
