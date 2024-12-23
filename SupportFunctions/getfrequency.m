function [freq] = getfrequency(A,B,OM)
%   Looks for bifurcation point of the network
%   Detailed explanation goes here

freq=[];

for ii = 1:length(OM)
  
    ei = eig(A-B*OM(ii));
    rei = real(ei);
    
    for i = 1:length(rei)
        if abs(rei(i))<1e-6  
            f=real(1i*ei(i));
            freq=[freq;f];           
        end
    end
    
end   

end

