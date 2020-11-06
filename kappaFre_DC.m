function [ei, kappa] = kappaFre_DC(A,B,C,lam,gamma) 

% Looks for k and omega for MHB EA case

ref = 3;
K = 0;
kappa = -0.5;
fre=[];
Kap_max = 2;

while ref > 10e-3
    
    ei = eig(A - B*C*(-kappa + gamma - gamma*lam));
    rei = real(ei);
    
    for i = 1:length(rei)
        if abs(rei(i))<1e-6 
            f=real(1i*ei(i));
            fre=[fre;f];           
        end
    end
  if ~isempty(fre)
     break
  end
    kappa = kappa + 0.000001;
 if kappa > Kap_max
     break
  end
end
end