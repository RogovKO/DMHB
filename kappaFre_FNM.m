function [fre, kappa] = kappaFre_FNM(A,B,C,lam,tau,om,gamma) 

% Looks for k and omega for MHB EA case

ref = 3;
K = 0;
kappa = -0.5;
fre=[];
Kap_max = 3;


while ref > 10e-4
    
    ei = eig(A - B*C*(-kappa + gamma - gamma*exp(-1i*om*tau)*lam));
    rei = real(ei);
    
    for i = 1:length(rei)
        if abs(rei(i))<1e-6 
            f=real(1i*ei(i));
            fre=[fre;f];           
        end
    end
    kappa = kappa + 0.000001;
  if ~isempty(fre)
     break
  end
  if kappa > Kap_max
     break
  end
  
end

ei
end