function [fre, kappa] = kappaFre_FNM(A,B,C,Z,lam,tau,om,gamma) 

% Looks for k and omega for MHB EA case

ref = 3;
K = 0;
kappa = 0.00001;
fre=[];

while ref > 10e-4
    
    ei = eig(A - B*(kappa*Z + gamma*C - gamma*exp(-1i*om*tau)*lam*C));
    rei = real(ei);
    
    for i = 1:length(rei)
        if abs(rei(i))<1e-6 
            f=real(1i*ei(i));
            fre=[fre;f];           
        end
    end
    kappa = kappa + 0.000001;
  if length(fre) > 1
     break
  end
end
end