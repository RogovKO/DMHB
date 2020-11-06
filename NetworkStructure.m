function [G,n,offset] = NetworkStructure(nn,flag)
% Set up network structure
%   Detailed explanation goes here

if flag==1

    I = eye(nn);
    Gring = I - I([2:end 1],:);
    Gring = Gring + Gring';
    G = [I+Gring -I; -I I];
    n = length(G);
    offset = 1:1:(n);
    figure
    plot(digraph(G-diag(diag(G))))
    title('Network Structure')
    axis on
    
elseif flag==0
  
    I = eye(nn);
    Gring = I - I([2:end 1],:);
    G = Gring + Gring';  
    n=nn;
    offset = 0:1:(n-1);
    figure
    plot(digraph(G-diag(diag(G))))
    title('Network Structure')
    axis on
    
elseif flag==2
  
    I = eye(nn);
    Gring = I([2:end 1],:);
    G = Gring; 
    n=nn;
    offset = 0:1:(n-1);
    figure
    plot(digraph(G))
    title('Network Structure')
    axis on
else
    
di3 ='Wrong flag: choose 0,1 or 2';
disp(di3)
   
end

end

