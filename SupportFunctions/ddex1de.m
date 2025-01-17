function dydt = ddex1de(t,y,Z)
% Differential equations function for DDEX1.
ylag1 = Z(:,1);
ylag2 = Z(:,1);
dydt = [ ylag1(1)
   ylag1(1) + ylag2(2)
   y(2)               ];
end

