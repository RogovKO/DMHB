function [time,qsim] = Oscillations_shifted(amp,B,omega,ifi,fi,offset,tend)
% Build approximation of oscillations
%   Detailed explanation goes here

qsim = [];
time = 0.0:0.1:tend;

for tau = 0.0:0.1:tend
     [qsim_step] = B + amp.*sin(omega*tau + ifi + fi) + fliplr(offset)';               % phases must be in RAD
      qsim = [qsim, qsim_step];
end


end

