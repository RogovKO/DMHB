function [time,qsim] = Oscillations(amp, omega,ifi,fi,offset,tend)
% Build approximation of oscillations
%   Detailed explanation goes here

qsim = [];
time = 0.0:0.1:tend;

for tau = 0.0:0.1:tend
     [qsim_step] =  amp.*sin(omega*tau + ifi + fi) + fliplr(offset)';               % phases must be in RAD
      qsim = [qsim, qsim_step];
end


end

