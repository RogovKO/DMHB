function s = ICnearzero(t)
    % Constant history function for DDEX1.
    global nn
    s = 1e-4*rand(2*nn,1);
end
