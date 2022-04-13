function[value,isterminal,direction]=eventfun(t,y,k,theta_step)
% The integration stops when theta reaches the closest multiple of theta_step
if y(4) - k*deg2rad(theta_step) < 1e-12
    value = y(4)-y(4);
    isterminal = 1;   % event stops the integration 
else
    value = y(4) - k*deg2rad(theta_step);
    isterminal = 0;   %  the integration does not stop
end

direction  = 0;
end