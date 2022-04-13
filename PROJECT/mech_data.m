% Hybrid stepper motor parameters

SP.L = 325e-3;      % solar panel length [m]
SP.h = 205e-3;      % solar panel heightBB [m]
SP.t = 5e-3;        % solar panel thickness [m]
SP.m = 0.5;         % solar panel mass [kg]
SP.J = 1/12*SP.m*SP.h^2;        % inertia of one panel [Kg m^2]

motor.step = 20;    % step per revolution
motor.R = 13.6;     % terminal resistance [Ohm] 
motor.L = 2e-3;     % terminal inductance [H]
motor.b_EMF = 0.53; % back-EMF amplitude [Vs]
motor.tau  = 4;     % gearhead reduction [-]
motor.V = 12;       % voltage [V]
motor.theta_step = 360/motor.step;  % theta step [deg]
motor.phases = 4;   % full step 
motor.Nr = 360/(motor.theta_step*motor.phases); % number of rotor's teeth 
motor.B = 8e-5;     % viscous friction [Nms/rad]
motor.T_L = 0;      % torque load
motor.J = 4.5e-5;   % rotor inertia [Kg m^2] 

yoke.rho = 2700;                            % density Al 6061 T6 [Kg/m^3] 
yoke.t = 5e-3;                              % thickness [m]
yoke.L =  198e-3;                           % length [m]
yoke.h = 1e-1;                              % height [m]
yoke.m = yoke.rho*yoke.L*yoke.h*yoke.t;     % mass [kg]
yoke.J = 1/12*yoke.m*(yoke.h^2);            % inertia [Kg m^2]
