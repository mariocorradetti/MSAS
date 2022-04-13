function [dy,par] = HSMeq(t,y,data,motor,SP,yoke)
% Function that describes the dynamics of the HSM
% INPUTS 
% t             time [s]
% y             vector [Ia Ib w th]
% data, motor, SP, yoke:    structures containing system parameters
%  
% OUTPUTS
% dy            derivatives of the variables [Ia_d Ib_d w_d th_d]
% par           vector of useful outputs [Va Vb Torque]

% mech_data;
% Initialize parameters
Ia = y(1);
Ib = y(2);
w = y(3);
th = y(4);

R = motor.R;
L = motor.L;
Nr = motor.Nr;
bEMF = motor.b_EMF;
B = motor.B;
J = motor.J;
J_sp = SP.J*2;
J_y = yoke.J;
tau = motor.tau;


% Define the voltages of HSM
if data.index == 4 
    Va=motor.V;
    Vb=0;
elseif data.index == 1 
    Va=0;
    Vb=motor.V;
elseif data.index == 2 
    Va=-motor.V;
    Vb=0;
else
    Va=0;
    Vb=-motor.V;
end
   
% Equations
dy = zeros(4,1);
% Ia_dot
dy(1) = (Va - R*Ia + w*Nr*bEMF*sin(Nr*th) )/L; 
% Ib_dot
dy(2) = (Vb - R*Ib - w*Nr*bEMF*cos(Nr*th) )/L; 
% w_dot
dy(3) = (-Ia*Nr*bEMF*sin(Nr*th) + Ib*Nr*bEMF*cos(Nr*th) - B*w)/(J + (J_sp+J_y)/tau);     
% theta_dot
dy(4) = w;

Torque = -Ia*Nr*bEMF*sin(Nr*th)+Ib*Nr*bEMF*cos(Nr*th)-B*w;

par = [Va Vb Torque];
end