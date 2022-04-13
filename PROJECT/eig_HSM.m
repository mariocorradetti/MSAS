function [] = eig_HSM(data)
% Function that computes and plots the eigenvalues of the mechanical system
% INPUT
% data      structure containing useful datas
mech_data

% initial conditions
Ia0 = data.Ia0;
Ib0 = data.Ib0;
w0 = data.w0;
th0 = data.theta0;

R = motor.R;
L = motor.L;
b_EMF = motor.b_EMF;
tau = motor.tau;
Nr = motor.Nr;
B = motor.B;
J = motor.J;
J_2sp = 2*SP.J;
J_y = yoke.J;
Jtot = (J + (J_2sp + J_y)/tau);

% system matrix
A = [-R/L, 0, Nr*b_EMF*sin(Nr*th0)/L, Nr^2*b_EMF*w0*cos(Nr*th0)/L;
       0, -R/L,  -Nr*b_EMF*cos(Nr*th0)/L, Nr^2*b_EMF*w0*sin(Nr*th0)/L;
    -Nr*b_EMF*sin(Nr*th0)/Jtot, Nr*b_EMF*cos(Nr*th0)/Jtot,  -B/Jtot, - Ia0*b_EMF/Jtot*Nr^2*cos(Nr*th0) - Ib0*b_EMF/Jtot*Nr^2*sin(Nr*th0);
       0, 0, 1, 0];

plot(real(eig(A)),imag(eig(A)),'o');
grid on;
xlabel('Real(\lambda))'); ylabel('Imag(\lambda))')
title('Mechanical system eigenvalues for t = 0')
end