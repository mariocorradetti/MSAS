function [] = eig_HEAT(T0)

% Function that computes the eigenvalues for the thermal system at initial
% conditions
% 
% INPUT
%   T0: vectors of temperature at time t=0 (initial conditions)
% OUTPUT
%   Plot of the eigenvalues
% 

%load of the data
mech_data
therm_data

%initialization
sigma = 5.67037e-8;
A_sq = body.h^2;                    
A_rect = body.h*body.L;             
A_p = SP.L*SP.h;       
m_sq = AL.rho*body.t*body.h^2;          
m_rect = AL.rho*body.t*body.h*body.L;  
m_p = AL.rho*SP.t*SP.L*SP.h;

A = [A_sq A_rect A_rect A_sq A_rect A_rect A_p A_p A_p A_p];        %areas vector
eps_sp = (SC.eps*SC.A + AL.eps*(A_p-SC.A) + back.eps*A_p);          %equivalent emissivity solar panels  
eps_coat = (AL.eps*48 + AU.eps*27 + AG.eps*25)/100;                 %body emissivity
C = [m_sq m_rect m_rect m_sq m_rect m_rect m_p m_p m_p m_p]*AL.c;   %heat capacities vector
eps = [eps_coat eps_coat eps_coat eps_coat eps_coat ...             %emissivity vector
    eps_sp eps_sp eps_sp eps_sp eps_sp];

R = [C(1) res.R_12 res.R_13 0 res.R_15 res.R_16 0 0 0 0;...
     res.R_12 C(2) res.R_23 res.R_24 0 res.R_26 res.R_27 0 0 0;...
     res.R_13 res.R_23 C(3) res.R_34 res.R_45 0 0 0 0 0;...
     0 res.R_24 res.R_34 C(4) res.R_45 res.R_46 0 0 0 0;...
     res.R_15 0 res.R_35 res.R_45 C(5) res.R_56 0 0 res.R_59 0;...
     res.R_16 res.R_26 0 res.R_46 res.R_56 C(6) 0 0 0 0;...
     0 res.R_27 0 0 0 0 C(7) res.R_78 0 0;...
     0 0 0 0 0 0 res.R_78 C(8) 0 0;...
     0 0 0 0 res.R_59 0 0 0 C(9) res.R_910;...
     0 0 0 0 0 0 0 0 res.R_910 C(10)];
 
%computation of the matrix H of the state-space representation of the
%linearized system
H = zeros(10); %preallocation
for i=1:10
        for j=1:10
            if i == j
                index=R(i,:)==0;
                H(i,j) = ( -4*eps(i)*A(i)*sigma*T0(i)^3 - sum(1./R(i,index==0))-1/R(i,j) )/R(i,j);

            elseif R(i,j) == 0
                H(i,j) = 0;

            else
                H(i,j) = 1/(R(i,j)*R(i,i));
            end
        end
end

%computation and plotting of the eigenvalues of the system
plot(real(eig(H)),imag(eig(H)),'x');
grid on;
xlabel('Real(\lambda)'); ylabel('Imag(\lambda)')
title('Thermal system eigenvalues for t = 0')
end

