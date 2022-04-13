function [T_dot] = therm_fun (time, T, interp_cos, data)
% Function that describes the thermal dynamic of the system  
%
% INPUT 
% time        instant of time [s]
% T           vector containing the temperature of the nodes [K]
% interp_cos  coefficient of interpolation 
% data        structure containing inital conditions and other useful data
%
% OUTPUTS
% T_dot       derivatives of the temperatures


% load useful data
mech_data
therm_data

% Set minimum temperatur for LUMIO's faces (from references)
T_min = -5+273.15; 

% Heater is on when at least one face has tempertaure lower than T_min or if the solar panel is broken 
if data.code ==1 && time < data.time_OD*3600 || data.code == 0
    if sum(T < T_min)>= 1
        heat = data.pow.heat_on;
    else
        heat = data.pow.heat_off;
    end
else
    heat = data.pow.heat_on;
end

Qdiss = 0.8*4.2;                    % dissipated power
Q = Qdiss + heat;                   % total power
q = Q/(2*body.h^2 + 4*body.h*body.L);
pow.Q_sq = body.h^2*q;              % power contribution on square face
pow.Q_rect = body.L*body.h*q;       % power contribution on rectangular face 

E_s = 1360.8;                       % mean solar power [W/m^2]
sigma = 5.67037e-8;                 % Stefan-Boltzmann constant [W/(m^2*K^4)]

alpha_coat = (AL.alpha*48 + AU.alpha*27 + AG.alpha*25)/100; %thermal coating absorptivity
eps_coat = (AL.eps*48 + AU.eps*27 + AG.eps*25)/100;         %thermal coating emissivity
c = AL.c;       % heat capacity [J/(kg*K)]

A_sq = body.h^2;                    % area of square face [m^2]
A_rect = body.h*body.L;             % area of rectangular face [m^2]
A_p = SP.L*SP.h;                    % area of solar panel face [m^2]

m_sq = A_sq*body.t*AL.rho;      % mass of square face [kg]
m_rect = A_rect*body.t*AL.rho;  % mass of rectangular face [kg]
m_p = A_p*SP.t*AL.rho;          % mass of solar panel [kg]

cos_alpha = ppval(interp_cos,time/3600);  % evaluations of spline

% check on sign of cos_alpha: do not consider the contribution of radiation 
% absorbed by the sides corresponding to negative cosine
for i = 1:length(cos_alpha)
    if cos_alpha(i) < 1e-10
        cos_alpha(i) = 0;
    end
end

T1 = T(1) ;     T2 = T(2) ;
T3 = T(3) ;     T4 = T(4) ;
T5 = T(5) ;     T6 = T(6) ;
T7 = T(7) ;     T8 = T(8) ;
T9 = T(9) ;     T10 = T(10) ;

Q1 = alpha_coat*A_sq*E_s*cos_alpha(1) + pow.Q_sq - eps_coat*A_sq*sigma*T1^4 ;
T1_dot = ((T5-T1)/res.R_15 + (T3-T1)/res.R_13 + (T2-T1)/res.R_12 + (T6-T1)/res.R_16 + Q1)/(m_sq*c);

Q2 = alpha_coat*A_rect*E_s*cos_alpha(2) + pow.Q_rect - eps_coat*A_rect*sigma*T2^4;
T2_dot = ((T4-T2)/res.R_24 + (T3-T2)/res.R_23 + (T6-T2)/res.R_26 + (T1-T2)/res.R_12 + (T7-T2)/res.R_27 + Q2) / (m_rect*c);

Q3 = alpha_coat*A_rect*E_s*cos_alpha(3) + pow.Q_rect - eps_coat*A_rect*sigma*T3^4;
T3_dot = ((T1-T3)/res.R_13 + (T2-T3)/res.R_23 + (T5-T3)/res.R_35 + (T4-T3)/res.R_34 +Q3)/(m_rect*c);

Q4 = alpha_coat*A_sq*E_s*cos_alpha(4) + pow.Q_sq - eps_coat*A_sq*sigma*T4^4;
T4_dot = ( (T2-T4)/res.R_24 + (T5-T4)/res.R_45 + (T6-T4)/res.R_46 + (T3-T4)/res.R_34 + Q4)/(m_sq*c);

Q5 = alpha_coat*A_rect*E_s*cos_alpha(5) + pow.Q_rect - eps_coat*A_rect*sigma*T5^4;
T5_dot = ( (T4-T5)/res.R_45 + (T1-T5)/res.R_15 + (T3-T5)/res.R_35 + (T6-T5)/res.R_56 + (T9-T5)/res.R_59 + Q5)/(m_rect*c);

Q6 = alpha_coat*A_rect*E_s*cos_alpha(6) + pow.Q_rect - eps_coat*A_rect*sigma*T6^4;
T6_dot = ( (T5-T6)/res.R_56 + (T1-T6)/res.R_16 + (T4-T6)/res.R_46 + (T2-T6)/res.R_26 + Q6)/(m_rect*c);

Q7 =  (((SC.alpha-SC.eta)*SC.A + AL.alpha*(A_p-SC.A))*cos_alpha(7) + (cos_alpha(9))*back.alpha*A_p)*E_s...       
    - (SC.eps*SC.A + AL.eps*(A_p-SC.A) + back.eps*A_p)*sigma*T7^4;
T7_dot = ( (T2-T7)/res.R_27 + (T8-T7)/res.R_78 + Q7)/(m_p*c);

Q8 =  (((SC.alpha-SC.eta)*SC.A + AL.alpha*(A_p-SC.A))*cos_alpha(7) + (cos_alpha(9))*back.alpha*A_p)*E_s...      
    - (SC.eps*SC.A + AL.eps*(A_p-SC.A) + back.eps*A_p)*sigma*T8^4;
T8_dot = ( (T7-T8)/res.R_78 + Q8)/ (m_p*c);

% panels 9-10 are subjected to off nominal design
if data.code == 0 % nominal conditions: i need cos_alpha(7) for front and cos_alpha(9) for back
    m = 7;
    n = 9;
elseif data.code ==1  % off design
    m = 8;
    n = 10;
end
    
Q9 = (((SC.alpha-SC.eta)*SC.A + AL.alpha*(A_p-SC.A))*cos_alpha(m) + (cos_alpha(n))*back.alpha*A_p)*E_s... 
    - (SC.eps*SC.A + AL.eps*(A_p-SC.A) + back.eps*A_p)*sigma*T9^4;
T9_dot = ( (T5-T9)/res.R_59 + (T10-T9)/res.R_910 + Q9)/(m_p*c);

Q10 = (((SC.alpha-SC.eta)*SC.A + AL.alpha*(A_p-SC.A))*cos_alpha(m) + (cos_alpha(n))*back.alpha*A_p)*E_s...  
    - (SC.eps*SC.A + AL.eps*(A_p-SC.A) + back.eps*A_p)*sigma*T10^4;
T10_dot = ( (T9-T10)/res.R_910 + Q10 )/(m_p*c); 


T_dot(1,1) = T1_dot;
T_dot(2,1) = T2_dot;
T_dot(3,1) = T3_dot;
T_dot(4,1) = T4_dot;
T_dot(5,1) = T5_dot;
T_dot(6,1) = T6_dot;
T_dot(7,1) = T7_dot;
T_dot(8,1) = T8_dot;
T_dot(9,1) = T9_dot;
T_dot(10,1) = T10_dot;

end
