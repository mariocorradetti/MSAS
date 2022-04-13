function [r_ECI, v_ECI] = kep2car(kep,mu)

a = kep(1);
e = kep(2);
i = kep(3);
Om = kep(4);
om = kep(5);
th = kep(6);


p = a*(1-e^2);
h = sqrt(p*mu);
r = p / (1+e*cos(th));

r_PF = r*[cos(th), sin(th), 0]';
v_PF = (mu/h) * [-sin(th), (e+cos(th)), 0]';

% Rotation matrices:
R_om = [cos(om)  sin(om)    0   ;
        -sin(om) cos(om)    0   ;
           0        0       1   ];
       
R_i  = [   1        0       0   ;
           0     cos(i)   sin(i);
           0    -sin(i)   cos(i)];
       
R_Om = [cos(Om)  sin(Om)    0   ;
        -sin(Om) cos(Om)    0   ;
           0        0       1   ];
    
R313 = R_om * R_i * R_Om; % ECI -> PF


% PF -> ECI
r_ECI = R313'* r_PF;
v_ECI = R313'* v_PF;


end

