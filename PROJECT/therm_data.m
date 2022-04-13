% thermal system's data

body.h = 226e-3;        %short side body length [m]
body.L = 340e-3;        %long side body length [m]
body.t = 1.5e-3;        %body thickness [m]
body.A_tot = 409.4e-3;  %total area of the body [m^2]
body.m = 2700*body.A_tot*body.t;    %body mass [kg]

SC.alpha = 0.90;         %solar cell absorptivity
SC.eps = 0.80;           %solar cell emissivity
SC.eta = 0.15;           %solar cells efficiency
SC.A = 16*26.51e-4;      %area occupied by solar cells [m^2]

back.alpha = 0.8;   %back of solar panel material absorptivity
back.eps = 0.5;     %back of solar panel material emissivity

SP.L = 325e-3;              % solar panel length [m]
SP.h = 205e-3;              % solar panel heigth [m]
SP.t = 5e-3;                % solar panel thickness [m]
SP.m = 0.5;                 % solar panel mass [kg]
SP.J = 1/12*SP.m*SP.h^2;    % inertia of one panel [Kg m^2]

AL.k = 160;         %aluminum Al 6061 T6 thermal conductivity [W/(m*K)]
AL.alpha = 0.031;   %aluminum Al 6061 T6 absorptivity
AL.eps = 0.039;     %aluminum Al 6061 T6 emissivity
AL.rho = 2700;      %aluminum Al 6061 T6 density [kg/m^3]
AL.c = 900;         %aluminum heat capacity [J/(kg*K)]

AU.alpha = 0.3;     %gold absorptivity
AU.eps = 0.03;      %gold emissivity

AG.alpha = 0.163;   %silvered teflon absorptivity
AG.eps = 0.8;       %silvered teflon emissivity


% contact resistance between surfaces [K/W]
R_long = body.L/(AL.k*body.h*body.t);     % contact resistance for the long side
R_short = body.h/(AL.k*body.L*body.t);    % contact resistance for the short side
R_square = body.h/(AL.k*body.h*body.t);   % contact resistance for the squared panel
R_yoke = yoke.L/(AL.k*yoke.h*yoke.t);     % contact resistance for the yolks
R_panel = SP.L/(2*AL.k*SP.h*SP.t);        % contact resistance for the solar panels


m_sq = AL.rho*body.t*body.h^2;          
m_rect = AL.rho*body.t*body.h*body.L;  
m_p = AL.rho*SP.t*SP.L*SP.h;


res.R_12 = R_square + R_long;
res.R_13 = res.R_12;
res.R_15 = res.R_12;
res.R_16 = res.R_12;
res.R_24 = res.R_12;
res.R_34 = res.R_12;
res.R_45 = res.R_12;
res.R_46 = res.R_12;

res.R_23 = 2*R_short;
res.R_26 = res.R_23;
res.R_35 = res.R_23;
res.R_56 = res.R_23;

res.R_27 = R_yoke + R_panel;
res.R_78 = 2*R_panel;

res.R_59 = res.R_27;
res.R_910 = res.R_78;
