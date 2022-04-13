%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   LUMIO Solar Array Drive Assembly                      %
%                                                                         %
%                   Balossi Claudia       952954                          %
%                   Corradetti Mario      952839                          %
%                   Al Naber Alex Yousef  939915                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clearvars; close all;

addpath('time'); 
addpath('Functions');
 
data.code = -1; 
data.time_OD = -1;
while data.code ~= 0 && data.code ~= 1
data.code = input('Choose 0 for nominal conditions or 1 for off design conditions (one solar array gets stuck): ');
if data.code == 1
    while data.time_OD < 0 || data.time_OD > 40 
    data.time_OD = input('\nChoose time at which the solar array stops moving, 0 days < t < 40 days: ');
    end
end
end
data.time_OD = data.time_OD*24; % in hours


%% ORBIT DEFINITION  
fprintf('\nOrbit definition ...')
% LUMIO ORBIT -------------------------------------------------------------
% orbit.mat is a 12-cols matrix containing date (6), position vector (3)
% and velocity vector (3) in ECEF
% [year, month, day, hour, min, sec, rx, ry, rz, vx, vy, vz]

load('orbit.mat')
date = orbit(:,1:6);            % [year, month, day, hour, min, sec]
R_Lumio_Earth = orbit(:,7:9);   % [rx, ry, rz]
% -------------------------------------------------------------------------

% EARTH AND MOON ORBITS ---------------------------------------------------
eps = deg2rad(23.4);            % obliquity of the ecliptic [rad]
R_M = astroConstants(30);       % Moon radius [km]
R_E = astroConstants(23);       % Earth radius [km] 

% vectors initialization
mjd2000 = zeros(size(date,1),1);
kep_Earth = zeros(size(date,1),6);  
R_Earth_Sun = zeros(length(mjd2000),3);
R_Moon_Earth = zeros(length(mjd2000),3);

for i = 1:size(date,1)
    % functions uplanet and ephmoon need date in modified Julian day 2000
    mjd2000(i,:) = date2mjd2000(date(i,:)); 
    % ---------------------------------------------------------------------
    % Earth orbital elements in Sun-centred ecliptic reference frame
    [kep, ksun] = uplanet(mjd2000(i,:), 3); % keplerian parameters
    [r_E, v_E] = kep2car(kep, ksun);        % cartesian coordinates 
    
    % rotate from ecliptic to equatorial plane to obtain Earth position in 
    % Sun-centred equatorial reference frame
    R_Earth_Sun(i,:) = (rot_x(-eps)*r_E)';  % [x_E, y_E, z_E]_sun
    % ---------------------------------------------------------------------
    % Moon position and velocity vector in ECEF
    [r_M, v_M] = ephMoon(mjd2000(i,:));
    R_Moon_Earth(i,:) = r_M;     % [x_M, y_M, z_M]_earth
    % ---------------------------------------------------------------------
end

% compute Lumio position wrt Moon 
R_Lumio_Moon = R_Lumio_Earth - R_Moon_Earth;
% compute Lumio position wrt Sun 
R_Lumio_Sun = R_Lumio_Earth + R_Earth_Sun;
% compute Moon position wrt Sun 
R_Moon_Sun = R_Moon_Earth + R_Earth_Sun;
% -------------------------------------------------------------------------

% Moon and Lumio orbits as seen from Earth 
figure(1)
subplot(1,2,1)
plot3(R_Moon_Earth(:,1), R_Moon_Earth(:,2), R_Moon_Earth(:,3),'LineWidth',1.2); 
hold on; grid on;
plot3(R_Lumio_Earth(:,1),R_Lumio_Earth(:,2), R_Lumio_Earth(:,3),'LineWidth',1.2);
Plot_Earth(0,0,0)
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
legend('Moon orbit', 'Lumio orbit','Location','SouthWest')
title('Moon and LUMIO orbit around Earth')
% Lumio orbit wrt Moon
subplot(1,2,2)
plot3(R_Lumio_Moon(:,1), R_Lumio_Moon(:,2), R_Lumio_Moon(:,3),'LineWidth',1.2); hold on
scatter3(R_Lumio_Moon(1,1),R_Lumio_Moon(1,2), R_Lumio_Moon(1,3),'Fill'); 
scatter3(R_Lumio_Moon(end,1),R_Lumio_Moon(end,2), R_Lumio_Moon(end,3), 'Fill'); 
Plot_Moon(0,0,0); grid on;
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
legend('Lumio orbit', 'Starting point','Final point','Location','SouthWest')
title('LUMIO orbit around Moon')


%% ATTITUDE & POINTING
fprintf('\nDefining LUMIO attitude ...')

% definition of Lumio's six faces
faces = [1 0 0;... % top (X_Lumio) (faccia 6 nel disegno)
         0 1 0;... % right with solar panel (Y_Lumio) (faccia 4 nel disegno)
         0 0 1;... % back (Z_Lumio) (faccia 2 nel disegno)
        -1 0 0;... % bottom (-X_Lumio) (faccia 5 nel disegno)
        0 -1 0;... % left with solar panel (-Y_Lumio) (faccia 3 nel disegno)
        0 0 -1];   % front (-Z_Lumio) (faccia 1 nel disegno)

days = mjd2000 - mjd2000(1);    % converts time vector from mjd2000 to days 

% vectors initialization
theta_panels = zeros(length(mjd2000),2);
eclipse = zeros(length(mjd2000),1);
alpha_faces = zeros(length(mjd2000),6);


% figure(2)
% In figure 2 the blue line represents X_Lumio and the red line points to
% % the Sun
% video1=VideoWriter('LUMIOorbit');
% video1.FrameRate=10;
% open(video1)
for i = 1: length(mjd2000)
    % ---------------------------------------------------------------------
    % definition of Lumio's body frame 
    moon_vect = R_Lumio_Moon(i,:)/norm(R_Lumio_Moon(i,:)); % unitary vector from Moon to Lumio
    sun_vect = R_Lumio_Sun(i,:)/norm(R_Lumio_Sun(i,:)); % unitary vector from Sun to Lumio
    earth_vect = R_Lumio_Earth(i,:)/norm(R_Lumio_Earth(i,:)); % unitary vector from Earth to Lumio
    % X_Lumio: unitary vector from Lumio to Moon
    X_Lumio = - moon_vect;   
    % Y_Lumio: unitary vector perpendicular to plane of Lumio_Moon and Lumio_Sun
    Y_Lumio = cross(X_Lumio, sun_vect); 
    % Z_Lumio: third component of Lumio body reference frame
    Z_Lumio = cross(X_Lumio, Y_Lumio)/norm(cross(X_Lumio, Y_Lumio));
    % ---------------------------------------------------------------------
    
    % transformation matrix from Lumio body frame to Sun centred frame
    A_Sun_Lumio = [X_Lumio', Y_Lumio', Z_Lumio']; 
    
    % (Sun coordinates in Lumio body frame) = (A_Sun_Lumio)\(Lumio coordinates in Sun body frame)
    A_Lumio_Sun = A_Sun_Lumio';  % R_L = A_LS * R_S ---> R_S = A_SL * R_L
    R_Sun_Lumio = A_Lumio_Sun * sun_vect'; % R_Sun_Lumio is directed from Sun to Lumio
    
    % ---------------------------------------------------------------------
    theta = atan2(-R_Sun_Lumio(1), -R_Sun_Lumio(3));  
    % angle between -R_Sun_Lumio and X_Lumio (R_Sun_Lumio is directed from Sun to Lumio)
    theta_panels(i,1) = theta; 
    
    theta = atan2(-R_Sun_Lumio(1), -R_Sun_Lumio(3)); % panel 2 has opposite Z axis
    theta_panels(i,2) = theta;
        
    % ---------------------------------------------------------------------
    % angle alpha between Sun direction and the normal to the surfaces of
    % Lumio (function vect-angle is on the bot of the script)
    alpha_faces(i,:) = vect_angle(faces, -R_Sun_Lumio); % minus sign because R_Sun_Lumio is directed from Sun to Lumio

    % ---------------------------------------------------------------------
    % ECLIPSE: check on the angle between Lumio-Sun and Lumio-Moon 
    rm = norm(R_Lumio_Moon(i,:));  % LUMIO - center Moon distance
    re = norm(R_Lumio_Earth(i,:)); % LUMIO - center Earth distance
    
%     % CONIC MODEL for eclipse --------------------------------------------
%     angle_moon_sun(i) = acos(dot(-moon_vect, -sun_vect)); % angle between LUMIO-Sun and LUMIO-Moon
%     angle_moon_surface(i) = acos(sqrt(rm^2 - R_M^2)/rm); % angle between LUMIO-center_of_Moon and LUMIO-tanget_Moon_surface
% 
%     if angle_moon_sun(i) <= angle_moon_surface(i)
%         eclipse(i) = 0;     % eclipse
%     else
%         eclipse(i) = 1;     % no eclipse
%     end
%    % ---------------------------------------------------------------------
%     
    % CYLINDRICAL MODEL: Eclipse wrt Moon
    cosphi = dot(moon_vect,-R_Moon_Sun(i,:))/(norm(-R_Moon_Sun(i,:)));
    if cosphi<0 && rm*sqrt(1-cosphi^2)<R_M
        eclipse(i) = 0;     % eclipse
    else
        eclipse(i) = 1;     % no eclipse
    end
    % CYLINDRICAL MODEL: Eclipse wrt Earth
    cosphi = dot(earth_vect,-R_Earth_Sun(i,:))/(norm(-R_Earth_Sun(i,:)));
    if cosphi<0 && re*sqrt(1-cosphi^2)<R_E
        eclipse(i) = 0;     % eclipse
    else
        eclipse(i) = 1;     % no eclipse
    end
      
% % %     ANIMATED PLOT of figure 2
%     x = R_Lumio_Moon(i,1);
%     y = R_Lumio_Moon(i,2);
%     z = R_Lumio_Moon(i,3);
%     
%     clf
%     p = plot3(R_Lumio_Moon(:,1), R_Lumio_Moon(:,2), R_Lumio_Moon(:,3),'LineWidth',0.5);
%     hold on
%     p = plot3(0,0,0, 'ok', 'MarkerSize',4, 'MarkerFaceColor','black');
%     p.HandleVisibility = 'Off';
%     hold on; grid on;
%     p = plot3(x,y,z, 'ob', 'MarkerSize',4, 'MarkerFaceColor','k');
%     p.HandleVisibility = 'On';
%     p = line([x X_Lumio(1)*1e4], [y, X_Lumio(2)*1e4], [z X_Lumio(3)*1e4]);
%     p.HandleVisibility = 'On';
%     l = line([0 -R_Moon_Sun(i,1)],[0 -R_Moon_Sun(i,2)],[0 -R_Moon_Sun(i,3)]);
%     l.HandleVisibility='On';
%     l.Color='r';
%     axis([-7e4 7e4 -7e4 7e4 -6e4 6e4]);
%     title('Animated LUMIO orbit around the Moon')
%     xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
%     frame1=getframe(gcf);
%     writeVideo(video1,frame1);
%     pause(0.00000001);

end
% close(video1)

figure(3)
plot(days, eclipse, 'LineWidth', 1.5); grid on
xlabel('Days'); ylabel('0: eclipse, 1: no eclipse')
title('Eclipse event')

% since the incident angle alpha is defined between [0 90] (positive and 
% negative alpha wrt face's normal are considered equivalent), angles > 90
% mean that the face is in shadow 

figure(4)
plot(days, rad2deg(alpha_faces(:,1)), days, rad2deg(alpha_faces(:,2)), days, rad2deg(alpha_faces(:,3)),...
    days, rad2deg(alpha_faces(:,4)), days, rad2deg(alpha_faces(:,5)), days, rad2deg(alpha_faces(:,6)), 'LineWidth',1.2);
hold on; grid on;
xlabel('Time [days]'); ylabel('Incident angle \alpha [deg]');
yl=ylim; xl=xlim;
Xbox=[xl(1) xl(1) xl(2) xl(2)];
Ybox=[90 yl(2) yl(2) 90];
patch(Xbox,Ybox,'black','FaceColor','black','FaceAlpha',0.1);
legend('top (X)','right (Y)','back (Z)','bottom (-X)','left (-Y)','front (-Z)','shadow zone');
title('Angle between the normal of LUMIO faces and the Sun')



%% STEPPER MOTOR
fprintf('\nRunning stepper motor analysis ...')

% load mechanical system's data
mech_data

time_hours = days.*24;  % convert time vector in hours 

data.index = 0; % rotation index of HSM for voltage switch
% Initial conditions for HSM
data.theta0 = 0;
data.Ia0 = 0;
data.Ib0 = 0;
data.w0 = 0;
% Integration time of 1 hour
data.ti = 0;
data.tf = 1*3600;

data.thpanel0=(18*pi/180)*round((theta_panels(1,1))/(18*pi/180)); % assume as initial inclination of the panel the closest multiplier of 18
tic
values = HSM(time_hours,theta_panels,data);
t = toc;
fprintf('\n integration time = %f s',t)
 
% Plot of current trend
figure(5)
sgtitle('Current trend')
    subplot(2,1,1)
        plot(values.T,values.Ia,'LineWidth',1);
        xlabel('Time [s]')
        ylabel('I_a [A]')
        set(gca,'xtick',[0:10:120])
        set(gca,'ytick',[-1:0.4:1])
        grid on
    subplot(2,1,2)
        plot(values.T,values.Ib,'r','LineWidth',1);
        xlabel('Time [s]')
        ylabel('I_b [A]')
        set(gca,'xtick',[0:10:120])
        set(gca,'ytick',[-1:0.4:1])
        grid on

%Plot of Voltage trend
figure(6)
sgtitle('Voltage trend')
    subplot(2,1,1)
        plot(values.T,values.Va,'LineWidth',1.2);
        xlabel('Time [s]')
        ylabel('V_a [V]')
        grid on
        ylim([-15 15]);
        set(gca,'xtick',[0:10:120])
        set(gca,'ytick',[-16:4:16])
    subplot(2,1,2)
        plot(values.T,values.Vb,'r','LineWidth',1.2);
        xlabel('Time [s]')
        ylabel('V_b [V]')
        grid on
        ylim([-15 15]);
        set(gca,'xtick',[0:10:120])
        set(gca,'ytick',[-12:6:12])


figure(7)
plot(time_hours/(24), rad2deg(theta_panels(:,1)),'LineWidth',0.8)
hold on; grid on
plot(time_hours/(24), rad2deg(values.thetaHSM(:,1)), time_hours/(24), rad2deg(values.thetaHSM(:,2)), 'LineWidth',1)
xlabel('Time [days]'); ylabel('Angle \theta [deg]');
legend('Nominal angle', 'Right solar panel', 'Left solar panel', 'Location', 'southeast')
title('Real and SA angle')

figure(8)
plot(rad2deg(values.w),values.Torque)
grid on
title('Torque vs angular velocity')
set(gca,'xtick',[-280:50:280])
set(gca,'ytick',[-3:0.4:3])
xlabel('Angular velocity [deg/s]');
ylabel('Torque [Nm]')

% EIGENVALUES OF MECHANICAL SYSTEM ----------------------------------------
figure(9)
eig_HSM(data)


%% Animated plot of solar arrays movement
pause(3) 
figure('units','normalized','outerposition',[0 0 1 1])
video2=VideoWriter('LUMIO_SA');
video2.FrameRate=10;
open(video2)
for i=1:5:length(values.thetaHSM)
    if i<length(values.thetaHSM)+1 % control to clear the figure
        clf
    end
    subplot(1,2,1)
        plot(time_hours/24, rad2deg(theta_panels(:,1)))
        hold on; grid on
        plot(time_hours(1:i)/24, rad2deg(values.thetaHSM(1:i,1)),time_hours(1:i)/24, rad2deg(values.thetaHSM(1:i,2)));
        xlabel('Time [days]'); ylabel('Angle \theta [deg]');
        legend('Nominal angle', 'Solar array 1','Solar array 2', 'Location', 'Best')

    subplot(1,2,2)
        sat(rad2deg(values.thetaHSM(i,1)),rad2deg(values.thetaHSM(i,2)));
        frame2=getframe(gcf);
        writeVideo(video2,frame2);
    pause(0.000001);
end
close(video2)
pause(3)


%% THERMAL MODEL
% takes as input the incidence angle alpha of LUMIO faces and panels
% LUMIO faces: alpha_faces [rad]
% Panels: dtheta_panels [rad]

fprintf('\nRunning thermal analysis ...')

% Definition of heater power ----------------------------------------------
data.pow.heat_off = 0;  % no heater
data.pow.heat_on = 15;  % heater power [W]
%--------------------------------------------------------------------------

% Definition of incidence angle on solar panels
% theta_panels:     vector of nominal angle to be reached by solar panels
% values.thetaHSM : vector of "stepped" angle actually reached by the panels
dtheta_panels = abs(values.thetaHSM - theta_panels);  % incidence angle of sun rays on front side of the panel
dtheta_panels_back = pi-dtheta_panels; % incidence angle of sun rays on back side of the panel
alpha_panels = [dtheta_panels, dtheta_panels_back];

% computing cosines
cos_f = cos(alpha_faces);
cos_p = cos(alpha_panels);
cos_mat = [cos_f, cos_p];

% Interpolation
time_hours = days.*24;
interp_cos = interp1(time_hours, cos_mat,'spline','pp');

T0 = [293, 293, 293, 293, 293, 293, 313, 313, 313, 313];    % initial guess
options = odeset('RelTol',1e-13,'OutputFcn',@odewbar);
tic;
[tt,T] = ode15s(@therm_fun, [0 time_hours(end)*3600], T0, options, interp_cos, data);
t=toc;
fprintf('\n integration time = %f s',t)



% Plot temperatures profile
figure('units','normalized','outerposition',[0 0 1 1])    
subplot(2,1,1)
    plot(tt/86400, T(:,1)-273.15,'LineWidth',1.3)
    hold on; grid on
    plot(tt/86400,T(:,2)-273.15,'LineWidth',1.3)
    plot(tt/86400,T(:,3)-273.15,'LineWidth',1.3)
    plot(tt/(86400),T(:,4)-273.15,'LineWidth',1.3)
    plot(tt/(86400),T(:,5)-273.15,'LineWidth',1.3)
    plot(tt/(86400),T(:,6)-273.15,'LineWidth',1.3)
    legend('Face 1', 'Face 2', 'Face 3', 'Face 4', 'Face 5', 'Face 6','Location', 'SouthWest')
    xlabel('Time [days]'); ylabel('T [C]')
    title('Temperature of LUMIO faces')
subplot(2,1,2)
    plot(tt/(86400),T(:,7)-273.15,'-','LineWidth',1.3)
    hold on; grid on
    plot(tt/86400,T(:,8)-273.15,'-','LineWidth',1.3)
    plot(tt/(86400),T(:,9)-273.15,'-','LineWidth',1.3)
    plot(tt/(86400),T(:,10)-273.15,'-','LineWidth',1.3)
    xlabel('Time [days]'); ylabel('T [C]')
    title('Temperature of solar panels')
    legend('SP dx: 7', 'SP dx: 8', 'SP sx: 9', 'SP sx: 10','Location', 'SouthWest')

    
% EIGENVALUES OF THERMAL SYSTEM ----------------------------------------    
figure(11)
eig_HEAT(T0)


%% SYSTEM PARAMETRIC RESPONSE
figure(12)
sgtitle('Current intensity over one single step with modified R and L')
for i=0:1
data.index = 0; % rotation index of HSM for voltage switch
% Initial conditions for HSM
data.theta0 = 0;
data.Ia0 = 0;
data.Ib0 = 0;
data.w0 = 0;
% Integration time of 1 hour
data.ti = 0;
data.tf = 1*3600;
data.thpanel0=(18*pi/180)*round((theta_panels(1,1))/(18*pi/180)); % assume as initial inclination of the panel the closest multiplier of 18
var='R';
tic
values = HSM(time_hours,theta_panels,data,i,var);
subplot(1,2,1)
    plot(values.T, values.Ia,'LineWidth',1)
    hold on; grid on
    axis([0 0.5 0 0.75])
subplot(1,2,2)
    plot(values.T, abs(values.Ib),'LineWidth',1)
    hold on; grid on
    axis([0 0.5 0 1.4])
end
for i=1:1
data.index = 0; % rotation index of HSM for voltage switch
% Initial conditions for HSM
data.theta0 = 0;
data.Ia0 = 0;
data.Ib0 = 0;
data.w0 = 0;
% Integration time of 1 hour
data.ti = 0;
data.tf = 1*3600;
data.thpanel0=(18*pi/180)*round((theta_panels(1,1))/(18*pi/180)); % assume as initial inclination of the panel the closest multiplier of 18
var='L';
tic
values = HSM(time_hours,theta_panels,data,i,var);
subplot(1,2,1)
    plot(values.T, values.Ia, '-.','LineWidth',1)
    hold on; grid on
    axis([0 0.5 0 0.75])
subplot(1,2,2)
    plot(values.T, abs(values.Ib), '-.','LineWidth',1)
    hold on; grid on
    axis([0 0.5 0 1.4])
end

subplot(1,2,1)
xlabel('Time [s]'); ylabel('I_a [A]')
legend('Nominal conditions','R = 10.6 \Omega','L = 3.2e-3 H','Location', 'Best')
subplot(1,2,2)
xlabel('Time [s]'); ylabel('current I_b [A]')
legend('Nominal conditions','R = 10.6 \Omega','L = 3.2e-3 H', 'Location', 'SouthEast')


%% PRESENTATION
c=['b','r','y','m'];
for i=0:3
data.index = 0; % rotation index of HSM for voltage switch
% Initial conditions for HSM
data.theta0 = 0;
data.Ia0 = 0;
data.Ib0 = 0;
data.w0 = 0;
% Integration time of 1 hour
data.ti = 0;
data.tf = 1*3600;
data.thpanel0=(18*pi/180)*round((theta_panels(1,1))/(18*pi/180)); % assume as initial inclination of the panel the closest multiplier of 18
var='R';
values = HSM(time_hours,theta_panels,data,i,var);
figure(1)
subplot(2,2,1)
plot(values.T, abs(rad2deg(values.w)),c(i+1),'LineWidth',1)
hold on; grid on
xlim([0 0.12])
subplot(2,2,1)
xlabel('time [s]')
ylabel(' \omega [rad/s]')
title('Angular velocity decreasing R')

subplot(2,2,3)
plot(rad2deg(values.w),values.Torque,c(i+1),'LineWidth',1)
hold on; grid on
xlim([0 360])
xlabel('\omega [rad/s]')
ylabel(' Torque [Nm]')
title('Torque decreasing R')
end

for i=0:3
data.index = 0; % rotation index of HSM for voltage switch
% Initial conditions for HSM
data.theta0 = 0;
data.Ia0 = 0;
data.Ib0 = 0;
data.w0 = 0;
% Integration time of 1 hour
data.ti = 0;
data.tf = 1*3600;
data.thpanel0=(18*pi/180)*round((theta_panels(1,1))/(18*pi/180)); % assume as initial inclination of the panel the closest multiplier of 18
var='L';
values = HSM(time_hours,theta_panels,data,i,var);
figure(1)
subplot(2,2,2)
plot(values.T, abs(rad2deg(values.w)),c(i+1),'LineWidth',1)
hold on; grid on
xlim([0 0.006])
xlabel('time [s]')
ylabel(' \omega [rad/s]')
title('Angular velocity increasing L')

subplot(2,2,4)
plot(rad2deg(values.w),values.Torque,c(i+1),'LineWidth',1)
hold on; grid on
xlim([0 360])
xlabel('\omega [rad/s]')
ylabel(' Torque [Nm]')
title('Torque increasing L')
end




%% USEFUL FUNCTIONS
function alpha_faces = vect_angle(faces, vect)
% vect_angle calculates the angle btw faces of Lumio and Sun direction
alpha_faces = zeros(1,6);
for i = 1: 6
    face = faces(i,:); 
    alpha_faces(i) = acos(dot(face, vect)/(norm(face)*norm(vect)));
end
end

function [vect2]=rotatevect(vect1,theta,k)
%%RODRIGUES'S FORMULA
vect2=vect1*cosd(theta)+cross(k,vect1)*sind(theta)+k*dot(k,vect1)*(1-cosd(theta));
end