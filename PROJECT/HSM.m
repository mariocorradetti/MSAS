function values = HSM(time_hours,theta_panels, data,i,var,varargin)
% Function to analyze the behaviour of the HSM during the assigned time
% INPUTS
% time_hours:         orbit time in hours
% theta_panels:       nominal angle to be reached by panels 
% data:
% 
% OUTPUTS
% values:     structure containing   - thetaHSM (values of the angle of the solar array)
%                                    - Ia & Ib
%                                    - Va & Vb
%                                    - w & th (angular velocity and angular displacement)
%                                    - T (time vector)

% Load parameters of motor and satellite
mech_data;
if nargin == 5
    if var == 'R'
    motor.R = motor.R-i*3;
    elseif var == 'L'
    motor.L = motor.L+i*1.2e-3;
    end
end
k=0; % Integer factor for rotation in event function
y0=[data.Ia0 data.Ib0 data.w0 data.theta0];
% Initialize field thpanel0
values.thpanel0(1,1)=data.thpanel0;
% Initialize output vectors
values.thetaHSM(1,1) = data.thpanel0;
values.thetaHSM(1,2) = data.thpanel0;
values.Ia=[];
values.Ib=[];
values.w=[];
values.th=[];
values.T=[];
values.Va=[];
values.Vb=[];
values.Torque=[];

row = 1; % Initialize the row index for values struct
j = 1;

wb = waitbar(0,'Integration of mechanical system ODEs','Name','HSM');

for i = 1:length(time_hours)-1  % end-1 because I need to compare 2 consecutive elements
    
    % Step forward , rotate if the angle between current position and the final
    % position is greater than half of step angle of the solar array (step angle of sp = 18/tau = 4.5)
    if theta_panels(i+1,1)-theta_panels(i,1) > 0 && theta_panels(i,1) - values.thpanel0(row,1) > 0.5*(18*pi/180)/motor.tau
        data.index = data.index + 1; % Update rotation index
        
        % Bound rotation index between 1 and 4 ---------------------------
        if data.index == 5
            data.index = 1;
        end
        if data.index == 0
            data.index = 4;
        end
        % -----------------------------------------------------------------
        k = k+1;
        
        % Integration HSM dynamics
        dt = [data.ti,data.tf]; % integration interval
        options = odeset('RelTol',1e-12,'Events',@(t,y) eventfun(t,y,k,motor.theta_step));
        [time,y] = ode23tb(@(t,y)HSMeq(t,y,data,motor,SP,yoke), dt, y0, options);
        
        for p = 1 : length(time)
            [~,par] = HSMeq(time(p),y(p,:),data,motor,SP,yoke);
            values.Va = [values.Va; par(1)];
            values.Vb = [values.Vb; par(2)];
            values.Torque = [values.Torque; par(3)];
        end
        
        % Update initial conditions
        y0 = y(end,:);
        % Update time interval
        data.ti = time(end);
        data.tf = time(end)*2;
        
        %output values
        values.Ia = [values.Ia; [y(:,1);NaN]];
        values.Ib = [values.Ib; [y(:,2);NaN]];
        values.w = [values.w; [y(:,3);NaN]];
        values.th = [values.th; [y(:,4);NaN]];
        values.T = [values.T; [time;NaN]];
        values.Va = [values.Va; NaN];
        values.Vb = [values.Vb; NaN];
        values.Torque = [values.Torque; NaN];
        row = row+1;
        j = j+1;
        
        values.thpanel0(row,1) = values.thpanel0(row-1,1)+(18*pi/180)/motor.tau;
        values.thetaHSM(j,1) = values.thetaHSM(j-1,1) + (18*pi/180)/motor.tau;
        
        % Control for thetaHSM for design and off-design conditions
        if data.code == 0
            values.thetaHSM(j,2) = values.thetaHSM(j-1,2) + (18*pi/180)/motor.tau;
        elseif data.code == 1
            if time_hours(j) <= data.time_OD
                values.thetaHSM(j,2) = values.thetaHSM(j-1,2) + (18*pi/180)/motor.tau;
            else
                values.thetaHSM(j,2) = values.thetaHSM(j-1,2);
            end
        end
        
    % Step backward, rotate if the angle between current position and the final
    % position is less than half of step angle of the solar array (step angle of sp = 18/tau = 4.5)
    elseif theta_panels(i+1,1)-theta_panels(i,1) < 0 && theta_panels(i,1) - values.thpanel0(row,1)< -0.5*(18*pi/180)/motor.tau
        
        data.index = data.index-1; % update rotation index
        
        % Bound rotation index between 1 and 4 ----------------------------
        if data.index == -1 % in case i start with index=0
            data.index = 3;
        end
        if data.index == 0
            data.index = 4;
        end
        % -----------------------------------------------------------------
        
        k = k-1;
        dt = [data.ti,data.tf]; %integration interval
        options = odeset('RelTol',1e-12,'Events',@(t,y) eventfun(t,y,k,motor.theta_step));
        [time,y] = ode23tb(@(t,y) HSMeq(t,y,data,motor,SP,yoke), dt, y0, options);
        
        for p = 1 :length(time)
            [~,par] = HSMeq(time(p),y(p,:),data,motor,SP,yoke);
            values.Va = [values.Va; par(1)];
            values.Vb = [values.Vb; par(2)];
            values.Torque = [values.Torque; par(3)];
        end
        
        % Update initial conditions
        y0 = y(end,:);
        % Update time interval
        data.ti = time(end);
        data.tf = time(end)*2;
        %output values
        values.Ia = [values.Ia; [y(:,1);NaN]];
        values.Ib = [values.Ib; [y(:,2);NaN]];
        values.w = [values.w; [y(:,3);NaN]];
        values.th = [values.th; [y(:,4);NaN]];
        values.T = [values.T; [time;NaN]];
        values.Va = [values.Va; NaN];
        values.Vb = [values.Vb; NaN];
        values.Torque = [values.Torque; NaN];

        row = row+1;
        j = j+1;
        values.thpanel0(row,1) = values.thpanel0(row-1,1)-(18*pi/180)/motor.tau;
        values.thetaHSM(j,1) = values.thetaHSM(j-1,1) - (18*pi/180)/motor.tau;
        
        if data.code == 0
            values.thetaHSM(j,2) = values.thetaHSM(j-1,2) - (18*pi/180)/motor.tau;
        elseif data.code == 1
            if time_hours(j) <= data.time_OD
                values.thetaHSM(j,2) = values.thetaHSM(j-1,2) - (18*pi/180)/motor.tau;
            else
                values.thetaHSM(j,2) = values.thetaHSM(j-1,2);
            end
        end
        
        
    else % case where my motor is not turn on
        j = j+1;
        values.thetaHSM(j,1) = values.thetaHSM(j-1,1);
        values.thetaHSM(j,2) = values.thetaHSM(j-1,2);
        
    end
    
    waitbar(i/(length(time_hours)-1))
end
delete(wb)
end
