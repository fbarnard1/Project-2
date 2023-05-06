%% Project 2 - Question 2
clc; clear; close all;

% given data values
gamma = 1.4; R = 287; r = 8.4; L = 12; S = 8; V0 = 50; b = 9; T1 = 300;
theta_s = 3*pi/2;  q_in = 2.8*10^6; Tw = 300;E = 8/(2*12); theta_b = pi;

% given span for the crank angle
tspan = [pi 3*pi];

% initial values that are given for the pressure and the work
y0 = [1.013*10^5 0];

% initial value given for the temperature
y2 = 300;

% ODE function to solve for the work and pressure when w=50
sol = ode45(@Q1,tspan,y0);

% ODE function to solve for the temperature when w=50
sol2 = ode45(@T2,tspan,y2);

% ODE function to solve for the work and pressure when w=100
sol3 = ode45(@Q2,tspan,y0);

% ODE function to solve for the temperature when w=100
sol4 = ode45(@T3,tspan,y2);

% creating the pressure vs. crank angle plot, and adding both lines when
% w = 50 and when w = 100
subplot(2,2,1)
plot(smooth(sol.y(1,:)),'r')
hold on
plot(smooth(sol3.y(1,:)),'b')
grid on
xlabel('Crank Angle (rad)')
ylabel('Pressure (Pa)')
legend('50 rad/s','100 rad/s')
title('Pressure vs. Crank Angle')

% creating the work vs. crank angle plot, and adding both lines when
% w = 50 and when w = 100
subplot(2,2,2)
plot(smooth(sol.y(2,:)),'r')
hold on
plot(smooth(sol3.y(2,:)),'b')
grid on
xlabel('Crank Angle (rad)')
ylabel('Work (J)')
legend('50 rad/s','100 rad/s')
title('Work vs. Crank Angle')

% creating the volume vs. crank angle plot, and adding both lines when
% w = 50 and when w = 100
subplot(2,2,3)
fplot(@V2,tspan,'r')
hold on
fplot(@V3,tspan,'b')
grid on
xlabel('Crank Angle (rad)')
ylabel('Volume (cm^3)')
legend('50 rad/s','100 rad/s')
title('Volume vs. Crank Angle')

% creating the temperature vs. crank angle plot, and adding both lines when
% w = 50 and when w = 100
subplot(2,2,4)
plot(sol2.y(1,:),'r')
hold on
plot(sol4.y(1,:),'b')
grid on
xlabel('Crank Angle (rad)')
ylabel('Temperature (K)')
legend('50 rad/s','100 rad/s')
title('Temperature vs. Crank Angle')


%% Function to solve for the pressure and work when w = 50
function [dP] = Q1(theta,p)
    gamma = 1.4;   b = 9; qin = 2.8*10^6; Tw = 300; w = 50; h = 50; Aw = (4*V2(theta))/b;
    dP = [((-gamma*(p(1)*dV2(theta))/V2(theta)))+((gamma-1)*((M2(theta)*qin)/V2(theta)*dx2(theta)))-((gamma-1)*h*(Aw)*(T2(theta,p)-Tw)/(V2(theta)*w))-((gamma*M2(theta))/(M2(theta)*w))*p(1);
          p(1)*dV2(theta)];
end

%% Function to solve for the pressure and work when w = 50
function [dP2] = Q2(theta,p)
    gamma = 1.4;  b = 9; qin = 2.8*10^6; Tw = 300; w2 = 100; h = 50; Aw = (4*V2(theta))/b;
    dP2 = [((-gamma*(p(1)*dV2(theta))/V2(theta)))+((gamma-1)*((M3(theta)*qin)/V2(theta)*dx2(theta)))-((gamma-1)*h*(Aw)*(T3(theta,p)-Tw)/(V2(theta)*w2))-((gamma*M3(theta))/(M3(theta)*w2))*p(1);
           p(1)*dV2(theta)];
end

%% Function for the change in x with respect to theta
function Derx = dx2(theta)
theta_s = 3*pi/2; theta_b = pi;
if (pi<= theta)&& (theta <= theta_s)
    Derx = 0;
elseif (theta_s <= theta) && (theta <= theta_s+theta_b)
    Derx = sin(theta - (3*pi)/2)/2;
elseif  (theta_s+theta_b <=  theta) && (theta<=3*pi)
    Derx = 0;
end
end

%% Function to find the change in volume
function deltaV = dV2(theta)
    deltaV = 185*sin(theta) + (185*cos(theta)*sin(theta))/(3*(1 - sin(theta)^2/9)^(1/2));
end

%% function to solve for the temperature when w = 50
function Temp = T2(theta,p)
    R = 287;  
    Temp = (p(1)*V2(theta))/((M2(theta))*(R));

end

%% function to solve for the temperature when w = 100
function Temp = T3(theta,p)
    R = 287; 
    Temp = (p(1)*V2(theta))/((M3(theta))*(R));
end

%% function to solve for the change in mass when w = 50

function mass = M2(theta)
    M0 = 58.8269; C = .8; w = 50;
    mass = M0*exp((-C/w)*(theta-pi));
end

%% function to solve for the change in mass when w = 100
function mass = M3(theta)
    M0 = 58.8269; C = .8; w2 = 100;
    mass = M0*exp((-C/w2)*(theta-pi));
end

%% function to solve for the change in volume when w = 50
function Volume = V2(theta)
    r = 8.4; V0 = 50; E = 8/(2*12);
    Volume = V0*(1+((r-1)/(2*E))*(1+E*(1-cos(theta))-sqrt(1-E^2*(sin(theta))^2)));
end

%% function to solve for the change in volumme when w = 100
function Volume = V3(theta)
    r = 8.4; V0 = 50; E = 8/(2*12);
    Volume = V0*(1+((r-1)/(2*E))*(1+E*(1-cos(theta))-sqrt(1-E^2*(sin(theta))^2)));
end
