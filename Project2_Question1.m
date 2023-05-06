%% Project 2 - Question 1

clc;clear;close all;

% initial values that are given
gamma = 1.4; R = 287; r = 8.4; L = 12; S = 8; V0 = 50; b = 9; T1 = 300;
theta_s = 3*pi/2; q_in = 2.8*10^6; Tw = 300;E = 8/(2*12); M0 = 58.8269;

theta2 = pi:3*pi;
% given values for the span of the crank angle
tspan = [pi 3*pi];

% initial values that are given for the pressure and the work
y0 = [1.013e5 0];

% ODE function to solve for the pressure and work with respect to theta
sol = ode45(@Q1,tspan,y0);

% initial value given for the temperature
y2 = 300;

% ODE function to solve for the temperature with respect to theta
sol2 = ode45(@T,tspan,y2);

% creating the pressure vs. crank angle plot and using subplot so that we can put all 4 of the different graphs on one screen
subplot(2,2,1)
plot(smooth(sol.y(1,:)),'r<-','MarkerSize',4,'LineWidth',1)
grid on
xlabel('Crank Angle (rad)')
ylabel('Pressure (Pa)')
title('Pressure vs. Crank Angle')

% creating the work vs. crank angle plot
subplot(2,2,2)
plot(smooth(sol.y(2,:)),'b<-','MarkerSize',4,'LineWidth',1)
grid on
xlabel('Crank Angle (rad)')
ylabel('Work (J)')
title('Work vs. Crank Angle')

% creating the volume vs. crank angle plot
subplot(2,2,3)
fplot(@V,tspan,'g<-','MarkerSize',4,'LineWidth',1)
grid on
xlabel('Crank Angle (rad)')
ylabel('Volume (cm^3)')
title('Volume vs. Crank Angle')

% creating the temperature vs. crank angle plot
subplot(2,2,4)
plot(sol2.y(1,:),'m<-','MarkerSize',4,'LineWidth',1)
grid on
xlabel('Crank Angle (rad)')
ylabel('Temperature (K)')
title('Temperature vs. Crank Angle')


%% function to solve for the change in pressure and work
function [dP] = Q1(theta,p)
    gamma = 1.4; 
    dP = [-gamma*(p(1)*dV(theta))/V(theta);
          p(1)*dV(theta)];
end

%% function to solve for the change in volume
function deltaV = dV(theta)
    deltaV = 185*sin(theta) + (185*cos(theta)*sin(theta))/(3*(1 - sin(theta)^2/9)^(1/2));
end

%% function to solve for the change in temperature
function [Temp] = T(theta,p)
    R = 287; M0 = 58.8269;
    Temp = (p(1)*V(theta))/((M0)*(R));
end

%% function to solve for the volume with respect to theta
function Volume = V(theta)
    r = 8.4;  V0 = 50; E = 8/(2*12);
    Volume = V0*(1+((r-1)/(2*E))*(1+E*(1-cos(theta))-sqrt(1-E^2*(sin(theta))^2)));
end


