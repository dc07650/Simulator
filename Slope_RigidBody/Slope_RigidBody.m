%% Lagrangian Mechanics MATLAB
% Author: Won Bin Choi (B4)
% E-mail: dc07650@naver.com
% Organization: Sogang University(Korea, Republic of), Mechanical Engineering
% Date: April 22, 2021
clear all; close all; clc;

%% Initilize
% Time Variables
Tfinal = 20;     % Simulation time length [sec]
T = 0.001;      % Sampling time [sec]
t = 0:T:Tfinal; % Time vector
N = length(t);
% Space Variables
x = zeros(1,N);     % Displacement [m]
dx = zeros(1,N);    % Speed [m/s]
ddx = zeros(1,N);   % Acceleration [m/s^2]
ph = zeros(1,N);    % Angular Displacement [deg]
dph = zeros(1,N);   % Angular Speed [deg/s]
ddph = zeros(1,N);  % Angular Acceleration [deg/s^2]
% Instance Variables
R = 5;           % Body radius [m]
m = 1;             % Body mass [m]
Ig = 1/2*m*R^2;    % Mass Moment of Inertia [kg*m^2]
theta = 45/180*pi; % Slope angle [deg -> rad]
% Environment Variables
g = 9.81;          % Gravitational Acceleration [m/s^2]
L = 100;            % Slope length [m]
tau = zeros(1,N);
%% Simulation
% Torque Emulation
Tfinal_torque = 0;
for k=1:N
    if k<round((Tfinal_torque)/(Tfinal)*N,0)
        tau(k)  = 20*(1 - cos(pi*t(k)));
    end
end
% Equation of Motion
for k = 1:N-1
    ddph(k) = (tau(k) - m*g*sin(theta)*R)/(m*R^2 + Ig);
    dph(k+1) = dph(k) + ddph(k)*T;
    ph(k+1) = ph(k) + dph(k)*T;
end

%% Coordinate Transformation
x_s = zeros(1,N); y_s = zeros(1,N); % Slope coordinate initilization
for k = 1:N
    x(k) = - R*ph(k);
    dx(k) = - R*dph(k);
    ddx(k) = - R*ddph(k);   
    x_s(k) = R*sin(theta) + x(k)*cos(theta);
    y_s(k) = L*sin(theta) + R*cos(theta) - x(k)*sin(theta);
end

%% Simulation
SimulationSpeed = 50;
for k = 1:SimulationSpeed:N
    if y_s(k) < R
       break; 
    end
    % Body Plot
    plot([0 L*cos(theta)], [L*sin(theta) 0], 'r', 'linewidth', 2); hold on;
    plot([L*cos(theta) L*cos(theta)+2*R], [0 0], 'r', 'linewidth', 2); hold on;
    plot([-R 0], [L*sin(theta) L*sin(theta)], 'r', 'linewidth', 2); hold on;
    circle(x_s(k),y_s(k),R); hold on;
    plot(x_s(k), y_s(k), 'ro', 'MarkerSize', 2); hold on;
    plot(x_s(k) + R*cos(ph(k)), y_s(k) + R*sin(ph(k)), 'bo'); hold on;
    axis([-R L*cos(theta)+2*R -1 2*R+L*sin(theta)])                         % Axis X, Y
    set(gca,'DataAspectRatio',[1 1 1])
    grid on;
    drawnow;
    hold off;
end

%% Plot
figure('color','w');

subplot(211); % Graphs of end point
plot(t,tau,'b','linewidth',2); hold on;
legend('\tau')
ylabel('Torque [N*m]'); xlabel('Time [sec]')

subplot(212); % Graphs of displacement
plot(t,ph(:)*180/pi,'r','linewidth',2); hold on;
plot(t,dph(:)*180/pi,'g','linewidth',2); hold on;
plot(t,ddph(:)*180/pi,'b','linewidth',2); hold on;
legend('\phi [deg]','\omega [deg/s]', '\alpha [deg/s^2]')
ylabel('Values'); xlabel('Time [sec]')

figure('color','w');

subplot(311); % Graphs of end point
plot(t,x_s(:),'b','linewidth',2); hold on;
ylabel('x [m]'); xlabel('Time [sec]')

subplot(312); % Graphs of joint angle
plot(t,y_s(:),'r','linewidth',2); hold on;
ylabel('y [m]'); xlabel('Time [sec]')

subplot(313); % Graphs of joint angle
plot(x_s(:), y_s(:),'b','linewidth',2); hold on;
ylabel('y [m]'); xlabel('x [m]')

%% Functions
function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang = 0:0.01:2*pi;
xp = r*cos(ang);
yp = r*sin(ang);
plot(x+xp,y+yp, 'g', 'Linewidth', 2);
end
