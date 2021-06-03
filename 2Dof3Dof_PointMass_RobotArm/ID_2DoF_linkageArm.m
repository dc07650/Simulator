close all; clear; clc;

%% Forward Kinematics
Tfinal = 2;    % Time range (sec)
T = 0.001;      % Sampling period (sec)
t = 0:T:Tfinal; % Time
N = length(t);   % Number of data

L1 = 0.65; % Length of pendulum
L2 = 0.65;
r = 0.2;

P = zeros(N,2); % EE position
P_dot = zeros(N,2); % EE velocity
theta = zeros(N,2); % Angle
dtheta = zeros(N,2); % Angle
ddtheta = zeros(N,2); % Angle
F = zeros(N,2); % Force on EE
m1 = 0;    % Mass of end points
m2 = 2.5;
g = 9.81;

%% EE trajection
% 
% for k = 1:N
%    th(k) = 2*pi*(k-1)/N; % 0~2*pi
%    P(k,:) = [0.5 + th(k)/15*cos(th(k))
%              0.25 + th(k)/15*sin(th(k))];
% end
% 
% for k = 1:N-1
%    P_dot(k+1,:) = (P(k+1,:) - P(k,:))/T; %% EE derivative
% end
% 
% theta(1,:) = [330, 75.8864, 88.9287] * pi/180;

% for k = 1:N
%     if k < N*(1/4)
%         xr(k) = 0.75;
%         yr(k) = (0.85-0.25)/(N/4)*(k)+0.25;
%     elseif k < N*(2/4)
%         xr(k) = (0.25-0.75)/(N/4)*(k-N*(1/4))+0.75;
%         yr(k) = (0.75-0.85)/(N/4)*(k-N*(1/4))+0.85;
%     elseif k < N*(3/4)
%         xr(k) = 0.25;
%         yr(k) = (0.5-0.75)/(N/4)*(k-N*(2/4))+0.75;
%     elseif k <= N*(4/4)
%         xr(k) = (0.75-0.25)/(N/4)*(k-N*(3/4))+0.25;
%         yr(k) = (0.25-0.5)/(N/4)*(k-N*(3/4))+0.5;
%     end
%     P(k,:) = [xr(k)
%               yr(k)];
% end
% for k = 1:N-1
%    P_dot(k+1,:) = (P(k+1,:) - P(k,:))/T; %% EE derivative
% end

theta(1,:) = [-34.0071, 105.049] * pi/180;

%% Force Generation
for k = 1:N
    if k < N*(1/4)
        F(k,1) = -5;    % Force applied on EE
        F(k,2) = 2.5*g; % Gravity Compensation
    elseif k < N*(2/4)
        F(k,1) = 5;
        F(k,2) = 2.5*g;
    elseif k < N*(3/4)
        F(k,1) = 0;
        F(k,2) = 2.5*g + 5;
    elseif k <= N*(4/4)
        F(k,1) = 0;
        F(k,2) = 2.5*g - 5;
    end
end

%% Jacobian Calculation & Forward Dynamics

J = zeros(N, 4);
Jacobian = zeros(3,3);
tau = zeros(N,2);

for k = 1:N-1
    J(k,1) = - L1*sin(theta(k,1)) - L2*sin(theta(k,1) + theta(k,2));
    J(k,2) = - L2*sin(theta(k,1) + theta(k,2));
    J(k,3) = L1*cos(theta(k,1)) + L2*cos(theta(k,1) + theta(k,2));
    J(k,4) = L2*cos(theta(k,1) + theta(k,2));
    
    Jacobian = [J(k,1) J(k,2);
                J(k,3) J(k,4)];
    
    tau(k,:) = (Jacobian')*F(k,:)';
    
    M = [m2*L2^2 + (m1+m2)*L1^2 + 2*m2*L1*L2*cos(theta(k,2)) , m2*L2^2 + m2*L1*L2*cos(theta(k,2));
        m2*L2^2 + m2*L1*L2*cos(theta(k,2))                  , m2*L2^2];
    V = [-m2*L1*L2*sin(theta(k,2))*dtheta(k,2)^2 - 2*m2*L1*L2*sin(theta(k,2))*dtheta(k,1)*dtheta(k,2),
        m2*L1*L2*sin(theta(k,2))*dtheta(k,1)*dtheta(k,1)];
    G = [m2*L2*g*cos(theta(k,1)+theta(k,2)) + (m1 + m2)*L1*g*cos(theta(k,1))
        m2*L2*g*cos(theta(k,1)+theta(k,2))];
    
    ddtheta(k,:) = inv(M)*(tau(k,:)' - V - G);
    dtheta(k+1,:) = dtheta(k,:) + ddtheta(k,:)*T;
    theta(k+1,:) = theta(k,:) + dtheta(k,:)*T;
end

%% Forward Kinematics

x = zeros(N,3); y = zeros(N,3);
for k = 1:N
   x(k,1) = 0; % Joint 0
   y(k,1) = 0; % Joint 0
   x(k,2) = L1*cos(theta(k,1)); % Joint 1
   y(k,2) = L1*sin(theta(k,1)); % Joint 1
   x(k,3) = L1*cos(theta(k,1)) + L2*cos(theta(k,1) + theta(k,2)); % Joint 2
   y(k,3) = L1*sin(theta(k,1)) + L2*sin(theta(k,1) + theta(k,2)); % Joint 2
end

%% Simulation

figure('color', 'w');
for k = 1:10:N
    plot(P(:,1),P(:,2),'y', 'linewidth',4); hold on; % Predefined trajectory
    plot(x(k,1),y(k,1),'ko'); hold on; % Joint 0
    plot(x(k,2),y(k,2),'bo'); hold on; % Joint 1
    plot(x(k,3),y(k,3),'rs'); hold on; % Joint 2
    plot([x(k,1) x(k,2)], [y(k,1) y(k,2)], 'b', 'linewidth', 2); hold on;
    plot([x(k,2) x(k,3)], [y(k,2) y(k,3)], 'r', 'linewidth', 2); hold on;
    plot(x(1:100:k,3), y(1:100:k,3), 'b', 'marker', '.', 'markersize', 10); % End effector trajectory
    axis([-0.5 1 -0.5 1]);
    grid on;
    drawnow;
    hold off;
end
%% Graph

figure('color', 'w');
subplot(121)
plot(t, theta(:,1)*180/pi, 'b', 'linewidth', 2); hold on;
plot(t, theta(:,2)*180/pi, 'r', 'linewidth', 2); hold on;
legend('\theta_1','\theta_2');
ylabel('Angle (Deg)'); xlabel('time (sec)');
axis([0 2 0 360]);

subplot(122);
plot(t,x(:,3), 'b', 'linewidth', 2); hold on;
plot(t,y(:,3), 'r', 'linewidth', 2); hold on;
legend('x_2', 'y_2');
ylabel('Distance(m)'); xlabel('time (sec)');
axis([0 2 -0.5 1]);

figure('color', 'w');
subplot(121)
plot(t, F(:,1), 'g', 'linewidth', 2); hold on;
plot(t, F(:,2), 'y', 'linewidth', 2); hold on;
legend('F_x','F_y');
ylabel('Force (N)'); xlabel('time (sec)');
low = min(F)*1.2;
low = min(low(1),low(2));
high = max(F)*1.2;
high = max(high(1),high(2));
axis([0 2 low(1) high(1)]);

subplot(122);
plot(t, tau(:,1), 'r', 'linewidth', 2); hold on;
plot(t, tau(:,2), 'b', 'linewidth', 2); hold on;
legend('\tau_1','\tau_2');
ylabel('Torque (Nm)'); xlabel('time (sec)');
low = min(tau)*1.2;
low = min(low(1),low(2));
high = max(tau)*1.2;
high = max(high(1),high(2));
axis([0 2 low high]);