close all; clear; clc;

%% Forward Kinematics
Tfinal = 2;    % Time range (sec)
T = 0.001;      % Sampling period (sec)
t = 0:T:Tfinal; % Time
N = length(t);   % Number of data

L1 = 0.5; % Length of pendulum
L2 = 0.4;
L3 = 0.3;
r = 0.2;

P = zeros(N,2); % EE position
P_dot = zeros(N,2); % EE velocity
theta = zeros(N,3); % Angle

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

for k = 1:N
    if k < N*(1/4)
        xr(k) = 0.75;
        yr(k) = (0.85-0.25)/(N/4)*(k)+0.25;
    elseif k < N*(2/4)
        xr(k) = (0.25-0.75)/(N/4)*(k-N*(1/4))+0.75;
        yr(k) = (0.75-0.85)/(N/4)*(k-N*(1/4))+0.85;
    elseif k < N*(3/4)
        xr(k) = 0.25;
        yr(k) = (0.5-0.75)/(N/4)*(k-N*(2/4))+0.75;
    elseif k <= N*(4/4)
        xr(k) = (0.75-0.25)/(N/4)*(k-N*(3/4))+0.25;
        yr(k) = (0.25-0.5)/(N/4)*(k-N*(3/4))+0.5;
    end
    P(k,:) = [xr(k)
              yr(k)];
end
for k = 1:N-1
   P_dot(k+1,:) = (P(k+1,:) - P(k,:))/T; %% EE derivative
end

theta(1,:) = [330, 60.2267 65.2493] * pi/180;

%% Jacobian Calculation

J = zeros(N, 9);
Jacobian = zeros(3,3);

for k = 1:N-1
   J(k,1) = - L1*sin(theta(k,1)) - L2*sin(theta(k,1) + theta(k,2)) - L3*sin(theta(k,1) + theta(k,2) + theta(k,3));
   J(k,2) = - L2*sin(theta(k,1) + theta(k,2)) - L3*sin(theta(k,1) + theta(k,2) + theta(k,3));
   J(k,3) = - L3*sin(theta(k,1) + theta(k,2) + theta(k,3));
   J(k,4) = L1*cos(theta(k,1)) + L2*cos(theta(k,1) + theta(k,2)) + L3*cos(theta(k,1) + theta(k,2) + theta(k,3));
   J(k,5) = L2*cos(theta(k,1) + theta(k,2)) + L3*cos(theta(k,1) + theta(k,2) + theta(k,3));
   J(k,6) = L3*cos(theta(k,1) + theta(k,2) + theta(k,3));
   
   Jacobian = [J(k,1) J(k,2) J(k,3);
               J(k,4) J(k,5) J(k,6)];
   invJacobian = pinv(Jacobian); % Apply pseudo inverse matrix
   
   theta(k+1,:) = theta(k,:) + (invJacobian * P_dot(k,:)')'*T;
    
    
end

%% Forward Kinematics

x = zeros(N,4); y = zeros(N,4);
for k = 1:N
   x(k,1) = 0; % Joint 0
   y(k,1) = 0; % Joint 0
   x(k,2) = L1*cos(theta(k,1)); % Joint 1
   y(k,2) = L1*sin(theta(k,1)); % Joint 1
   x(k,3) = L1*cos(theta(k,1)) + L2*cos(theta(k,1) + theta(k,2)); % Joint 2
   y(k,3) = L1*sin(theta(k,1)) + L2*sin(theta(k,1) + theta(k,2)); % Joint 2
   x(k,4) = L1*cos(theta(k,1)) + L2*cos(theta(k,1) + theta(k,2)) + L3*cos(theta(k,1) + theta(k,2) + theta(k,3)); % Joint 3
   y(k,4) = L1*sin(theta(k,1)) + L2*sin(theta(k,1) + theta(k,2)) + L3*sin(theta(k,1) + theta(k,2) + theta(k,3)); % Joint 3
end

%% Simulation

figure('color', 'w');
for k = 1:10:N
    plot(P(:,1),P(:,2),'y', 'linewidth',4); hold on; % Predefined trajectory
    plot(x(k,1),y(k,1),'ko'); hold on; % Joint 0
    plot(x(k,2),y(k,2),'bo'); hold on; % Joint 1
    plot(x(k,3),y(k,3),'rs'); hold on; % Joint 2
    plot(x(k,4),y(k,4),'gs'); hold on; % Joint 3
    plot([x(k,1) x(k,2)], [y(k,1) y(k,2)], 'b', 'linewidth', 2); hold on;
    plot([x(k,2) x(k,3)], [y(k,2) y(k,3)], 'r', 'linewidth', 2); hold on;
    plot([x(k,3) x(k,4)], [y(k,3) y(k,4)], 'g', 'linewidth', 2); hold on;
    plot(x(1:100:k,4), y(1:100:k,4), 'b', 'marker', '.', 'markersize', 10); % End effector trajectory
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
plot(t, theta(:,3)*180/pi, 'g', 'linewidth', 2); hold on;
legend('\theta_1','\theta_2','\theta_3');
ylabel('Angle (Deg)'); xlabel('time (sec)');
axis([0 2 0 360]);

subplot(122);
plot(t,x(:,4), 'b', 'linewidth', 2); hold on;
plot(t,y(:,4), 'r', 'linewidth', 2); hold on;
legend('x_3', 'y_3');
ylabel('Distance(m)'); xlabel('time (sec)');
axis([0 2 -0.5 1]);