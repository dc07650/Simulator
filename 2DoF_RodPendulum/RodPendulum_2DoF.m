% dynamic simulation
clear all; clc; close all;
T=0.001;            % Sampling Period [sec]
Tfinal=10;          % Simulation Time [sec]
t=0:T:Tfinal;       % Time matrix [sec]
N=length(t);        % Sample Units [Samples]
L1 = 0.4; L2 = 0.4; % Rod Length [m]
m1 = 5;  m2 = 3;   % Mass of Rod [kg]
g = 9.81;           % Gravitational Acceleration [m/s^2]

%% Generalized coordinate
% th: Angle [rad], dth: Angular Velocity [rad/s], ddth: Angular Acceleration [rad/s^2]
th = zeros(N,2); 
dth = zeros(N,2); 
ddth = zeros(N,2);
tau1 = zeros(N,1); 
tau2 = zeros(N,1); 
tau = zeros(N,2);

th(1,1) = 45*pi/180;  th(1,2) = 30*pi/180; % Initial Values

Tfinal_torque = 0;
for k=1:N
    % Input Torque
    if k<round((Tfinal_torque)/(Tfinal)*N,0)
        tau1(k)  = 5*sin(pi*t(k));
        tau2(k)  = 2*(1 - cos(pi*t(k)));
        tau(k,:) = [tau1(k) tau2(k)];
    end
end

%% Simulation
for k=1:N-1
    % Inertia
    M11(k) = 1/4*m1*L1^2 + 1/12*m1*L1^2 + m2*L1^2;
    M12(k) = 1/4*m2*L1*L2*cos(th(k,1) - th(k,2));
    M21(k) = M12(k);
    M22(k) = 1/4*m2*L2^2;
    M = [M11(k) M12(k); M21(k) M22(k)];
    M_inv = inv(M);   
    % Centrifugal and Coriolois Force
    C1(k) = 1/4*m2*L1*L2*sin(th(k,1) - th(k,2))*(dth(k,2))^2;
    C2(k) = -1/4*m2*L1*L2*sin(th(k,1) - th(k,2))*(dth(k,1))^2;
    C = [C1(k) C2(k)];
    % Gravitation
    G1(k) = 1/2*m1*g*L1*sin(th(k,1)) + m2*g*L1*sin(th(k,1));
    G2(k) = 1/2*m2*g*L2*sin(th(k,2));
    G = [G1(k) G2(k)];
    % Forward Dynamics (Torque realm -> Angular Acceleration)
    C_input = tau(k,:) - C - G;
    ddth(k,:) = (M_inv*C_input');       % Computation of Angular Acceleration of time step k
    dth(k+1,:) = dth(k,:)+ddth(k,:)*T;  % Linear Integration
    th(k+1,:) = th(k,:)+dth(k+1,:)*T;   
end

%% Computation from Angle to Cartesian
x=zeros(N,3); y=zeros(N,3); % Initialization
for k=1:N
    x(k,1) = 0;
    y(k,1) = 0;
    x(k,2) =  L1*sin(th(k,1));
    y(k,2) = -L1*cos(th(k,1));
    x(k,3) = x(k,2) + L2*sin(th(k,2));
    y(k,3) = y(k,2) - L2*cos(th(k,2));
end

%% Animation
figure('color','w');
T_persistant = 0.5;
kp = round(N*T_persistant/Tfinal,0);
for k=1:10:N
    plot(x(k,1),y(k,1),'ko'); hold on; % Joint 1
    plot(x(k,2),y(k,2),'bo'); hold on; % Joint 2
    plot(x(k,3),y(k,3),'rs'); hold on; % Joint 3
    plot([x(k,1) x(k,2)],[y(k,1) y(k,2)],'b','linewidth',2); hold on; % Link 1
    plot([x(k,2) x(k,3)],[y(k,2) y(k,3)],'r','linewidth',2); hold on; % Link 2
    if k < kp
        plot(x(1:k,3),y(1:k,3),'c.');            % End Point
    else
        plot(x((1+k-kp):k,3),y((1+k-kp):k,3),'c.');            % End Point
    end
    axis([-1 1 -1 1])                         % Axis X, Y
    set(gca,'DataAspectRatio',[1 1 1])
    grid on;
    drawnow;
    hold off;
end
figure('color','w');
subplot(211); % Graphs of end point
plot(t,tau1,'b','linewidth',2); hold on;
plot(t,tau2,'r','linewidth',2); hold on;
legend('\tau_1','\tau_2')
ylabel('Torque(Nm)'); xlabel('time(sec)')
subplot(212); % Graphs of joint angle
plot(t,th(:,1)*180/pi,'b','linewidth',2); hold on;
plot(t,th(:,2)*180/pi,'r','linewidth',2); hold on;
legend('\theta_1','\theta_2')
ylabel('Angle(deg)'); xlabel('time(sec)')
figure('color','w');
subplot(311); % Graphs of end point
plot(t,x(:,3),'b','linewidth',2); hold on;
ylabel('x3(m)'); xlabel('time(sec)')
subplot(312); % Graphs of joint angle
plot(t,y(:,3),'r','linewidth',2); hold on;
ylabel('y3(m)'); xlabel('time(sec)')
subplot(313); % Graphs of joint angle
plot(x(:,3),y(:,3),'b','linewidth',2); hold on;
ylabel('y3(m)'); xlabel('x3(m)')