%% Main

R0 = 3000;
theta0 = 30*pi/180;
delta0 = 0*pi/180;
delta_des = 0*pi/180;
alphap0 = theta0+delta0;
xp0 = 0;
yp0 = 0;
alphat0 = 170*pi/180;
xt0 = xp0 + R0*cos(theta0);
yt0 = yp0 + R0*sin(theta0);

nu = 1.2;
Vt0 = 200;
Vp0 = nu*Vt0;

X0 = [R0,theta0,alphap0,xp0,yp0,alphat0,xt0,yt0,Vp0,Vt0]';

options = odeset('Events', @(t, X) event_terminal(t, X));
tspan = linspace(0,500,1000*500);
[t,X] = ode45(@(t,X)system(t,X,delta_des),tspan,X0,options);
R = X(:,1);
theta = X(:,2);
alphap = X(:,3);
xp = X(:,4);
yp = X(:,5);
alphat = X(:,6);
xt = X(:,7);
yt = X(:,8);
A_latax = (1./R).*(Vt0.*sin(alphat-theta) - Vp0.*sin(alphap-theta))*Vp0 - 10*(alphap-theta-delta_des);

Vr = Vt0*cos(alphat-theta) - Vp0*cos(alphap-theta);
Vtheta = Vt0*sin(alphat-theta) - Vp0*sin(alphap-theta);
%% Plots

figure(1); clf;
title('Trajectory')
xlabel('X')
ylabel('Y')
hold on
legend('on')
plot(xp,yp,'b',LineWidth=2,DisplayName='Pursuer')
plot(xt,yt,'r',LineWidth=2,DisplayName='Target')

figure(2); clf;
title('Heading Angles')
xlabel('t')
ylabel('Heading Angles')
hold on
legend('on')
plot(t,theta*180/pi,'b',LineWidth=2,DisplayName='\theta')
plot(t,alphap*180/pi,'r',LineWidth=1,DisplayName='\alpha_P')
plot(t,(alphap-theta)*180/pi,'g',LineWidth=1,DisplayName='\delta')

figure(3); clf
title('R')
xlabel('t')
ylabel('R')
hold on
legend('on')
plot(t,R,'b',LineWidth=2,DisplayName='R')

figure(4); clf
% title('V_r V_\theta space')
xlabel('V_\theta')
ylabel('V_r')
hold on
plot(Vtheta,Vr,'b',LineWidth=2)

figure(5); clf
title('R vs \psi')
xlabel('\psi')
ylabel('R')
hold on
legend('on')
plot(alphat-theta,R,'b',LineWidth=2,DisplayName='R')

figure(6); clf;
title('\delta')
xlabel('t')
ylabel('\delta')
hold on
legend('on')
plot(t,(alphap-theta)*180/pi,'g',LineWidth=1,DisplayName='\delta')

figure(7); clf;
title('Latax')
xlabel('t')
ylabel('Latax')
hold on
legend('on')
plot(t,A_latax,'r',LineWidth=1,DisplayName='Latax')

%% System Dynamics: X = [1_R,2_theta,3_alphap,4_xp,5_yp,6_alphat,7_xt,8_yt,9_Vp,10_Vt]

function dX = system(t,X,delta_des)
    % delta_des = 8*pi/180;
    dX = zeros(length(X),1);
    dX(1) = X(10)*cos(X(6)-X(2)) - X(9)*cos(X(3)-X(2));
    dX(2) = (1/X(1))*(X(10)*sin(X(6)-X(2)) - X(9)*sin(X(3)-X(2)));
    A_latax = dX(2)*X(9) - 10*(X(3)-X(2)-delta_des); % Thetadot*Vp for PP
    dX(3) = A_latax/X(9); % PP, alphap = theta
    dX(4) = X(9)*cos(X(3));
    dX(5) = X(9)*sin(X(3));
    dX(6) = 0*dX(2); % For Non Maneuvering, dalphat = 0
    dX(7) = X(10)*cos(X(6));
    dX(8) = X(10)*sin(X(6));
    dX(9) = 0;
    dX(10) = 0;
end

%% Terminating Condition

function [value, isterminal, direction] = event_terminal(t, X)
    value = X(1)-50;
    isterminal = 1;
    direction = 0;
end