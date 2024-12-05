%% Main

R0 = 7000;

bsx = 0;
bsy = 0;

theta0 = 30*pi/180;
alphap0 = theta0;

xp0 = bsx;
yp0 = bsy;

alphat0 = 170*pi/180;
xt0 = xp0 + R0*cos(theta0);
yt0 = yp0 + R0*sin(theta0);

nu = 0.6;
Vp0 = 200;
Vt0 = nu*Vp0;

A_target = 0*-3;

X0 = [R0,theta0,alphap0,xp0,yp0,alphat0,xt0,yt0,Vp0,Vt0]';

options = odeset('Events', @(t, X) event_terminal(t, X));
tspan = linspace(0,500,1000*500);
[t,X] = ode45(@(t,X)system(t,X,[bsx,bsy],A_target),tspan,X0,options);
R = X(:,1);
theta = X(:,2);
alphap = X(:,3);
xp = X(:,4);
yp = X(:,5);
alphat = X(:,6);
xt = X(:,7);
yt = X(:,8);
A_cmd = zeros(1,length(t));

for i = 1:length(t)
    % Vr = X(10)*cos(X(6)-(X(2))) - X(9)*cos(X(3)-X(2));
    thetam = atan2(yp(i)-bsy,xp(i)-bsx);
    thetat = atan2(yt(i)-bsy,xt(i)-bsx);
    Rm = norm([yp(i)-bsy,xp(i)-bsx]);
    A_cmd(i) = 1*Rm*(thetat-thetam);
end

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

figure(7); clf
title('A_{cmd}')
xlabel('t')
ylabel('A_{cmd}')
hold on
legend('on')
plot(t,A_cmd,'b',LineWidth=2,DisplayName='Latax')

%% System Dynamics: X = [1_R,2_theta,3_alphap,4_xp,5_yp,6_alphat,7_xt,8_yt,9_Vp,10_Vt]

function dX = system(t,X,bs,A_target)

    % Vr = X(10)*cos(X(6)-(X(2))) - X(9)*cos(X(3)-X(2));
    thetam = atan2(X(5)-bs(2),X(4)-bs(1));
    thetat = atan2(X(8)-bs(2),X(7)-bs(1));
    Rm = norm([X(5)-bs(2),X(4)-bs(1)]);
    
    dX = zeros(length(X),1);
    dX(1) = X(10)*cos(X(6)-X(2)) - X(9)*cos(X(3)-X(2));
    dX(2) = (1/X(1))*(X(10)*sin(X(6)-X(2)) - X(9)*sin(X(3)-X(2)));
    A_cmd = 1*Rm*(thetat-thetam);
    dX(3) = A_cmd/X(9);
    dX(4) = X(9)*cos(X(3));
    dX(5) = X(9)*sin(X(3));
    dX(6) = A_target/X(10); 
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