%% Assignment 2 - LOS Guidance:
%%

clc; clear; close all;

%% Constants:

VP = 200;
VT = 0.6*VP;
K = 1;

%% Initial Conditions:

R0 = 7000;%sqrt(10000^2 + 1000^2);
theta0 = pi/6; %atan(0.1);

alphaT0 = 17*pi/18;
alphaP0 = theta0;

VR0 = VT*cos(alphaT0-theta0) - VP*cos(alphaP0-theta0);
Vtheta0 = VT*sin(alphaT0-theta0) - VP*sin(alphaP0-theta0);
% aP0 = VP*Vtheta0/R0;

xp0 = 0;
yp0 = 0;

xt0 = R0*cos(theta0);
yt0 = R0*sin(theta0);

RP0 = sqrt(xt0^2 + yt0^2);


%% tspan:

tspan = linspace(0, 500, 5000);

%% Options: Event Terminal

options = odeset('Events', @(t, y) event_terminate(t, y));

%% ODE45:

%y0 = [xP, yP, alphaP, thetaP, RP, xT, yT, thetaT, RT]';
%y0 = [xP, yP, alphaP, RP, xT, yT, RT]';
y0 = [xp0, yp0, theta0, RP0, 10000, 1000, R0]';

%y0 = [R0, theta0, alphaP0, alphaT0, xp0, yp0, xt0, yt0, VR0, Vtheta0]';
y20 = [R0, theta0, theta0, alphaT0, xp0, yp0, xt0, yt0, VR0, Vtheta0]';

%[t, y] = ode45(@(t, y)BR_los(t, y, VT, VP, K, alphaT0), tspan, y0, options);
[t, y] = ode45(@(t, y)los_BR(t, y, VT, VP, 0, K), tspan, y20, options);

R = y(:, 1);
theta = y(:, 2);
VR = y(:, 9);
Vtheta = y(:, 10);
%aP = y(:, 5);
alphaT = y(:, 4);
alphaP = y(:, 3);

xP = y(:, 5);
yP = y(:, 6);

xT = y(:, 7);
yT = y(:, 8);

thetaP = atan2(yP, xP);
thetaP(1) = (theta0);
thetaT = atan2(yT, xT);

Rp = sqrt(xP.^2 + yP.^2);
aP = K*VP*Rp.*(thetaT - thetaP);

% xP = y(:, 1);
% yP = y(:, 2);
% alphaP = y(:, 3);
% RP = y(:, 4);
% xT = y(:, 5);
% yT = y(:, 6);
% RT = y(:, 7);

%% Plotting Block:

% t-R:
figure(1)
plot(t, R);
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('R - t')
xlabel("t -->")
ylabel("R -->")
axis padded
%saveas(gcf, 'Shear Flow.png

% t-theta:
figure(2)
plot(t, theta)
hold on
plot(t, alphaP)
plot(t, alphaP-theta)
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('\theta - t')
xlabel("t -->")
ylabel("\theta -->")
axis padded
%saveas(gcf, 'Shear Flow.png')

% t-Vtheta:
figure(3)
plot(t, Vtheta);
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('V_\theta - t')
xlabel("t -->")
ylabel("V_\theta -->")
axis padded

% t-VR:
figure(4)
plot(t, VR);
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('V_R - t')
xlabel("t -->")
ylabel("V_R -->")
axis padded

% t-aP:
figure(5)
plot(t, aP);
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('a_P - t')
xlabel("t -->")
ylabel("a_P -->")
axis padded

% P, T:
figure(6)
plot(xP, yP, 'DisplayName', 'Pursuer');
hold on
plot(xT, yT, 'DisplayName', 'Target');
%axis equal
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('Trajectory')
xlabel("x -->")
ylabel("y -->")
axis padded


%% ode_function:

%y0 = [R0, theta0, alphaP0, alphaT0, xp0, yp0, xt0, yt0, VR0, Vtheta0]';

function dydt = los_BR(t, y, VT, VP, cT, K)

    dydt = zeros(length(y), 1);

    dydt(1) = y(9);
    dydt(2) = y(10)/y(1);
    
    Rp = sqrt(y(5)^2 + y(6)^2);
    thetaT = atan2(y(8), y(7));

    if y(2)==0 && y(1)==0
        thetaP = thetaT;
    else
        thetaP = atan2(y(6), y(5));
    end

    dydt(3) = K*Rp*(thetaT - thetaP);
    dydt(4) = -3/VT; %0*cT*dydt(2);

    dydt(5) = VP*cos(y(3));
    dydt(6) = VP*sin(y(3));
    
    dydt(7) = VT*cos(y(4));
    dydt(8) = VT*sin(y(4));

    dydt(9) = -(dydt(4) - dydt(2))*VT*sin(y(4)-y(2)) + (dydt(3) - dydt(2))*VP*sin(y(3)-y(2));
    dydt(10) = (dydt(4) - dydt(2))*VT*cos(y(4)-y(2)) - (dydt(3) - dydt(2))*VP*cos(y(3)-y(2));
    
end
