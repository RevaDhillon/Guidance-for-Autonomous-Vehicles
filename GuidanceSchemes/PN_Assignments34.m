%% Assignment 3,4 - PN:

clc; clear; close all;

%% Constants:

% Both 0.8 -- Important!!!!!!!!!!!!!! -- Tolerance = 0, Zoom in.

%  ν = 0.8, 0.9, 1, 1.5, 3; k = 1,5,10.
nu = 1.2;%input("nu : ");
VT = 240;
VP = 400;

%% Initial Conditions:

R0 = 7000;
theta0 = pi/6;

alphaT0 = 60*pi/180;
alphaP0 = 40*pi/180; %theta0 + delta;

VR0 = VT*cos(alphaT0-theta0) - VP*cos(alphaP0-theta0);
Vtheta0 = VT*sin(alphaT0-theta0) - VP*sin(alphaP0-theta0);

xp0 = 0;
yp0 = 0;

xt0 = R0*cos(theta0);
yt0 = R0*sin(theta0);

N = 3;
phi0 = alphaP0 - N*theta0;


%% tspan:

tspan = linspace(0, 50, 500000);

%% Options: Event Terminal

options = odeset('Events', @(t, y) event_terminate(t, y));

%% ODE45:

y0 = [R0, theta0, alphaP0, alphaT0, xp0, yp0, xt0, yt0, VR0, Vtheta0]';

[t, y] = ode45(@(t, y)ode_PPN_NM(t, y, VT, VP, N, phi0), tspan, y0, options);

R = y(:, 1);
theta = y(:, 2);
VR = y(:, 9);
Vtheta = y(:, 10);
%aP = y(:, 5);
alphaT = y(:, 4);
alphaP = y(:, 3);

xp = y(:, 5);
yp = y(:, 6);

xt = y(:, 7);
yt = y(:, 8);

aP = VP*N*Vtheta./R;

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
% plot(t, alphaP)
% plot(t, alphaP-theta)
% %axis equal
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
plot(xp, yp, 'DisplayName', 'Pursuer');
hold on
plot(xt, yt, 'DisplayName', 'Target');
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

function dydt = ode_PPN_NM(t, y, VT, VP, N, phi0)

    dydt = zeros(length(y), 1);
    k = N-1;

    dydt(1) = y(9);
    dydt(2) = y(10)/y(1);

    dydt(3) = N*dydt(2);
    dydt(4) = 0;%cT*dydt(2);

    dydt(5) = VP*cos(y(3));
    dydt(6) = VP*sin(y(3));
    
    dydt(7) = VT*cos(y(4));
    dydt(8) = VT*sin(y(4));

    dydt(9) = dydt(2)*(y(10) + N*VP*sin(k*y(2) + phi0));
    dydt(10) = -dydt(2)*(y(9) + N*VP*cos(k*y(2) + phi0));
    
end

%% ode_function:

%y0 = [R0, theta0, alphaP0, alphaT0, xp0, yp0, xt0, yt0, VR0, Vtheta0]';

function dydt = ode_PPN_M(t, y, VT, VP, N, phi0)

    dydt = zeros(length(y), 1);
    k = N-1;

    dydt(1) = y(9);
    dydt(2) = y(10)/y(1);

    dydt(3) = N*dydt(2);
    dydt(4) = 0;%cT*dydt(2);

    dydt(5) = VP*cos(y(3));
    dydt(6) = VP*sin(y(3));
    
    dydt(7) = VT*cos(y(4));
    dydt(8) = VT*sin(y(4));

    dydt(9) = dydt(2)*(y(10) + N*VP*sin(k*y(2) + phi0));
    dydt(10) = -dydt(2)*(y(9) + N*VP*cos(k*y(2) + phi0));
    
end