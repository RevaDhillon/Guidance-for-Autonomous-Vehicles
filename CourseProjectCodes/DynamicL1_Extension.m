
clc; close all; clear;

%% Initial Conditions:

xP0 = -10;
yP0 = 1;
alphaP0 = pi/3;
V = 10;
dt = 0.01;
R = 3;

[aP, L1] = return_aPL1(alphaP0, xP0, yP0, V, R);%2*V*V*sin(eta)./L1;
aPmin = min(aP);
L1min = L1(aP == aPmin);

xP = xP0;
yP = yP0;
alphaP = alphaP0;

while L1min(end)>2*R
    [xP1, yP1, alphaP1] = next_point(xP(end), yP(end), V, alphaP(end), dt, min(aP(end, :)));
    [aP2, L12] = return_aPL1(alphaP1, xP1, yP1, V, R);
    aP = [aP; aP2];
    L1 = [L1; L12];
    aPmin = [aPmin, min(aP2)];
    L1min = [L1min, L12(aP2 == min(aP2))];
    xP = [xP, xP1];
    yP = [yP, yP1];
    alphaP = [alphaP, alphaP1];

end

% [xP1, yP1, alphaP1] = next_point(xP0, yP0, V, alphaP0, dt, min(aP));
% [aP2, L12] = return_aPL1(alphaP1, xP1, yP1, V);
% 
% [xP2, yP2, alphaP2] = next_point(xP1, yP1, V, alphaP1, dt, min(aP2));
% [aP3, L13] = return_aPL1(alphaP2, xP2, yP2, V);

%% Plotting Block

figure(1)
plot(L1(1, :), aP(1, :), 'LineWidth', 2, 'DisplayName', "Start")
hold on;
for i=2:1:length(L1min)-1
    if(rem(i, 3)==0)
    plot(L1(i, :), aP(i, :), 'LineWidth', 1, 'HandleVisibility','off')
    end
end
plot(L1(end, :), aP(end, :), 'LineWidth', 2, 'DisplayName', "End")
plot(L1min, aPmin, 'k', 'LineWidth', 2, 'DisplayName', "Minimum a_P envelope")
scatter(L1min, aPmin, 'k', 'filled', 'HandleVisibility','off')
%plot(L12, aP2, 'LineWidth', 2, 'DisplayName', "2")
%plot(L13, aP3, 'LineWidth', 2, 'DisplayName', "3")
legend("show")
title("a_P vs L_1")
xlabel("L_1 (m)")
ylabel("a_P (m/s^2)")
grid on;

figure(2)
plot(L1(1, :), aP(1, :), 'LineWidth', 2, 'DisplayName', "Start")
hold on;
legend("show")
title("a_P vs L_1: Initial")
xlabel("L_1 (m)")
ylabel("a_P (m/s^2)")
grid on;

figure(3)
t_list = 0:1:length(L1min)-1;
t_list = dt*t_list;
plot(t_list, L1min, 'LineWidth', 2, 'DisplayName', "Start")
hold on;
title("L_1 vs t")
ylabel("L_1 (m)")
xlabel("t (s)")
grid on;

figure(4)
plot(t_list, aPmin, 'LineWidth', 2, 'DisplayName', "Start")
hold on;
title("a_P vs t")
ylabel("a_P (m/s^2)")
xlabel("t (s)")
grid on;



%% function for aP:

function [aP, L1] = return_aPL1(alphaP0, xP0, yP0, V, R)
    theta = linspace(-pi, pi, 500);
    
    xR0 = R*(cos(theta));
    yR0 = R*(sin(theta));
    
    eta = alphaP0 - atan2(yR0-yP0, xR0-xP0);
    L1 = sqrt((yR0-yP0).^2 + (xR0-xP0).^2);
    
    aP = 2*V*V*sin(eta)./L1;
end

function [xP1, yP1, alphaP1] = next_point(xP0, yP0, V, alphaP0, dt, aP)
 
    xP1 = xP0 + V*cos(alphaP0)*dt;
    yP1 = yP0 + V*sin(alphaP0)*dt;
    alphaP1 = alphaP0 - (aP)*dt/V;

end

