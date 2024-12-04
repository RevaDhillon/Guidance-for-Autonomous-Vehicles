%% Circular Trajectory:

clc; clear; close all;

%% Constants:

V = 70.0;

xP0 = -100;
yP0 = 0;
xR0 = 498-498*(cos(0.001));
yR0 = 498*sin(0.001);
xC = [498, 1026, 674, 1061, 1449];
yC = [0, 530, 886, 496, 886];

theta_list = [39, 140, 125, 37]*pi/180;%

xC0 = xC(1);
yC0 = yC(1);

R = sqrt((xR0-xC0)^2 + (yR0-yC0)^2);
L1 = sqrt((xR0-xP0)^2 + (yR0-yP0)^2);

theta0 = atan2(yR0-yC0, xR0-xC0);
alpha0 = atan2(yR0-yP0, xR0-xP0);
alphaT0 = theta0 - pi/2;

alphaP0 = 80*pi/180;
VT0 = V*cos(alphaP0 - alpha0)/sin(theta0 - alpha0);

%% Initial Conditions:

% alpha0 = 0*pi/180;
% 
% eta0 = 0*pi/180;
% beta0 = 0*pi/180;
% 
% VT0 = V*cos(eta0)/cos(beta0);
% 
% alphaT0 = 90*pi/180;
% alphaP0 = eta0 + alphaT0 - beta0;
% 
% xR0 = R*cos(alpha0);
% yR0 = R*sin(alpha0);
% xP0= xR0 - L1*cos(alphaP0 - eta0);
% yP0 = yR0 - L1*sin(alphaP0 - eta0);

%% tspan:

tspan = linspace(0, 100, 500*3);

%% Options: Event Terminal

%options = odeset('Events', @(t, y) event_terminate(t, y));

%% ODE45:
%y0 = [x, y, xR, yR, alphaP, VT]';

y0 = [xP0, yP0, xR0, yR0, alphaP0]';

options = odeset('RelTol',1e-9,'AbsTol',1e-12);

[t, y] = ode45(@(t, y)ode_complete(t, y, V, L1, xC, yC, theta_list), tspan, y0, options);

xP_sol = y(:, 1);
yP_sol = y(:, 2);
xR_sol = y(:, 3);
yR_sol = y(:, 4);
alphaP_sol = y(:, 5);
% theta = y(:, 6);
% alphaP = y(:, 7);
% alphaT = y(:, 8);
%VT_sol = y(:, 6);


%% Plotting Block:

%L1, t:


figure;
hold on;
grid on;
xlabel('x (m)');
ylabel('y (m)');
% Initialize the plots for the animation
h1 = plot(xP_sol(1), yP_sol(1), 'bo', 'MarkerSize', 8, 'DisplayName', 'Pursuer Trajectory'); % Moving point for the first trajectory
path1 = plot(xP_sol(1), yP_sol(1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Pursuer Path'); % Path line for the first trajectory
h2 = plot(xR_sol(1), yR_sol(1), 'ro', 'MarkerSize', 8, 'DisplayName', 'Target Trajectory'); % Moving point for the second trajectory
path2 = plot(xR_sol(1), yR_sol(1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Target Path'); % Path line for the second trajectory
legend('show');
% Determine the minimum length of time vectors for the loop
minLength = min(length(xP_sol), length(xR_sol));
% Animate the results for both trajectories
for i = 1:minLength
    % Update the moving point's position for the first trajectory
    set(h1, 'XData', xP_sol(i), 'YData', yP_sol(i));
    % Update the path line data for the first trajectory
    set(path1, 'XData', xP_sol(1:i), 'YData', yP_sol(1:i));

    % Update the moving point's position for the second trajectory
    set(h2, 'XData', xR_sol(i), 'YData', yR_sol(i));
    % Update the path line data for the second trajectory
    set(path2, 'XData', xR_sol(1:i), 'YData', yR_sol(1:i));

    % Update the title to display the current time
    % title(['Time = ', num2str(t1(i), '%.2f'), ' s']);

    % Pause for a short duration to control the speed of the animation
    pause(0.0001); % Adjust this value to speed up or slow down the animation
end
axis("equal")
hold off;


% figure(1)
% plot(tspan, sqrt((xR_sol-xP_sol).^2 + (yR_sol-yP_sol).^2), 'r-','DisplayName', 'L1');
% axis equal
% legend("show")
% ax = gca;
% ax.FontSize = 16;
% grid on
% title('Trajectory')
% xlabel("x -->")
% ylabel("y -->")
% axis padded
% 
% %P, T:
% figure(2)
% plot(xP_sol(1), yP_sol(1), 'r-','DisplayName', 'Pursuer');
% hold on
% plot(xR_sol(1), yR_sol(1), 'k-', 'DisplayName', 'Target');
% for i=2:1:length(xP_sol)
%     scatter(xP_sol(i), yP_sol(i), 'r', 'HandleVisibility','off');
%     scatter(xR_sol(i), yR_sol(i), 'k', 'HandleVisibility','off');
%     %plot([xP_sol(i), xR_sol(i)], [yP_sol(i), yR_sol(i)])
%     pause(0.00001)
% end
% axis equal
% legend("show")
% ax = gca;
% ax.FontSize = 16;
% grid on
% title('Trajectory')
% xlabel("x -->")
% ylabel("y -->")
% axis padded


%% Circular trajectory:
%y0 = [x, y, xR1, yR1, xR2, yR2, xR3, yR3 alphaP, VT]';

function dydt = ode_complete(t, y, V, L1, xC, yC, theta_list)

    dydt = zeros(length(y), 1);
    
    persistent theta %= atan2(y(4)-yC(1), y(3)-xC(1));

    if isempty(theta)
        theta = atan2(y(4)-yC(1), y(3)-xC(1));
    end

    persistent xCt

    if isempty(xCt)
        xCt = xC(1);
    end

    persistent yCt % = yC(1);

    if isempty(yCt)
        yCt = yC(1);
    end
    

    alpha = atan2(y(4)-y(2), y(3)-y(1));
    alphaT = theta - pi/2;
    eps = 0.1;

    VT = V*cos(y(5) - alpha)/sin(theta - alpha);

    if abs(alphaT - theta_list(1) + pi/2)<eps && xCt == xC(1) % && yCt == yC(1)
        xCt = xC(2);
        yCt = yC(2);
        theta = atan2(y(4)-yC(1), y(3)-xC(1));

    elseif  abs(alphaT - theta_list(2) + pi/2)<eps && xCt == xC(2)  && yCt == yC(2) %y(3)<1350 && y(4) < 800
        xCt = xC(3);
        yCt = yC(3);
        theta = atan2(y(4)-yC(3), y(3) - xC(3));

    elseif abs(alphaT - theta_list(3) + pi/2)<eps && xCt == xC(3)  && yCt == yC(3)
        xCt = xC(4);
        yCt = yC(4);
        theta = atan2(y(4)-yC(4), y(3) - xC(4));

    elseif abs(alphaT - theta_list(4) + pi/2)<eps && xCt == xC(4)  && yCt == yC(4)
        xCt = xC(5);
        yCt = yC(5);
        theta = atan2(y(4)-yC(5), y(3) - xC(5));
    else 
        theta = atan2(y(4)-yCt, y(3) - xCt);
    
    end
    %disp(theta)

    dydt(1) = V*cos(y(5));
    dydt(2) = V*sin(y(5));

    dydt(4) = VT*sin(alphaT);
    dydt(3) = VT*cos(alphaT);

    dydt(5) = -2*V*sin(y(5)-alpha)/L1;

end