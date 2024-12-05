clc; clear;

%% Initial Conditions:

R = 10;

alphap0 = 2*pi/4;
Vp = 7;
alphat0 = 2*pi/4;

Xc = 4;
Yc = 4;

xt0 = Xc + R;
yt0 = Yc + 0;

thetadash = -1*pi/4;
thetadashdash = pi/12;

xp0 = xt0; % + L1*cos(thetadash);
yp0 = yt0; % + L1*sin(thetadash);

xt0 = Xc + R*cos(thetadashdash);
yt0 = Yc + R*sin(thetadashdash);
alphat0 = alphat0+thetadashdash;

r = [xt0;yt0] - [xp0;yp0];
alpha0 = atan2(r(2),r(1));
eta0 = alpha0-alphap0;
options = optimoptions('fsolve','Display','none');
Vt0 = fsolve(@(Vt) Vt*cos(alpha0-alphat0) - Vp*cos(eta0),Vp,options);

ExtraPoints = 1;
extra_array = zeros(1,ExtraPoints*3);
extra_thetadashdash = linspace(thetadashdash,0,ExtraPoints+2);
extra_thetadashdash = extra_thetadashdash(2:end-1);
for i = 1:ExtraPoints
    extra_array(3*(i-1)+1) = Xc + R*cos(extra_thetadashdash(i));
    extra_array(3*(i-1)+2) = Xc + R*sin(extra_thetadashdash(i));
    extra_array(3*(i-1)+3) = pi/2+extra_thetadashdash(i);
end


Switch_Matrix = [R, 0.5*R; pi/2, 1.25*pi; 1, -1];

%Switch_Matrix = [R; pi/2; 1];
                 

%% ODE45:
t_end = 5;
tspan = linspace(0, t_end, t_end*100);
y0 = [xp0,yp0,xt0,yt0, extra_array, alphap0,alphat0];
options = odeset('RelTol',1e-15,'AbsTol',1e-15);

[t, y] = ode45(@(t, y)modified_ode(t,y,Vp,R,Switch_Matrix,ExtraPoints),tspan,y0);

xP = y(:, 1);
yP = y(:, 2);
xT = y(:, 3);
yT = y(:, 4);
alphaP = y(:,end-1);
alphaT = y(:,end);

xT1 = y(:, 5);
% xT2 = y(:, 8);
% xT3 = y(:, 11);

yT1 = y(:, 6);
% yT2 = y(:, 9);
% yT3 = y(:, 12);

L = zeros(length(xP),1);
for i = 1:length(xP)
    L(i) = norm([xP(i);yP(i)]-[xT(i);yT(i)]);
end

%% Plotting Block:

figure(1); clf
hold on
title('L Variation')
plot(t,L)


figure(3); clf
hold on;
grid on;
xlabel('Time (s)');
ylabel('Position');
title('Extension')
plot(xP, yP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Pursuer');
plot(xT, yT, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Pseudotarget');
%plot(xT1, yT1, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Extra');
%scatter(Xc+R*cos(0.75*pi),Yc + R*sin(0.75*pi),'filled')

% % Initialize the plots for the animation
% h1 = plot(xP(1), yP(1), 'bo', 'MarkerSize', 8, 'DisplayName', 'Pursuer'); % Moving point for the first trajectory
% path1 = plot(xP(1), yP(1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Pursuer Path'); % Path line for the first trajectory
% 
% h2 = plot(xT(1), yT(1), 'ro', 'MarkerSize', 8, 'DisplayName', 'Pseudotarget'); % Moving point for the second trajectory
% path2 = plot(xT(1), yT(1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Pseudotarget Path'); % Path line for the second trajectory
% 
% h_ = plot(xT1(1), yT1(1), 'ro', 'MarkerSize', 8, 'DisplayName', 'extra');
% h_p= plot(xT1(1), yT1(1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'extra');
% % h__ = plot(xT2(1), yT2(1), 'go', 'MarkerSize', 8, 'DisplayName', 'extra');
% % h__p = plot(xT2(1), yT2(1), 'g-', 'LineWidth', 1.5, 'DisplayName', 'extra');
% % h___ = plot(xT3(1), yT3(1), 'bo', 'MarkerSize', 8, 'DisplayName', 'extra');
% % h___p = plot(xT3(1), yT3(1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'extra');
% legend('show');
% 
% % Determine the minimum length of time vectors for the loop
% minLength = min(length(xP), length(xT));
% axis equal
% 
% % Animate the results for both trajectories
% for i = 1:minLength-1
%     % Update the moving point's position for the first trajectory
%     set(h1, 'XData', xP(i), 'YData', yP(i));
%     % Update the path line data for the first trajectory
%     set(path1, 'XData', xP(1:i), 'YData', yP(1:i));
% 
%     % Update the moving point's position for the second trajectory
%     set(h2, 'XData', xT(i), 'YData', yT(i));
%     % Update the path line data for the second trajectory
%     set(path2, 'XData', xT(1:i), 'YData', yT(1:i));
% 
%     set(h_, 'XData', xT1(i), 'YData', yT1(i));
%     % set(h__, 'XData', xT2(i), 'YData', yT2(i));
%     % set(h___, 'XData', xT3(i), 'YData', yT3(i));
% 
%     set(h_p, 'XData', xT1(1:i), 'YData', yT1(1:i));
%     % set(h__p, 'XData', xT2(1:i), 'YData', yT2(1:i));
%     % set(h___p, 'XData', xT3(1:i), 'YData', yT3(1:i));
% 
%     % Update the title to display the current time
%     % title(['Time = ', num2str(t1(i), '%.2f'), ' s']);
% 
%     % Pause for a short duration to control the speed of the animation
%     pause(0.1*(t(i+1)-t(i))); % Adjust this value to speed up or slow down the animation
% end
% 
% hold off;

%% Functions

function Y = modified_ode(t,y,V,R,Switch_Matrix,ExtraPoints)
% [xp1, yp2, xt3, yt4, extra_array, alphap5, alphat6]

Y = zeros(length(y),1);

% persistent i j;
% 
% if isempty(i)
%     i = 1;
% end
% 
% if isempty(j)
%     j = ones(1,ExtraPoints);
% end

e = 0;
[ai,bi] = size(Switch_Matrix);

% if i < bi
%     if Switch_Matrix(3,i) == 1
%         if y(end)>Switch_Matrix(2,i+1)-e
%             i = i + 1;
%         end
%     else
%         if y(end)<Switch_Matrix(2,i+1)+e
%             i = i + 1;
%         end
%     end
% end

latax_array = zeros(1,ExtraPoints+1);

for check = 1:ExtraPoints

    % if j(check) < bi
    % if Switch_Matrix(3,j(check)) == 1
    %     if y(4 + 3*(check-1)+3)>Switch_Matrix(2,j(check)+1)-e
    %         j(check) = j(check) + 1;
    %     end
    % else
    %     if y(4 + 3*(check-1)+3)<Switch_Matrix(2,j(check)+1)+e
    %         j(check) = j(check) + 1;
    %     end
    % end
    % end

    r = [y(4 + 3*(check-1)+1);y(4 + 3*(check-1)+2)] - [y(1);y(2)];
    alpha = atan2(r(2),r(1));
    eta = alpha-y(end-1);
    % options = optimoptions('fsolve','Display','none');
    Vt = V*cos(eta)/cos(alpha-y(4 + 3*(check-1)+3)); %fsolve(@(Vt) Vt*cos(alpha-y(4 + 3*(check-1)+3)) - V*cos(eta),V,options);
    latax_array(check) = 2*V^2*sin(eta)/norm(r);
    
    Y(4 + 3*(check-1)+1) = Vt*cos(y(4 + 3*(check-1)+3));
    Y(4 + 3*(check-1)+2) = Vt*sin(y(4 + 3*(check-1)+3));
    if t<2.18699
        Y(4 + 3*(check-1)+3) = Vt/R;
    else
        Y(4 + 3*(check-1)+3) = -2*Vt/R;
    end
    
end


r = [y(3);y(4)] - [y(1);y(2)];
alpha = atan2(r(2),r(1));
eta = alpha-y(end-1);
% options = optimoptions('fsolve','Display','none');
Vt = V*cos(eta)/cos(alpha-y(end)); %fsolve(@(Vt) Vt*cos(alpha-y(end)) - V*cos(eta),V,options);
latax_array(end) = 2*V^2*sin(eta)/norm(r);


weights = [50;1];
latax = latax_array*weights/sum(weights);

Y(1) = V*cos(y(end-1));
Y(2) = V*sin(y(end-1));
Y(3) = Vt*cos(y(end));
Y(4) = Vt*sin(y(end));



Y(end-1) = (1/V)*latax;
if t<2
    Y(end) = Vt/R;
else
    Y(end) = -2*Vt/R;
end
% Y(end) = sign(Switch_Matrix(3,i))*Vt/Switch_Matrix(1,i); % + randn(1);
% plot(y(1),y(2));
end