% Define the parameters
alpha = [0.9, 0.92, 0.94];
eta_a = [2, 2, 2];
eta_b = [1, 1, 1];
rho = [0.25, 0.2, 0.15];
gamm = [0.1, 0.15, 0.2];
a = 35; b = 3; c = 28;
mu = 0.5;
eta = 0.5;
k = 1;

T = settlingTime(rho(1), eta);
xi_num = [-1, -2, -6].^T;

first_term = Caputo(alpha(1)-1)*sign(xi_num(1));
second_term = 2*Caputo(alpha(1)-1)*xi_num(1);
third_term = 1*Caputo(alpha(1)-2)*(xi_num(1)+abs(xi_num(1)).^mu*sign(xi_num(1)));

si = first_term + second_term + third_term; % Eq. (9)

f = @(x) sin(x);
first_term = (sign(si)/(sign(xi_num(1))+2));
second_term = ((1+2)*(abs(f(xi_num(1)))) + rho(1) + gamm(1));
third_term = (1*abs(Caputo(alpha(1)-1)*(xi_num(1)+abs(xi_num(1).^mu)*sign(xi_num(1)))));

ui_num = -((first_term * second_term) + third_term + k); % Eq. (21)

tMin = 0;
tMax = 200;
Nt = 100;

t = linspace(tMin, tMax, Nt);

% plot x1
D_1 = a*(xi_num(2)-xi_num(1))+((0.25*cos(2*t))-(0.1*sin(t))) + ui_num;

figure(1)
plot(t, D_1)
xlabel('t(s)')
ylabel('X1')


% %% plot x2
% first_term = Caputo(alpha(2)-1)*sign(xi_num(2));
% second_term = 2*Caputo(alpha(2)-1)*xi_num(2);
% third_term = 1*Caputo(alpha(2)-2)*(xi_num(2)+abs(xi_num(2)).^mu*sign(xi_num(2)));
% 
% si = first_term + second_term + third_term; % Eq. (9)
% 
% f = @(x) sin(x);
% first_term = (sign(si)/(sign(xi_num(2))+2));
% second_term = ((1+2)*(abs(f(xi_num(2)))) + rho(2) + gamm(2));
% third_term = (1*abs(Caputo(alpha(2)-1)*(xi_num(2)+abs(xi_num(2).^mu)*sign(xi_num(2)))));
% 
% ui_num = -((first_term * second_term) + third_term + k); % Eq. (21)
% 
% disp(ui_num);
% 
% D_2 = ((c - a)*xi_num(1)) + (c*xi_num(2)) - (xi_num(1) * xi_num(3)) + ((-0.2*cos(3*t)) + 0.15*sin(2*t)) + ui_num;
% 
% figure(2)
% plot(t, D_2)
% xlabel('t(s)')
% ylabel('X2')
% 
% %% plot x3
% first_term = Caputo(alpha(3)-1)*sign(xi_num(3));
% second_term = 2*Caputo(alpha(3)-1)*xi_num(3);
% third_term = 1*Caputo(alpha(3)-2)*(xi_num(3)+abs(xi_num(3)).^mu*sign(xi_num(3)));
% 
% si = first_term + second_term + third_term; % Eq. (9)
% 
% f = @(x) sin(x);
% first_term = (sign(si)/(sign(xi_num(3))+2));
% second_term = ((1+2)*(abs(f(xi_num(3)))) + rho(3) + gamm(3));
% third_term = (1*abs(Caputo(alpha(3)-1)*(xi_num(3)+abs(xi_num(3).^mu)*sign(xi_num(3)))));
% 
% ui_num = -((first_term * second_term) + third_term + k); % Eq. (21)
% 
% D_3 = (-b*xi_num(3)) + (xi_num(1) + xi_num(2)) + (0.15*cos(4*t)) - (0.2*sin(3*t)) + ui_num;
% 
% % plot results
% figure(3);
% plot(t, D_3)
% xlabel('t(s)')
% ylabel('X3')
% 
