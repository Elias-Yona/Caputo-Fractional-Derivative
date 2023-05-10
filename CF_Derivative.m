alpha = [0.9, 0.92, 0.94];
a = 35; b = 3; c = 28;
rho = [0.25, 0.2, 0.15];
gamma = [0.1, 0.15, 0.2];
eta_a = [2, 2, 2];
eta_b = [1, 1, 1];
upsilon = 0.5;
T = 1.^3;
x_0 = [-1, -2, -6].^T; % raised to T_3
disp(x_0);
t = 10;
u = 10;


tMin = 0;
tMax = 100;
Nt = 200;

t = linspace(tMin, tMax, Nt);

% Equation (31)
D_1 = a * (x_0(2) - x_0(1)) + ((0.25 * cos(2*t)) - (0.1 * sin(t))); % + (u * t);

% plot results
figure(1)
plot(t, D_1)
xlabel('t(s)')
ylabel('X1')

D_2 = ((c - a)*x_0(1)) + (c*x_0(2)) - (x_0(1) * x_0(3)) + ((-0.2*cos(3*t)) + 0.15*sin(2*t)); % + (u * t);

% plot results
figure(2)
plot(t, D_2)
xlabel('t(s)')
ylabel('X2')

D_3 = (-b*x_0(3)) + (x_0(1) + x_0(2)) + (0.15*cos(4*t)) - (0.2*sin(3*t)); % + (u * t);

% plot results
figure(3);
plot(t, D_3)
xlabel('t(s)')
ylabel('X3')
% axis([tMin, tMax, 0, 15]))
