function ST = settlingTime(p, eta)
% Define the system of differential equations
f = @(t, x) [-x(1)^3 + x(2); -x(1) + x(2)^3];

% Define the initial condition
x0 = [1; 1];

% Define the time interval to solve the ODE over
tspan = [0 200];

% Define the function V and its derivative
syms x1 x2
V = x1^2 + x2^2;
dVdx = jacobian(V, [x1, x2]) * f(tspan, [x1; x2]);


% Define the constants p and eta
p = p;
eta = eta;

% Define the settling time function T
T = @(x0) (subs(V^(1-eta) / (p*(1-eta)), {x1, x2}, {x0(1), x0(2)}));

% Solve the ODE using ode45
[t, x] = ode45(f, tspan, x0);

% Calculate the settling time for the initial condition
settling_time = T(x0);

ST = settling_time;
end