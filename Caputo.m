% clc; clear; close all;

function L1_Caputo_Der = Caputo(alpha)
% Inputs

h = 0.1; % Stepsize

x(1) = 0; % initial value

xlast = 1; % final value

x = x(1):h:xlast; % interval

n = ceil((xlast - x(1)) / h); % No of steps

% alpha = 0.5; % fractional order

f = @(x) sin(x); % the function to be differentiated

% Algorithm
k = 1:n - 1;

A = (h.^(-alpha))./(gamma(2-alpha));

Bnk = ((n-k).^(1-alpha)-(n-k-1).^(1-alpha));

L1_Caputo_Der = A * sum(Bnk.*(f(x(k+1))-f(x(k))));

% Exact Solution

% x=xlast;
% 
% Exact=(cos(x)*fresnelc(sqrt(x)*sqrt(2)/sqrt(pi)) + sin(x)*fresnels(sqrt(x)*sqrt(2)/sqrt(pi))*sqrt(2));
% 
% Error = abs(Exact-L1_Caputo_Der); % absolute error
% 
% disp('      Steps       Stepsize        Exact       Approximate     Error');
% disp('---------------------------------------------------------------------');
% Results = [n  h   Exact       L1_Caputo_Der       Error]

end
