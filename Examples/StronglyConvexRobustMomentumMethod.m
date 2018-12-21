clear all; clc;
% In this example, we use the robust momentum method for solving the 
% L-smooth mu-strongly convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x).
% We show how to compute the rate of the Lyapunov function developped in
%(**) Cyrus, S., Hu, B., Van Scoy, B., & Lessard, L. "A robust accelerated 
%     optimization algorithm for strongly convex functions." In 2018 Annual
%     American Control Conference (ACC) (pp. 1376-1381). IEEE.

% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
L  = 1;  param.L=L;      % Smoothness parameter
mu = .2; param.mu=mu;    % Strong convexity parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
xm1=P.StartingPoint();            % xm1 is some starting point
x0=P.StartingPoint();             % x0 is some starting point
[xs,fs]=F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)

% (3) Algorithm' parameters
kappa = param.L/param.mu;
lam   = .1;
rho   = lam*(1-1/kappa) + (1-lam)*(1-1/sqrt(kappa));

alpha  = kappa*(1-rho)^2*(1+rho)/L;
beta   = kappa*rho^3/(kappa-1);
gamma  = rho^3/((kappa-1)*(1-rho)^2*(1+rho));
lambda = mu^2*(kappa-kappa*rho^2-1)/(2*rho*(1-rho));
nnu    = (1+rho)*(1-kappa+2*kappa*rho-kappa*rho^2)/(2*rho);

% (4) Perform two iterations according to the notations in (**)
y0      = x0+gamma*(x0-xm1);
[g0,F0] = F.oracle(y0);
x1      = x0 + beta*(x0-xm1) - alpha*g0;
y1      = x1 + gamma*(x1-x0);
[g1,F1] = F.oracle(y1);
x2      = x1+beta*(x1-x0)-alpha*g1;

z0      = (x0-rho^2*xm1)/(1-rho^2);
z1      = (x1-rho^2*x0)/(1-rho^2);
z2      = (x2-rho^2*x1)/(1-rho^2);

% (5) Evaluate the Lyapunov function at the first and second iterations
q0 = (L-mu)*(F0-fs-mu/2*(y0-xs)^2)-1/2*(g0-mu*(y0-xs))^2;
q1 = (L-mu)*(F1-fs-mu/2*(y1-xs)^2)-1/2*(g1-mu*(y1-xs))^2;

initLyapunovValue   = lambda*(z1-xs)^2+q0;
finalLyapunovValue  = lambda*(z2-xs)^2+q1;

% (6) Normalize the initial value of the Lyapunov function to one and observe
% that the value of the second indeed corresponds to the rate.
P.InitialCondition(initLyapunovValue <= 1);     % (initial condition)
P.PerformanceMetric(finalLyapunovValue);   % (performance measure)

% (5) Solve the PEP
P.solve(1);
% the value of final should match rho^2
[double(finalLyapunovValue) rho^2]