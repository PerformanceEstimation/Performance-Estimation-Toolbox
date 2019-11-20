clear all; clc;
% In this example, we use a fast gradient method for solving the L-smooth 
% mu-strongly convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L   = 1;                % Smoothness parameter
param.mu  = .1;               % Strong convexity parameter
kappa     = param.L/param.mu; % condition number

% F is the objective function
F = P.DeclareFunction('SmoothStronglyConvex',param); 

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();         % x0 is some starting point
[xs,fs] = F.OptimalPoint();          % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 5;		% number of iterations

x    = cell(N+1,1);% store iterates in a cell
y    = cell(N+1,1);% store iterates in a cell
x{1} = x0;
y{1} = x0;

beta = (sqrt(kappa)-1)/(sqrt(kappa)+1); % momentum
for i=1:N
    x{i+1} = gradient_step(y{i}, F, 1/param.L);
    y{i+1} = (1+beta) * x{i+1} - beta * x{i};
end

% (4) Set up the performance measure
fN = F.value(x{N+1});         % fN=F(xN)
P.PerformanceMetric(fN-fs);   % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% Should be better than the standard guarantee
% f(xN)-f(xs)<= param.L*(1-1/sqrt(kappa))^N*double((x0-xs)^2)
