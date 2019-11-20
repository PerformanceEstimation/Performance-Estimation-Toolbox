clear all; clc;
% In this example, we use the heavy ball method for solving the L-smooth 
% convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x);
% where F(x) is L-smooth and convex.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L = 1;      % Smoothness parameter
F = P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();             % x0 is some starting point
[xs,fs] = F.OptimalPoint();              % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1);     % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10; % number of iterations

% parameters of the heavy-ball: 0 <= beta < 1 and 0 < alpha < 2*(1+beta)
alpha = 1;      
beta  = 1/2;

x = cell(N+1,1);% we store iterates in a cell 
g = cell(N+1,1);% we store gradients in a cell

x{1}         = x0;
[g{1}, f{1}] = F.oracle(x{1});
x{2}         = x{1} - alpha/param.L * g{1};
[g{2}, f{2}] = F.oracle(x{2});
for i = 2:N
    x{i+1}           = x{i} - alpha/param.L * g{i} + beta * (x{i} - x{i-1});
    [g{i+1}, f{i+1}] = F.oracle(x{i+1});
end

% (4) Set up the performance measure
P.PerformanceMetric(f{N+1}-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(f{N+1}-fs)   % worst-case objective function accuracy
