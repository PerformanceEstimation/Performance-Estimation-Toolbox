clear all; clc;
% In this example, we use a fast proximal gradient method
% for solving the composite convex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) } 
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is L-smooth and convex and where f_2(x) is a convex
% indicator function.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L = 1;      % Smoothness parameter

f1 = P.DeclareFunction('SmoothStronglyConvex',param); 
f2 = P.DeclareFunction('ConvexIndicator'); 
F  = f1 + f2;

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();            % x0 is some starting point
[xs,fs] = F.OptimalPoint();             % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1);    % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gamma = 1/param.L;  % stepsize
N     = 20;         % number of iterations

x    = cell(N+1,1); % we store the iterate in a cell for convenience
x{1} = x0;
y    = x0;
for i = 1:N
    xint    = gradient_step(y,f1, gamma);
    x{i+1}  = projection_step(xint, f2);
    y       = x{i+1} + (i-1)/(i+2) * (x{i+1}-x{i});
end

% (4) Set up the performance measure
fN  = F.value(x{N+1});
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% Result should be 2/(N^2+5*N+2)
% see Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization." SIAM Journal on Optimization (2017)

