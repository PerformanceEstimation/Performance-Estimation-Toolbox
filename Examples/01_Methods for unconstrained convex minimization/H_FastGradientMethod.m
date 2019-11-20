clear all; clc;
% In this example, we use a fast gradient method for solving the L-smooth 
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
F = P.DeclareFunction('SmoothStronglyConvex',param); % objective function

% (2) Set up the starting point and initial condition
x0       = P.StartingPoint();        % x0 is some starting point
[xs,fs ] = F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2 <= 1);  % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10;  % number of iterations

x        = cell(N+1,1);% store iterates in a cell
x{1}     = x0;
y        = x0;
theta    = cell(N+1,1);
theta{1} = 1;
for i = 1:N
    x{i+1}      = gradient_step(y, F, 1/param.L);
    theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
    y           = x{i+1} + (theta{i}-1)/theta{i+1} * (x{i+1}-x{i});
end

% (4) Set up the performance measure
fN = F.value(x{N+1});         % fN=F(xN)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% Should be better than the standard guarantee from FISTA: 2/(N+1)^2
%
% See Beck, Amir, and Marc Teboulle. 
%     "A fast iterative shrinkage-thresholding algorithm 
%     for linear inverse problems." 
%     SIAM journal on imaging sciences 2.1 (2009): 183-202.
