clear all; clc;
% In this example, we use a greedy first-order method (GFOM), or conjugate
% gradient, for solving the L-smooth (possibly mu-strongly) 
% convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the gradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.
%
% Exact worst-case results on GFOM can be found in:
% [1] Y. Drori and A. Taylor. Efficient first-order methods for convex
% minimization: a constructive approach. Mathematical Programming (2019)

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L     = 1;      % Smoothness parameter
param.mu    = 0.0;    % Strong convexity parameter

% F is the objective function
F=P.DeclareFunction('SmoothStronglyConvex',param); 

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10;		% number of iterations

x       = cell(N+1,1);% store iterates in a cell
g       = cell(N+1,1);% store gradients in a cell
x{1}    = x0;
g{1}    = F.gradient(x{1});
dirs{1} = g{1}; % search directions in a cell
for i = 1:N
    % span search in all directions contained in the cell dirs{}
    [x{i+1}, g{i+1}, f{i+1}] = exactlinesearch_step(x{i},F,dirs);
    dirs{2+(i-1)*2}  = x{i+1} - x{1}; % add search directions
    dirs{3+(i-1)*2}  = g{i+1};        % add search directions
end

% (4) Set up the performance measure
P.PerformanceMetric(f{N+1}-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(f{N+1}-fs)   % worst-case objective function accuracy
% The bounds should be 1/2/theta_N^2 (see [1]) with
theta    = cell(N+1,1);
theta{1} = 1;
for i = 1:N-1
    theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
end
theta{N+1}  = (1+sqrt(8*theta{N}^2+1))/2;
1/theta{N+1}^2/2