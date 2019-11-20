clear all; clc;
% In this example, we use the optimized gradient method (OGM) with 
% exact line-search for solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of OGM-LS starting with an initial iterate
% satisfying ||x0-xs||<=1.
%
% The method originates from:
% [1] Y. Drori and A. Taylor. Efficient first-order methods for convex
% minimization: a constructive approach. Mathematical Programming (2019)


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L = 1;      % Smoothness parameter

% F is the objective function:
F = P.DeclareFunction('SmoothStronglyConvex',param);

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10;		% number of iterations

x            = cell(N+1,1);% store iterates in a cell
g            = cell(N+1,1);% store iterates in a cell
f            = cell(N+1,1);% store iterates in a cell
x{1}         = x0;
[g{1}, f{1}] = F.oracle(x{1});
theta        = cell(N+1,1);
theta{1}     = 1;

agglomerate  = theta{1} * g{1};
for i = 1:N
    if i<N
        theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
    else
        theta{i+1}  = (1+sqrt(8*theta{i}^2+1))/2;
    end
    y      = (1 - 1/theta{i+1}) * x{i} + 1/theta{i+1} * x{1};
    d      = (1 - 1/theta{i+1}) * g{i} + 2/theta{i+1} * agglomerate;
    
    [x{i+1}, g{i+1}, f{i+1}] = exactlinesearch_step(y,F,d);
    agglomerate              = agglomerate + theta{i+1} * g{i+1};
end
          
P.PerformanceMetric(f{N+1}-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(f{N+1}-fs)   % worst-case objective function accuracy
% The result should be 1/2/theta{N+1}^2 (see [1])
