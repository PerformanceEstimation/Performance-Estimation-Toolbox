clear all; clc;
% In this example, we use a subgradient method for
% solving the non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
% where F(x) satisfies a Lipschitz condition; i.e., it has a bounded
% gradient ||g||<=R for all g being a subgradient of F at some point.
%
% We show how to compute the worst-case value of min_i F(xi)-F(xs) when xi is
% obtained by doing i steps of a subgradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.R = 1;	% 'radius'-type constraint on the subgradient norms: ||g||<=1

% F is the objective function
F = P.DeclareFunction('ConvexBoundedGradient',param); 

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();         % x0 is some starting point
[xs,fs] = F.OptimalPoint();          % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm and (4) performance measure
N = 5;                     % number of iterations
h = @(k)(1/sqrt(N+1)); % step sizes

x = x0;

% Note: the worst-case performance measure used in the PEP is the 
%       min_i (PerformanceMetric_i) (i.e., the best value among all
%       performance metrics added into the problem. Here, we use it
%       in order to find the worst-case value for min_i [F(x_i)-F(xs)]

% we create an array to save all function values (so that we can evaluate
% them afterwards)
f_saved = cell(N+1,1);
for i=1:N
    [g,f]       = F.oracle(x);
    f_saved{i}  = f;
    P.PerformanceMetric(f-fs,'min');
    x = x - h(i) * g;
end

[g,f]        = F.oracle(x);
f_saved{N+1} = f;
P.PerformanceMetric(f-fs,'min');

% (5) Solve the PEP
P.solve(3)

% (6) Evaluate the output
for i=1:N+1
    f_saved{i} = double(f_saved{i});
end
f_saved
% The result should be 1/sqrt(N+1).
