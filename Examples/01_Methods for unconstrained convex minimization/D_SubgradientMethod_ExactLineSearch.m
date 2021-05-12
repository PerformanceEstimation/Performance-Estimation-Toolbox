clear all; clc;
% In this example, we use a subgradient method with exact line search for
% solving the non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
% where F(x) satisfies a Lipschitz condition; i.e., it has a bounded
% gradient ||g||<=R for all g being a subgradient of F at some point.
%
% The method originates from:
% [1] Y. Drori and A. Taylor. Efficient first-order methods for convex
% minimization: a constructive approach. Mathematical Programming (2019)
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the specific subgradient method starting
% with an initial iterate satisfying ||x0-xs||<=1.

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
x = cell(N+1,1);
g = cell(N+1,1);
f = cell(N+1,1);

x{1}        = x0;
[g{1},f{1}] = F.oracle(x{1});
d           = g{1};
for i=1:N
    y                        = i/(i+1) * x{i} + 1/(i+1) * x{1};
    [x{i+1}, g{i+1}, f{i+1}] = exactlinesearch_step(y,F,d);
    d                        = d + g{i+1};
end

P.PerformanceMetric(f{N+1}-fs);

% (5) Solve the PEP
out=P.solve()

% (6) Evaluate the output
double(f{N+1}-fs)

% The result should be 1/sqrt(N+1), see [1].
