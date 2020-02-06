clear all; clc;
% In this example, we use a variant of Polyak steps for solving the
% L-smooth m-strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of ||x_{k+1}-x_*||^2 when 
% x_{k+1} is obtained by doing 1 steps of the method starting from x_k
% satisfying ||xk-xs||==1.
%
% Result to be compared with that of
% [1] Mathieu Barre, Adrien Taylor, Alexandre d'Aspremont (2020). 
%     "Complexity Guarantees for Polyak Steps with Momentum." 

% (0) Initialize an empty PEP
P = pep();

L = 1; m = .1;
% (1) Set up the objective function
param.mu = m;     % Strong convexity parameter
param.L  = L;     % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gamma = 2/(L+m);              % step size parameter
% (gamma must be contained in the interval [1/L,1/m]; worst-case achieved
% for gamma = 2/(L+m). Outside this interval, worst-case is zero
% (we are at an optimal point))

[g0,f0] = F.oracle(x0);
x1      = x0 - gamma * g0;

P.AddConstraint( gamma * g0^2 == 2*(f0 - fs));

% (4) Set up the performance measure
obj = (x1-xs)^2;
P.PerformanceMetric(obj);     

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output; the result should be (see [1])
% (gamma * L - 1) * (1 - gamma * m)/(gamma*(L+m)-1)
[double(obj) (gamma * L - 1) * (1 - gamma * m)/(gamma*(L+m)-1)]