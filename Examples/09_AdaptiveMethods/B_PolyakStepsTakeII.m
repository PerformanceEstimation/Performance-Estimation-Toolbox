clear all; clc;
% In this example, we use a variant of Polyak steps for solving the
% L-smooth m-strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of f(x_{k+1})-f_* when 
% x_{k+1} is obtained by doing 1 steps of the method starting from x_k
% satisfying f(x_{k})-f_*==1.
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

% (3) Algorithm
gamma = 2/(L+m);              % step size parameter
% (gamma must be contained in the interval [1/L,(2-m/L)/L]; worst-case
% achieved for gamma = 2/(L+m). Outside this interval, worst-case is zero
% (we are at an optimal point))

[g0,f0] = F.oracle(x0);
x1      = x0 - gamma * g0;
[g1,f1] = F.oracle(x1);

P.InitialCondition( f0-fs <= 1); % Initial condition
P.AddConstraint(  g0^2 == 2*L*(2-gamma)*(f0 - fs));

% (4) Set up the performance measure
obj = f1-fs ;
P.PerformanceMetric(obj);     

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output; the result should be (see [1])
% (gamma * L - 1) * (L*gamma*(3-gamma*(L+m))-1)
[double(obj) (gamma * L - 1) * (L*gamma*(3-gamma*(L+m))-1)]