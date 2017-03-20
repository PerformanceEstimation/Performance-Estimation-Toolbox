clear all; clc;
% In this example, we use a fixed-step gradient method for
% solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the gradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=P.AddObjective('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
h=1/param.L;		% step size
N=10;		% number of iterations

x=x0;
for i=1:N
    x=x-h*F.gradient(x);
    % % Alternate - shorter - form:
    % x=gradient_step(x,F,gam);
end
xN=x;

% (4) Set up the performance measure
fN=F.value(xN);                % g=grad F(x), f=F(x)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()
