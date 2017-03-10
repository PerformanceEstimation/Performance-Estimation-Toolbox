clear all; clc;

% (0) Initialize an empty PEP
P=pet();

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
    [g,f]=F.oracle(x);  % g=grad F(x), f=F(x)
    x=x-h*g;
    % % Alternative - shorter - form:
    % x=gradient_step(x,F,gam);
end

% (4) Set up the performance measure
[g,f]=F.oracle(x);                % g=grad F(x), f=F(x)
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()
