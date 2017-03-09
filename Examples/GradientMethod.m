clear all; clc;

% (0) Initialize an empty PEP
P=pet();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=P.AddComponentObjective('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.AddInitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/param.L;		% step size
N=10;		% number of iterations
[g0,f0]=F.oracle(x0,'xz');
g=g0;
x=x0;
for i=1:N
% %     x=gradient_step(x{i},F,gam);%x=x-gam/L*grad(x)
%     % This is the short form for:
    		% g=grad F(x), f=F(x)
    x=x-gam/param.L*g;
    [g,f]=F.oracle(x);
end

% (4) Set up the performance measure
% [g,f]=F.oracle(x);                % g=grad F(x), f=F(x)
P.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()
