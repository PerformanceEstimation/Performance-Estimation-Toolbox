clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.mu=.1;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=my_pep.AddComponentObjective('SmoothStronglyConvex',param); 
% F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
[g0, f0]=F.oracle(x0);               
my_pep.AddInitialCondition(f0-fs-1<=0); % Add an initial condition f0-fs<= 1

% (3) Algorithm
N=1;
eps=0.1;
x=x0;
for i=1:N
    d=inexactsubgradient(x,F,eps,0);
    x=exactlinesearch_step(x,F,d);
    [g,f]=F.oracle(x);
end

% (4) Set up the performance measure
my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as ||g||^2

% (5) Solve the PEP
my_pep.solve()

