clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.M1=Inf;	% 
param.M2=1;	% 
param.L=Inf;      % Smoothness parameter

F=my_pep.AddComponentObjective('SmoothConvexBoundedGradient',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)

my_pep.AddConstraint((x0-xs)^2<=1);

% (3) Algorithm
N=3;		% number of iterations
gam=1/sqrt(N+1);		% step size

x=cell(N+1,1);
x{1}=x0;
for i=1:N
    [g,f]=F.oracle(x{i});
    my_pep.AddPerformanceConstraint(f-fs);
    x{i+1}=x{i}-1/sqrt(i+1)*g;
end

[g,f]=F.oracle(x{N+1});
my_pep.AddPerformanceConstraint(f-fs);

% (5) Solve the PEP
my_pep.solve()

