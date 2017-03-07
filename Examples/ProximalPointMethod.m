clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
F=my_pep.AddComponentObjective('Convex'); % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
my_pep.AddInitialCondition(x0-xs,1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=4;		% number of iterations
gam=1;		% step size

x=cell(N+1,1);
x{1}=x0;
for i=1:N
    x{i+1}=proximal_step(x{i},F,gam);
end

% (5) Set up the performance measure
[g,f]=F.oracle(x{N+1});
my_pep.AddPerformanceConstraint(f-fs);

% (5) Solve the PEP
my_pep.solve()

