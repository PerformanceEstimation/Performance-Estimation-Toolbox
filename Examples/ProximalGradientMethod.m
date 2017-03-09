clear all; clc;

% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
paramf1.mu=.1;	% Strong convexity parameter
paramf1.L=1;      % Smoothness parameter
f1=my_pep.AddComponentObjective('SmoothStronglyConvex',paramf1);
f2=my_pep.AddComponentObjective('Convex');
F=f1+f2; % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
[g0,f10]=f1.oracle(x0);
[s0,f20]=f2.oracle(x0);
my_pep.AddInitialCondition((g0+s0)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/paramf1.L;		% step size
N=1;		% number of iterations

x=cell(N+1,1);
x{1}=x0;
for i=1:N
    xint=gradient_step(x{i},f1,gam);
    x{i+1}=proximal_step(xint,f2,gam);
    s=(xint-x{i+1})/gam;
end

% (4) Set up the performance measure
[gN,f1N]=f1.oracle(x{N+1});
[sN,f2N]=f1.oracle(x{N+1});
my_pep.AddPerformanceConstraint((gN+s)^2);

% (5) Solve the PEP
my_pep.solve()




