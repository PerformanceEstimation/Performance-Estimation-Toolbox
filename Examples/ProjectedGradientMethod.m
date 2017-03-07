clear all; clc;

% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
paramf1.mu=.1;	% Strong convexity parameter
paramf1.L=1;      % Smoothness parameter
paramf2.M1=Inf; % No bounded domain restriction
paramf2.M2=Inf; % No bounded domain restriction

f1=my_pep.AddComponentObjective('SmoothStronglyConvex',paramf1);
f2=my_pep.AddComponentObjective('ConvexIndicator',paramf2); 
%f2 is the indicator function for the set Q
F=f1+f2; % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
my_pep.AddInitialCondition(x0-xs,1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/paramf1.L;		% step size
N=1;		% number of iterations

x=cell(N+1,1);
x{1}=x0;
for i=1:N
    xint=gradient_step(x{i},f1,gam);
    x{i+1}=proximal_step(xint,f2,gam);
end

% (4) Set up the performance measure
my_pep.AddPerformanceConstraint(x{N+1}-xs);

% (5) Solve the PEP
my_pep.solve()


