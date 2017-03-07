clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

f1=my_pep.AddComponentObjective('SmoothStronglyConvex',param); 
f2=my_pep.AddComponentObjective('Convex'); 
F=f1+f2;

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
my_pep.AddInitialCondition(x0-xs,1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/param.L;		% step size
N=3;		% number of iterations

x=cell(N+1,1);
x{1}=x0;
y=x0;
for i=1:N
    xint=gradient_step(y,f1,gam);
    x{i+1}=proximal_step(xint,f2,gam);
    y=x{i+1}+(i-1)/(i+2)*(x{i+1}-x{i});
end

% (4) Set up the performance measure
[g,f]=F.oracle(x{N+1});                % g=grad F(x), f=F(x)
my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
my_pep.solve()

fprintf('Output should be LR^2/2 * 4/(N^2+5*N+2)=%5.3e\n',2/(N^2+5*N+2));


