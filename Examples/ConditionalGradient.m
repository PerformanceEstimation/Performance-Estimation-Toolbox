clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter
paramf2.M1=1;
paramf2.M2=Inf;
f1=my_pep.AddComponentObjective('SmoothStronglyConvex',param); 
f2=my_pep.AddComponentObjective('ConvexIndicator',paramf2);
F=f1+f2;

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
[g0, f0]=f1.oracle(x0);
[s0, h0]=f2.oracle(x0);

x=x0;
N=2;
for i=1:N
    [dir,~]=f1.oracle(x);
    y=linearoptimization_step(dir,f2);
    lambda=2/(1+i);
    x=(1-lambda)*x+lambda*y;
end

% (4) Set up the performance measure
[g,f]=F.oracle(x);                % g=grad F(x), f=F(x)
my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
my_pep.solve()

