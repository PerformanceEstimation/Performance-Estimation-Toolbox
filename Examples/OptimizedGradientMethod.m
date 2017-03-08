clear all; clc;
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=my_pep.AddComponentObjective('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
my_pep.AddInitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/param.L;		% step size
N=5;		% number of iterations

x=cell(N+1,1);
x{1}=x0;
y=x0;
theta=1;
for i=1:N
    x{i+1}=gradient_step(y,F,gam);
    theta_prev=theta;
    if i<N
        theta=(1+sqrt(4*theta^2+1))/2;
    else
        theta=(1+sqrt(8*theta^2+1))/2;
    end
    y=x{i+1}+(theta_prev-1)/theta*(x{i+1}-x{i})+theta_prev/theta*(x{i+1}-y);
end

% (4) Set up the performance measure
[g,f]=F.oracle(y);                % g=grad F(x), f=F(x)
my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
my_pep.solve()