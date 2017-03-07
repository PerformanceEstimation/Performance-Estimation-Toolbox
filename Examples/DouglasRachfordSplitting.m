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
w0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
[as,~]=f1.oracle(xs,1);
bs=-as;
lambda=0.9; ws=xs+lambda*bs;

my_pep.AddInitialCondition(w0-ws,1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=5;            % number of iterations
gam=lambda;		% step size

w=cell(N+1,1);
w{1}=w0;

for i=1:N
    x=proximal_step(w{i},f2,gam);
    y=proximal_step(2*x-w{i},f1,gam);
    w{i+1}=y-x+w{i};
end

% (5) Set up the performance measure
my_pep.AddPerformanceConstraint(w{i+1}-ws);

% (5) Solve the PEP
my_pep.solve()

