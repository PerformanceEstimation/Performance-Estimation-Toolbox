clear all; clc;
% (0) Initialize an empty PEP
P=pet();

% (1) Set up the objective function
param.R=1;	% 'radius'-type constraint on the subgradient norms: ||g||<=1

% F is the objective function
F=P.AddObjective('SmoothConvexBoundedGradient',param); 

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();            % x0 is some starting point
[xs,fs]=F.OptimalPoint();        % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1);% Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm and (4) performance measure
N=3; % number of iterations
h=ones(N,1)*1/sqrt(N+1); % step sizes

x=x0;

% Note: the worst-case performance measure used in the PEP is the 
%       min_i (PerformanceMetric_i) (i.e., the best value among all
%       performance metrics added into the problem. Here, we use it
%       in order to find the worst-case value for min_i [F(x_i)-F(xs)]
for i=1:N
    [g,f]=F.oracle(x);
    P.PerformanceMetric(f-fs);
    x=x-h(i)*g;
end

[g,f]=F.oracle(x);
P.PerformanceMetric(f-fs);

% (5) Solve the PEP
P.solve()

