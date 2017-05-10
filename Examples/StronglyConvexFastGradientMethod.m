clear all; clc;
% In this example, we use a fast gradient method for solving the L-smooth 
% mu-strongly convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter
param.mu=.1;    % Strong convexity parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();             % x0 is some starting point
[xs,fs]=F.OptimalPoint();              % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=10;		% number of iterations

x=cell(N+1,1);% we store the iterates in a cell for convenience
x{1}=x0;
y=x0;

kappa=param.mu/param.L;
coef=(1-sqrt(kappa))/(1+sqrt(kappa));
for i=1:N
    x{i+1}=gradient_step(y,F,1/param.L);
    y=(1+coef)*x{i+1}-coef*x{i};
end

% (4) Set up the performance measure
fN=F.value(x{N+1});         % fN=F(xN)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% Should be better than the standard guarantee
% f(xN)-f(xs)<= param.L*(1-sqrt(kappa))^N*double((x0-xs)^2)
