clear all; clc;
% In this example, we use a fast gradient method for solving the L-smooth 
% convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x);
% where F(x) is L-smooth and convex.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pet();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter

F=P.AddObjective('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();             % x0 is some starting point
[xs,fs]=F.OptimalPoint();              % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=15;		% number of iterations

x=cell(N+1,1);% we store the iterates in a cell for convenience
x{1}=x0;
y=x0;
theta=1;
for i=1:N
    x{i+1}=gradient_step(y,F,1/param.L);
    theta_prev=theta;
    theta=(1+sqrt(4*theta^2+1))/2;
    y=x{i+1}+(theta_prev-1)/theta*(x{i+1}-x{i});
end

% (4) Set up the performance measure
[g,f]=F.oracle(x{N+1});                % g=grad F(x), f=F(x)
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% Should be better than the standard guarantee from FISTA:
% 2/(N+1)^2
%
% See Beck, Amir, and Marc Teboulle. 
%     "A fast iterative shrinkage-thresholding algorithm 
%     for linear inverse problems." 
%     SIAM journal on imaging sciences 2.1 (2009): 183-202.