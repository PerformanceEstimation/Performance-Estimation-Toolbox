clear all; clc;
% In this example, we use a proximal point method for solving the 
% non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of a proximal method starting with an initial
% iterate satisfying ||x0-xs||<=1.

% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
F=P.DeclareFunction('Convex'); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=4;		% number of iterations
gam=1;		% step size

x=x0;
for i=1:N
    x=proximal_step(x,F,gam);
end
xN=x;

% (4) Set up the performance measure
fN=F.value(xN);
P.PerformanceMetric(fN-fs);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% The result should be (and is) 1/(4*N*gam)
% see Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization." to appear in SIAM Journal on Optimization
%     (2017)