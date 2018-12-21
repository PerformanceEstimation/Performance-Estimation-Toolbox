clear all; clc;
% In this example, we use a greedy first-order method (GFOM), or conjugate
% gradient, for solving the L-smooth (possibly mu-strongly) 
% convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the gradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter
param.mu=0.0;   % Strong convexity parameter

% F is the objective function
F=P.DeclareFunction('SmoothStronglyConvex',param); 

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=5;		% number of iterations

x=cell(N+1,1);%we store all the x's in a cell (for convenience)
g=cell(N+1,1);%we store all the g's in a cell (for convenience)
x{1}=x0;
g{1}=F.gradient(x{1});
dirs{1}=g{1};
for i=1:N
    [x{i+1}, g{i+1}] = exactlinesearch_step(x{i},F,dirs);
    dirs{2+(i-1)*2}  = x{i+1} - x{1};
    dirs{3+(i-1)*2}  = g{i+1};
end

% (4) Set up the performance measure
fN=F.value(x{N+1});                % g=grad F(x), f=F(x)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy
% The results are the same as those for the optimized gradient method.