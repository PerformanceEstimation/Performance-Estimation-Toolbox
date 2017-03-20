clear all; clc;
% In this example, we use the optimized gradient method (OGM) for
% solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the gradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.
%
% Note that OGM is developped in the following two works:
%(1)Drori, Yoel, and Marc Teboulle. 
%   "Performance of first-order methods for smooth convex minimization:
%   a novel approach." Mathematical Programming 145.1-2 (2014): 451-482.
%
%(2)Kim, Donghwan, and Jeffrey A. Fessler. 
%   "Optimized first-order methods for smooth convex minimization."
%   Mathematical programming 159.1-2 (2016): 81-107.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter

F=P.AddObjective('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/param.L;		% step size
N=5;		% number of iterations

x=cell(N+1,1);%we store all the x's in a cell (for convenience)
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
fN=F.value(y);                % g=grad F(x), f=F(x)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% The result should be 1/2/theta^2
% see: Kim, Donghwan, and Jeffrey A. Fessler. 
%      "Optimized first-order methods for smooth convex minimization."
%      Mathematical programming 159.1-2 (2016): 81-107.