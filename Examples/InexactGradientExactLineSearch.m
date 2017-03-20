clear all; clc;
% In this example, we use an inexact gradient method with exact line search
% for solving the L-smooth mu-strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
%   Starting from an iterate x0, the method performs at each iteration
%   an exact line search step in a direction d satisfying a relative 
%   accuracy criterion
%   (gi is the gradient of F at xi)
%       ||d-gi||<=eps*||gi|| (**)
%   that is, the method evaluates
%       gamma=argmin_gamma F(xi-gamma*d) for d satisfying (**)
%   and performs the update
%       x{i+1}=xi-gamma*d.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying F(x0)-F(xs)<=1.
%
% The full approach (based on convex relaxations) is available in 
%  De Klerk, Etienne, François Glineur, and Adrien B. Taylor. 
%  "On the worst-case complexity of the gradient method with exact 
%  line search for smooth strongly convex functions." 
%  Optimization Letters (2016).


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.mu=.1;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=P.AddObjective('SmoothStronglyConvex',param); 
% F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 	 % xs is an optimal point, and fs=F(xs)
[g0, f0]=F.oracle(x0);               
P.InitialCondition(f0-fs<=1); % Add an initial condition f0-fs<= 1

% (3) Algorithm
N=2;
eps=0.1;
x=x0;
for i=1:N
    d=inexactsubgradient(x,F,eps,0);
    x=exactlinesearch_step(x,F,d);
end
fN=F.value(x);
% (4) Set up the performance measure
P.PerformanceMetric(fN-fs); % Worst-case evaluated as ||g||^2

% (5) Solve the PEP
P.solve()

% The result should be
%((param.L*(1+eps)-param.mu*(1-eps))/(param.L*(1+eps)+param.mu*(1-eps)))^(2*N)


