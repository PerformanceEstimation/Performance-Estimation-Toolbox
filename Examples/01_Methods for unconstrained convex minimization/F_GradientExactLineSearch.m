clear all; clc;
% In this example, we use a gradient method with exact line search
% for solving the L-smooth mu-strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
%   Starting from an iterate x0, the method performs at each iteration
%   an exact line search step in the steepest descent direction
%       gamma=argmin_gamma F(xi-gamma*gi), with gi the gradient of F at xi,
%   and performs the update
%       x{i+1}=xi-gamma*gi.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying F(x0)-F(xs)<=1.
%
% The detailed approach (based on convex relaxations) is available in 
%  De Klerk, Etienne, François Glineur, and Adrien B. Taylor. 
%  "On the worst-case complexity of the gradient method with exact 
%  line search for smooth strongly convex functions." 
%  Optimization Letters (2017).

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.mu = .1;	% Strong convexity parameter
param.L  = 1;      % Smoothness parameter

F = P.DeclareFunction('SmoothStronglyConvex',param); 
% F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
[g0,f0] = F.oracle(x0);              
P.InitialCondition(f0-fs<=1);        % Initial condition f0-fs<= 1

% (3) Algorithm
N = 2;
x = x0;
for i = 1:N
    g = F.gradient(x);
    
    % exact line search on F from point x and in direction g:
    x = exactlinesearch_step(x,F,g);
end

% (4) Set up the performance measure
[g,f] = F.oracle(x);
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(f-fs)   % worst-case objective function accuracy

% The result should be
%((param.L-param.mu)/(param.L+param.mu))^(2*N)


