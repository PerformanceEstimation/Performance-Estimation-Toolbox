clear all; clc;
% In this example, we use the optimized gradient method (OGM) for
% solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of OGM starting with an initial iterate
% satisfying ||x0-xs||<=1.
%
% Note that OGM is developped in the following two works:
%[1] Drori, Yoel, and Marc Teboulle.
%    "Performance of first-order methods for smooth convex minimization:
%    a novel approach." Mathematical Programming 145.1-2 (2014): 451-482.
%
%[2] Kim, Donghwan, and Jeffrey A. Fessler.
%    "Optimized first-order methods for smooth convex minimization."
%    Mathematical programming 159.1-2 (2016): 81-107.


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L = 1;      % Smoothness parameter

% F is the objective function:
F = P.DeclareFunction('SmoothStronglyConvex',param);

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10;		% number of iterations

x        = cell(N+1,1);% store iterates in a cell
x{1}     = x0;
y        = x0;
theta    = cell(N+1,1);
theta{1} = 1;
for i = 1:N
    x{i+1} = gradient_step(y,F,1/param.L);
    if i<N
        theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
    else
        theta{i+1}  = (1+sqrt(8*theta{i}^2+1))/2;
    end
    y = x{i+1}+(theta{i}-1)/theta{i+1}*(x{i+1}-x{i})+...
        theta{i}/theta{i+1}*(x{i+1}-y);
end

% (4) Set up the performance measure
fN = F.value(y);              
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% The result should be 1/2/theta{N+1}^2
% see: Kim, Donghwan, and Jeffrey A. Fessler.
%      "Optimized first-order methods for smooth convex minimization."
%      Mathematical programming 159.1-2 (2016): 81-107.