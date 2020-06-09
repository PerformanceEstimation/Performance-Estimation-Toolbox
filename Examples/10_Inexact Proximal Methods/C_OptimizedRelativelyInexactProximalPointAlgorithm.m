clear all; clc;
% In this example, we use the optimized relatively inexact proximal point 
% algorithm (ORIP) for solving the non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of ORIP starting with an initial
% iterate satisfying ||x0-xs||<=1.
%
% The method originates from:
% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
f = P.DeclareFunction('Convex');

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=f.OptimalPoint(); 	 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((xs-x0)^2<=1)

% (3) Set up the method0
N               = 10;
lambda          = 2;
sigma           = 3;
opt.criterion   = 'Orip-style';

x           = cell(N+1,1);
fx          = cell(N+1,1);
v           = cell(N,1);
y           = cell(N,1);
z           = cell(N+1,1);
theta       = zeros(1,N);
x{1}        = x0;
z{1}        = x0;

for i = 1:N
    theta(i+1)              = (1+sqrt(4*theta(i)^2+1))/2;
    y{i}                    = (1-1/theta(i+1)) * x{i} + 1/theta(i+1) * z{i};
    [x{i+1},~,fx{i+1},~,v{i},~,epsVar] = inexact_proximal_step(y{i},f,lambda,opt);
    z{i+1}                  = z{i} - 2*lambda/(1+sigma) * theta(i+1) * v{i};
    P.AddConstraint( epsVar <= sigma/(1+sigma) * v{i}^2);
end


% (4) Set up the objective
P.PerformanceMetric(fx{N+1}-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output and compare it with [1, Theorem 4.2]
[double(fx{N+1}-fs) (1+sigma)/4/lambda/theta(N+1)^2]

