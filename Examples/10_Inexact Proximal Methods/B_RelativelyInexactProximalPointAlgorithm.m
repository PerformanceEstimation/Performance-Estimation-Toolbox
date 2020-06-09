clear all; clc;
% In this example, we use a relatively inexact proximal point algorithm
% for solving the non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.
%
% The inexact proximal point algorithm produces a x_{k+1} such that 
% gamma* (h(x_{k+1}) + h^*(v) - <v;x_{k+1}>) <= sigma^2/2*||x_{k+1}-x_k||^2
% with v = (x_k-x_{k+1})/gamma.
%
% The example is used in [1, Section 3].

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

% (3) Set up the method
N               = 8;
lambda          = 10;
sigma           =  sqrt(.4);
opt.criterion   = 'PD_gapIII';

x           = cell(N+1,1);
fx          = cell(N+1,1);
theta       = zeros(1,N);
x{1}        = x0;

for i = 1:N
    [x{i+1},~,fx{i+1},~,~,~,epsVar] = inexact_proximal_step(x{i},f,lambda,opt);
    P.AddConstraint(epsVar <= sigma^2/2 * (x{i}-x{i+1})^2);
end


% (4) Set up the objective
P.PerformanceMetric(fx{N+1}-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output
double(fx{N+1}-fs)
