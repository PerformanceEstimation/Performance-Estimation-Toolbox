clear all; clc;
% In this example, we use a fast gradient method for solving the L-smooth
% convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to verify the inequality (see [1] below)
%
%  lambda_{k}^2 (F(x_{k+1}) - F(xs)) + L/2 * ||z_{k+1}-xs||^2 
%                   <= lambda_{k-1}^2 (F(x_k) - F(xs)) + L/2 * ||z_k-xs||^2
% 
% when lambda_{k+1} = 1/2 * (1+sqrt(4*lambda_k^2+1)) and
%
%   y_{k}   = (1-tau_k) x_k + tau_k z_k
%   x_{k+1} = y_k - 1/L F'(y_k)
%   z_{k+1} = z_k - eta_k F'(y_k)
%
% with tau_k = 1/lambda_k and eta_k = (lambda_k^2 - lambda_{k-1}^2) / L
%
% [1] Nikhil Bansal, and Anupam Gupta.  "Potential-function proofsfor
%     first-order methods." (2019)

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.mu = 0;     % Strong convexity parameter
L        = 1;
param.L  = L;     % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point 
xk      = P.StartingPoint();		 % xk is some starting point
zk      = P.StartingPoint();		 % zk is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)


% (3) Algorithm
lambda   = 10;                          % lambda_k
h        = 1/param.L;                   % 1/L
lambdap1 = (1+sqrt(4*lambda^2+1))/2;    % lambda_{k+1}
tau      = 1/lambdap1;
eta      = (lambdap1^2 - lambda^2)/L;

yk            = (1-tau) * xk + tau * zk;
gyk           = F.gradient(yk);
xkp1          = yk - h * gyk;           % x_{k+1}
zkp1          = zk - eta * gyk;         % z_{k+1}

[  gk, fk ]   = F.oracle(xk);
[gkp1, fkp1 ] = F.oracle(xkp1);

% (4) Set up the performance measure:
% find the largest possible value for phi_{k+1}-phi_k (should be <= 0)


phi=@(lam,f,x)( lam*f+L/2*x^2);
objective = phi(lambdap1^2,fkp1-fs,(zkp1-xs))-phi(lambda^2,fk-fs,(zk-xs));
P.PerformanceMetric(objective);    


% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(objective)   % worst-case objective function accuracy

% The result should be <= 0 (see [1])
