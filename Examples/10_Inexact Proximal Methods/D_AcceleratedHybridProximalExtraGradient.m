clear all; clc;
% In this example, we use an accelerated hybrid proximal extragradient 
% for solving the non-smooth (strongly) convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
%
% We show how to verify the potential function of [2, Theorem 4.3]: 
% A_{k+1} (F(x_{k+1})-F_*) + (1+mu*A_{k+1})/2 * ||z_{k+1}-x_*||^2 
%             <=  A_{k} (F(x_{k})-F_*) + (1+mu*A_{k})/2 * ||z_{k}-x_*||^2
% By denoting 
%   Phi_k = A_{k} (F(x_{k})-F_*) + (1+mu*A_{k})/2 * ||z_{k}-x_*||^2,
% we, instead, verify that
%   max ( Phi_{k+1} - Phi_k ) <= 0.
%
% The method originates from [1] and was adapted to deal with strong
% convexity in [2]
% [1] R. D. Monteiro and B. F. Svaiter. An accelerated hybrid proximal
%     extragradient method for convex optimization and its implications
%     to second-order methods, SIAM Journal on Optimization (2013).
%
% [2] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators (2020).


% (0) Initialize an empty PEP
P = pep();
m = 1;      % strong convexity
L = Inf;    % smoothness (no smoothness)

param.mu = m;
param.L  = L;

% (1) Set up the objective function
f = P.DeclareFunction('SmoothStronglyConvex',param);

% (2) Set up the starting point
xk      = P.StartingPoint();    % xk is a previous iterate
zk      = P.StartingPoint();    % zk is a previous iterate
[xs,fs] = f.OptimalPoint(); 	% xs is an optimal point, and fs=F(xs)
fxk     = f.value(xk);

% (3) Set up the method
sigma           = 1;
opt.criterion   = 'PD_gapI';

lambda          = 1;
Ak              = 10;
ak              = (lambda + 2* Ak*lambda*m+sqrt(4*lambda*Ak*(Ak*m+1)*(lambda*m+1)+lambda^2))/2;
Ak1             = Ak + ak;


x           = cell(2,1); x{1} = xk;
z           = cell(2,1); z{1} = zk;
fx          = cell(2,1); fx{1}= fxk;

yk          = xk + (Ak1-Ak)*(Ak*m+1)/(Ak*m*(2*Ak1-Ak)+Ak1)*(zk-xk);
[x{2},~,fx{2},w,v,~,epsVar] = inexact_proximal_step(yk,f,lambda,opt);
P.AddConstraint(epsVar <= sigma^2/2 * (yk-x{2})^2);
z{2}        = zk + (Ak1-Ak)/(Ak1*m+1)*(m*(w-zk)-v);

phi         = @(A,fx,z)(A*(fx-fs)+(1+m*A)/2*(z-xs)^2);
phik        = phi(Ak,fx{1},z{1});
phik1        = phi(Ak1,fx{2},z{2});

% (4) Set up the objective...
obj = phik1-phik;
P.PerformanceMetric(obj); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output, should be <= 0 (up to numerical precision) to
% verify the statement.Tolerance: 1e-7.
fprintf('Potential inequality verified within prescribed tolerance?[0/1]: %d\n', double(obj)<=1e-7);




