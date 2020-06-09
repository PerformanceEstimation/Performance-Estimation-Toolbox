clear all; clc;
% In this example, we use an accelerated inexact forward-backward method
% for solving the non-smooth (strongly) convex minimization problem
%   min_x {F(x) = f(x)+g(x)}; for notational convenience we denote
%   xs=argmin_x f(x)+g(x);
%
% We show how to verify the potential function of [2, Theorem 4.6]: 
% A_{k+1} (F(x_{k+1})-F_*) + (1+mu*A_{k+1})/2 * ||z_{k+1}-x_*||^2 
%             <=  A_{k} (F(x_{k})-F_*) + (1+mu*A_{k})/2 * ||z_{k}-x_*||^2
% By denoting 
%   Phi_k = A_{k} (F(x_{k})-F_*) + (1+mu*A_{k})/2 * ||z_{k}-x_*||^2,
% we, instead, verify that
%   max ( Phi_{k+1} - Phi_k ) <= 0.
%
% The exact method is from [1, Section 4.3]
%
% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators (2020).


% (0) Initialize an empty PEP
P = pep();
m = 1;      % strong convexity
L = 2;      % smoothness (no smoothness)

paramf.mu = 0;
paramf.L  = L;

paramg.mu = m;
paramg.L  = Inf;

% (1) Set up the objective function
f = P.DeclareFunction('SmoothStronglyConvex',paramf);
g = P.DeclareFunction('SmoothStronglyConvex',paramg);
F = f + g;

% (2) Set up the starting point
xk      = P.StartingPoint();    % xk is a previous iterate
zk      = P.StartingPoint();    % xk is a previous iterate
[xs,Fs] = F.OptimalPoint(); 	% xs is an optimal point, and fs=F(xs)
fxk     = f.value(xk);
gxk     = g.value(xk);

% (3) Set up the method
sigma           = .2;
zeta            = .9;
xi              = 3;

opt.criterion   = 'PD_gapI';

lambda          = (1-sigma^2)/L; % should be between 0<= lambda <= (1-sigma^2)/L
eta             = (1-zeta^2)*lambda;
Ak              = 1;
ak              = (eta + 2* Ak*eta*m+sqrt(4*eta*Ak*(Ak*m+1)*(eta*m+1)+eta^2))/2;
Ak1             = Ak + ak;


x           = cell(2,1); x{1} = xk;
z           = cell(2,1); z{1} = zk;
gx          = cell(2,1); gx{1}= gxk;
fx          = cell(2,1); fx{1}= fxk;

yk          = xk + (Ak1-Ak)*(Ak*m+1)/(Ak*m*(2*Ak1-Ak)+Ak1)*(zk-xk);
[dfyk,fyk]   = f.oracle(yk);
[x{2},~,gx{2},w,v,~,epsVar] = inexact_proximal_step(yk-lambda*dfyk,g,lambda,opt);
P.AddConstraint(epsVar <= sigma^2/2 * (yk-x{2})^2+lambda^2*zeta^2/2*(v+dfyk)^2+xi/2);
fx{2}       = f.value(x{2});
z{2}        = zk + (Ak1-Ak)/(Ak1*m+1)*(m*(w-zk)-(v+dfyk));

phi         = @(A,fx,z)(A*(fx-Fs)+(1+m*A)/2*(z-xs)^2);
phik        = phi(Ak,fx{1}+gx{1},z{1});
phik1       = phi(Ak1,fx{2}+gx{2},z{2});

% (4) Set up the objective...
obj = phik1-phik-Ak1/2/lambda*xi;
P.PerformanceMetric(obj); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output, should be <= 0 (up to numerical precision) to
% verify the statement. Tolerance: 1e-7.
fprintf('Potential inequality verified within prescribed tolerance?[0/1]: %d\n', double(obj)<=1e-7);




