clear all; clc;
% In this example, we use a partially inexact Douglas-Rachford splitting
% for solving the non-smooth (strongly) convex minimization problem
%   min_x {F(x) = f(x)+g(x)}; for notational convenience we denote
%   xs=argmin_x f(x)+g(x); where f(.) is assumed to be smooth and stronly
%   convex, and g(.) is convex.
%
% We show how to compute the convergence rate to a fixed point of the
% operator, which is denoted by z_*.
%
%
% The exact method is from [1], its PEP formulation and solution from [2].
% The precise formulation we used is described in [2, Section 4.4]
%
% [1] J. Eckstein and W. Yao, Relative-error approximate versions of
%     Douglas–Rachford splitting and special cases of the ADMM.
%     Mathematical Programming (2018).
%
% [2] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators (2020).


% (0) Initialize an empty PEP
P = pep();
m = 1;      % strong convexity
L = 5;      % smoothness (no smoothness)

paramf.mu = m;
paramf.L  = L;

% (1) Set up the objective function
f = P.DeclareFunction('SmoothStronglyConvex',paramf);
g = P.DeclareFunction('Convex');
F = f + g;

% (2) Set up the starting point
zk      = P.StartingPoint();        % zk is a previous iterate
[xs,Fs] = F.OptimalPoint('xs'); 	% xs is an optimal point, and fs=F(xs)


% (3) Set up the method
N               = 5;
sigma           = rand;
opt.criterion   = 'PD_gapII';

lambda          = rand*4;

z           = cell(2,1); z{1} = zk;

for i = 1:N
    [x,df,~,~,~,~,epsVar] = inexact_proximal_step(z{i},f,lambda,opt);
    y                     = proximal_step(x-lambda*df,g,lambda);
    P.AddConstraint(epsVar <= sigma^2/lambda^2 * (y-z{i}+lambda*df)^2);
    z{i+1}                  = z{i} + y - x;
end

% (4) Set up the objective and initial condition

zs   = xs + lambda * f.gradient('xs');
init = (z{1}-zs)^2;
obj  = (z{N+1}-zs)^2;
P.InitialCondition(init <= 1);
P.PerformanceMetric(obj); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output, and compare to what it should be, as exposed in
% [2, Theorem 4.9]
theoretical_expr = (max( ((1-sigma+lambda*m*sigma)/(1-sigma+lambda*m))^2,...
    ((sigma+(1-sigma)*lambda*L)/(1+(1-sigma)*lambda*L))^2))^N;
[double(obj) theoretical_expr]



